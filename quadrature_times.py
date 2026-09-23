#!/usr/bin/env python3
"""
Propose spectroscopic observation times at quadrature for an eclipsing binary.

Quadratures:
  Q1 = photometric phase 0.25
  Q2 = photometric phase 0.75

T0 is the epoch of primary minimum (phase 0.0) in the system given by
--time-sys (default BJD_TDB). Period P is in days in that same system.
Each quadrature instant T = T0 + (n + phase)*P is converted to UTC at
the observatory (topocentric light-travel time for BJD/HJD).

Start date is local civil noon at the observatory (timezone from
longitude). The end date is inclusive.

A time is kept only if:
  - the star is at least --min-alt degrees above the horizon, and
  - the Sun is below --sun-alt degrees (default -18, astronomical night).

Lunar phase is the illuminated fraction of the Moon (0 = new, 1 = full).
Lunar separation is the sky angle between the Moon and the target.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, GCRS, SkyCoord, get_body, get_sun
from astropy.time import Time


TIME_SYS_CHOICES = ("bjd_tdb", "hjd", "jd_utc")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Quadrature times for spectroscopic observations of an EB"
    )
    p.add_argument("--name", default="", help="System name printed in the header")
    p.add_argument("--lat", type=float, required=True, help="Observatory latitude (deg, +N)")
    p.add_argument("--lon", type=float, required=True, help="Observatory longitude (deg, +E)")
    p.add_argument("--elev", type=float, default=0.0, help="Observatory elevation (m)")
    p.add_argument("--ra", required=True, help="Star RA (e.g. 12:34:56.7 or 188.736 deg)")
    p.add_argument("--dec", required=True, help="Star Dec (e.g. +12:34:56 or 12.582 deg)")
    p.add_argument("--t0", type=float, required=True, help="Epoch of primary minimum")
    p.add_argument("--period", type=float, required=True, dest="p", help="Orbital period (days)")
    p.add_argument(
        "--time-sys",
        choices=TIME_SYS_CHOICES,
        default="bjd_tdb",
        help="Time system of T0 (and P): bjd_tdb (default), hjd (HJD_UTC), jd_utc",
    )
    p.add_argument("--start", required=True, help="Start date YYYY-MM-DD (local noon)")
    p.add_argument("--end", required=True, help="End date YYYY-MM-DD")
    p.add_argument("--min-alt", type=float, default=30.0, help="Minimum star altitude (deg)")
    p.add_argument(
        "--sun-alt",
        type=float,
        default=-18.0,
        help="Sun must be below this altitude (deg); default -18",
    )
    p.add_argument(
        "--ra-unit",
        choices=["auto", "hour", "deg"],
        default="auto",
        help="How to parse --ra (auto: ':' or 'h' => hourangle)",
    )
    return p.parse_args()


def parse_coord(ra_str: str, dec_str: str, ra_unit: str) -> SkyCoord:
    ra_looks_sexagesimal = (":" in ra_str) or ("h" in ra_str.lower())
    if ra_unit == "hour" or (ra_unit == "auto" and ra_looks_sexagesimal):
        return SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg))
    try:
        ra_val = float(ra_str)
        return SkyCoord(ra=ra_val * u.deg, dec=dec_str, unit=u.deg)
    except ValueError:
        return SkyCoord(ra=ra_str, dec=dec_str, unit=(u.deg, u.deg))


def local_noon_utc(date_str: str, lon_deg: float) -> Time:
    """Civil 12:00 at approximate local mean time from longitude -> UTC Time."""
    dt_naive = datetime.strptime(date_str, "%Y-%m-%d").replace(hour=12, minute=0, second=0)
    offset_hours = lon_deg / 15.0
    dt_utc = dt_naive - timedelta(hours=offset_hours)
    dt_utc = dt_utc.replace(tzinfo=timezone.utc)
    return Time(dt_utc)


def _ltt_kind(time_sys: str) -> str | None:
    if time_sys == "bjd_tdb":
        return "barycentric"
    if time_sys == "hjd":
        return "heliocentric"
    return None


def native_from_utc(
    t_utc: Time, star: SkyCoord, location: EarthLocation, time_sys: str
) -> np.ndarray:
    """UTC Time -> JD values in the ephemeris time system."""
    kind = _ltt_kind(time_sys)
    if kind is None:
        return np.atleast_1d(t_utc.utc.jd).astype(float)

    if time_sys == "bjd_tdb":
        t = t_utc.tdb
        native_jd = t.jd + t.light_travel_time(star, kind=kind, location=location).to(u.day).value
        return np.atleast_1d(native_jd).astype(float)

    # HJD_UTC: heliocentric correction applied to UTC JD
    t = t_utc.utc
    native_jd = t.jd + t.light_travel_time(star, kind=kind, location=location).to(u.day).value
    return np.atleast_1d(native_jd).astype(float)


def utc_from_native(
    native_jd: np.ndarray, star: SkyCoord, location: EarthLocation, time_sys: str
) -> Time:
    """Ephemeris JD values -> topocentric UTC Time (iterated light-time)."""
    native_jd = np.atleast_1d(np.asarray(native_jd, dtype=float))
    kind = _ltt_kind(time_sys)

    if kind is None:
        return Time(native_jd, format="jd", scale="utc")

    if time_sys == "bjd_tdb":
        t = Time(native_jd, format="jd", scale="tdb")
        for _ in range(2):
            ltt = t.light_travel_time(star, kind=kind, location=location)
            t = Time(native_jd, format="jd", scale="tdb") - ltt
        return t.utc

    # HJD_UTC
    t = Time(native_jd, format="jd", scale="utc")
    for _ in range(2):
        ltt = t.light_travel_time(star, kind=kind, location=location)
        t = Time(native_jd, format="jd", scale="utc") - ltt
    return t.utc


def quadrature_native(
    t0: float, period: float, native_start: float, native_end: float
) -> list[tuple[float, str]]:
    """All Q1 / Q2 times in the native time system with start <= T < end."""
    out: list[tuple[float, str]] = []
    for phase, label in ((0.25, "Q1"), (0.75, "Q2")):
        n_min = np.ceil((native_start - t0) / period - phase)
        n_max = np.floor((native_end - t0) / period - phase - 1e-12)
        if n_max < n_min:
            continue
        ns = np.arange(int(n_min), int(n_max) + 1)
        ts = t0 + (ns + phase) * period
        for t in ts:
            if native_start <= t < native_end:
                out.append((float(t), label))
    out.sort(key=lambda x: x[0])
    return out


def format_ra(star: SkyCoord) -> tuple[str, str]:
    deg = f"{star.ra.deg:10.5f}"
    hms = star.ra.to_string(unit=u.hourangle, sep=":", precision=2, pad=True)
    return deg, hms


def format_dec(star: SkyCoord) -> tuple[str, str]:
    deg = f"{star.dec.deg:+10.5f}"
    dms = star.dec.to_string(unit=u.deg, sep=":", precision=1, alwayssign=True, pad=True)
    return deg, dms


def moon_illumination(sun: SkyCoord, moon: SkyCoord) -> np.ndarray:
    """Illuminated fraction of the Moon (0=new, 1=full)."""
    elongation = sun.separation(moon).rad
    phase_angle = np.pi - elongation
    return 0.5 * (1.0 + np.cos(phase_angle))


def time_sys_label(time_sys: str) -> str:
    return {"bjd_tdb": "BJD_TDB", "hjd": "HJD_UTC", "jd_utc": "JD_UTC"}[time_sys]


def main() -> None:
    args = parse_args()
    location = EarthLocation(lat=args.lat * u.deg, lon=args.lon * u.deg, height=args.elev * u.m)
    star = parse_coord(args.ra, args.dec, args.ra_unit)

    t_start = local_noon_utc(args.start, args.lon)
    t_end = local_noon_utc(args.end, args.lon) + 1.0 * u.day

    native_start = float(native_from_utc(t_start, star, location, args.time_sys)[0])
    native_end = float(native_from_utc(t_end, star, location, args.time_sys)[0])

    events = quadrature_native(args.t0, args.p, native_start, native_end)
    name = args.name.strip() or "Unnamed system"

    ra_deg, ra_hms = format_ra(star)
    dec_deg, dec_dms = format_dec(star)

    print(name)
    print(f"RA   {ra_deg} deg    {ra_hms}")
    print(f"Dec  {dec_deg} deg   {dec_dms}")
    print(f"T0={args.t0:.6f} {time_sys_label(args.time_sys)}   P={args.p:.6f} d")
    print(
        f"Site lat={args.lat:.5f}  lon={args.lon:.5f}  elev={args.elev:.0f} m"
        f"   {args.start} → {args.end}"
    )
    print()

    if not events:
        print("No quadrature times in the requested window.")
        return

    native_jds = np.array([e[0] for e in events])
    labels = [e[1] for e in events]
    times = utc_from_native(native_jds, star, location, args.time_sys)

    altaz = AltAz(obstime=times, location=location)
    star_alt = star.transform_to(altaz).alt.deg

    sun = get_sun(times)
    sun_alt = sun.transform_to(altaz).alt.deg

    moon = get_body("moon", times)
    star_gcrs = star.transform_to(GCRS(obstime=times))
    moon_sep = star_gcrs.separation(moon).deg
    illum = moon_illumination(sun, moon)

    keep = (star_alt >= args.min_alt) & (sun_alt <= args.sun_alt)

    hdr = f"{'Q':<4} {'UTC':^22} {'JD':^14} {'Elev':>6} {'Moon%':>6} {'MoonSep':>8}"
    print(hdr)
    print("-" * len(hdr))

    n_keep = 0
    for lab, t, alt, frac, sep, ok in zip(labels, times, star_alt, illum, moon_sep, keep):
        if not ok:
            continue
        n_keep += 1
        utc = t.strftime("%Y-%m-%d %H:%M:%S")
        print(f"{lab:<4} {utc:<22} {t.utc.jd:14.6f} {alt:6.1f} {100.0 * frac:6.0f} {sep:8.1f}")


    print()
    print(
        f"{n_keep} of {len(events)} quadratures with "
        f"star ≥ {args.min_alt:.0f}° and Sun ≤ {args.sun_alt:.0f}°"
    )


if __name__ == "__main__":
    main()
