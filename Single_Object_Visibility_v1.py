###
# Single_Object_Visibility_v1.py
# John Phillips, john.d.phillips@comcast.net
# Creation date: 2026-10-29
###

# Finds the next observing run for one fixed celestial object using astronomical
# darkness and the configured altitude limit. Reports first, maximum, and last
# visibility, summing separate observing windows without counting their gaps.

import json
import math
import os
import sys
import warnings
from dataclasses import dataclass
from datetime import date, datetime, time, timedelta
from pathlib import Path

try:
    import astropy.units as u
    import numpy as np
    from astroplan import FixedTarget, Observer
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    from astropy.utils import iers
    from dotenv import load_dotenv
except ImportError as error:
    raise SystemExit(
        f"Missing dependency: {error.name}. Run with the project's .venv Python."
    ) from error

BASE_DIR = Path(__file__).resolve().parent
GRID_POINTS = 360
SEARCH_CYCLE_NIGHTS = 400
CROSSING_PRECISION_SECONDS = 0.1


@dataclass
class NightVisibility:
    night_date: date
    intervals: list

    # --- hours ---
    # Returns elapsed observing hours, excluding any gaps between intervals.
    @property
    def hours(self):
        return sum((end - start).to_value(u.hour) for start, end in self.intervals)


@dataclass
class VisibilityRun:
    first: NightVisibility | None
    nights: list
    underway: bool
    complete: bool
    checked: int


# --- load_configuration ---
# Reads the observing site beside this script and validates numeric settings.
# Environment values remain available when a setting is absent from .env.
def load_configuration():
    load_dotenv(BASE_DIR / '.env', override=True)
    values = {}
    for name in ('LATITUDE', 'LONGITUDE', 'ALTITUDE', 'TIMEZONE', 'ALTITUDE_LIMIT'):
        raw = os.getenv(name, '30' if name == 'ALTITUDE_LIMIT' else '')
        values[name] = raw.split('#', 1)[0].strip()
        if not values[name]:
            raise ValueError(f'Missing {name} in .env.')
    latitude, longitude, altitude, limit = (
        float(values[name]) for name in ('LATITUDE', 'LONGITUDE', 'ALTITUDE', 'ALTITUDE_LIMIT')
    )
    if not all(math.isfinite(value) for value in (latitude, longitude, altitude, limit)):
        raise ValueError('Location and altitude-limit values must be finite.')
    if not -90 <= latitude <= 90 or not -180 <= longitude <= 180:
        raise ValueError('Latitude must be -90 to 90; longitude must be -180 to 180.')
    if not 0 <= limit <= 90:
        raise ValueError('ALTITUDE_LIMIT must be between 0 and 90 degrees.')
    location = EarthLocation.from_geodetic(longitude * u.deg, latitude * u.deg, altitude * u.m)
    site_name = os.getenv('SITE_NAME', '').strip() or 'Unnamed site'
    return Observer(location=location, timezone=values['TIMEZONE'], name=site_name), limit

# --- make_target ---
# Creates a fixed ICRS target from RA in hours and declination in degrees.
def make_target(name, ra, dec):
    coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    if not name.strip() or not np.all(np.isfinite([coord.ra.deg, coord.dec.deg])):
        raise ValueError('A name and finite RA/Dec coordinates are required.')
    return FixedTarget(coord=coord, name=name.strip())

# --- select_target ---
# Prompts for a custom object, an online catalog name, or manual coordinates.
# Invalid selections can be corrected without restarting the program.
def select_target():
    while True:
        print('\nObject selection: 1. custom_objects.json  2. catalog name  3. manual RA/Dec')
        choice = input('Enter choice [1]: ').strip() or '1'
        try:
            if choice == '1':
                with (BASE_DIR / 'custom_objects.json').open(encoding='utf-8') as source:
                    objects = json.load(source)
                if not isinstance(objects, list) or not objects:
                    raise ValueError('custom_objects.json must contain a nonempty list.')
                for index, obj in enumerate(objects, 1):
                    if not isinstance(obj, dict) or not all(key in obj for key in ('name', 'ra', 'dec')):
                        raise ValueError(f'Custom object {index} needs name, ra, and dec.')
                    print(f'{index:3}. {obj["name"]}')
                selection = input('Enter object number or exact name: ').strip()
                matches = [obj for obj in objects if str(obj['name']).casefold() == selection.casefold()]
                if not matches and selection.isdigit() and 1 <= int(selection) <= len(objects):
                    matches = [objects[int(selection) - 1]]
                if len(matches) != 1:
                    raise ValueError('Select one listed object by its number or unique name.')
                obj = matches[0]
                return make_target(obj['name'], obj['ra'], obj['dec'])
            if choice == '2':
                name = input('Catalog name (for example M31): ').strip()
                if not name:
                    raise ValueError('Enter an object name.')
                return FixedTarget.from_name(name)
            if choice == '3':
                name = input('Object name: ').strip()
                ra = input('RA (for example 15h08m09.14s): ').strip()
                dec = input('Dec (for example +39d58m12.9s): ').strip()
                return make_target(name, ra, dec)
            raise ValueError('Choose 1, 2, or 3.')
        except Exception as error:
            print(f'Could not select object: {error}')

# --- valid_time ---
# Rejects absent, masked, or nonfinite scalar astronomical event times.
def valid_time(value):
    return (
        isinstance(value, Time)
        and value.isscalar
        and not bool(np.any(value.mask))
        and bool(np.isfinite(value.jd))
    )

# --- night_bounds ---
# Returns astronomical twilight boundaries for a local evening date.
# Missing twilight is an explicit error rather than a zero-visibility night.
def night_bounds(observer, night_date):
    local_noon = observer.timezone.localize(datetime.combine(night_date, time(12)), is_dst=None)
    evening = observer.twilight_evening_astronomical(
        Time(local_noon), which='next', n_grid_points=GRID_POINTS
    )
    if not valid_time(evening):
        raise ValueError(f'No evening astronomical twilight for {night_date}; this site/date is unsupported.')
    morning = observer.twilight_morning_astronomical(
        evening, which='next', n_grid_points=GRID_POINTS
    )
    if not valid_time(morning) or not evening < morning:
        raise ValueError(f'No morning astronomical twilight for night {night_date}; this site/date is unsupported.')
    if evening.to_datetime(timezone=observer.timezone).date() != night_date:
        raise ValueError(f'Evening twilight does not fall on local date {night_date}.')
    return evening, morning

# --- select_start ---
# Accepts a local evening date or uses now, including an ongoing night's remainder.
# Explicit dates include the entire observing night, even when they are in the past.
def select_start(observer, now=None):
    now = Time.now() if now is None else now
    today = now.to_datetime(timezone=observer.timezone).date()
    while True:
        entered = input(f'\nStart evening date YYYY-MM-DD [today, {today}]: ').strip()
        if entered:
            try:
                if len(entered) != 10:
                    raise ValueError
                return date.fromisoformat(entered), None
            except ValueError:
                print('Enter a date in YYYY-MM-DD format.')
                continue
        previous_date = today - timedelta(days=1)
        evening, morning = night_bounds(observer, previous_date)
        return (previous_date if evening <= now <= morning else today), now

# --- threshold_crossing ---
# Bisects one bracketed altitude crossing to a small fraction of a second.
# Exact endpoint contacts remain zero-duration events rather than disappearing.
def threshold_crossing(altitude_at, start, end, start_height, end_height, limit):
    if start_height == limit:
        return start
    if end_height == limit:
        return end
    start_above = start_height >= limit
    while (end - start).to_value(u.s) > CROSSING_PRECISION_SECONDS:
        middle = start + (end - start) / 2
        if (altitude_at(middle) >= limit) == start_above:
            start = middle
        else:
            end = middle
    return start + (end - start) / 2

# --- intervals_above_limit ---
# Splits a night at upper/lower transits, brackets crossings on each segment,
# and merges touching qualifying intervals while retaining isolated instants.
def intervals_above_limit(knots, heights, altitude_at, limit):
    intervals = []
    for start, end, start_height, end_height in zip(knots, knots[1:], heights, heights[1:]):
        start_above, end_above = start_height >= limit, end_height >= limit
        if not start_above and not end_above:
            continue
        if start_above != end_above:
            crossing = threshold_crossing(altitude_at, start, end, start_height, end_height, limit)
            if start_above:
                end = crossing
            else:
                start = crossing
        if intervals and start <= intervals[-1][1]:
            intervals[-1] = (intervals[-1][0], end)
        else:
            intervals.append((start, end))
    return intervals

# --- calculate_night ---
# Applies v3's twilight/transit altitude checks, then totals actual dark windows.
# A lower-transit split preserves separate evening and morning observing windows.
def calculate_night(observer, target, limit, night_date):
    evening, morning = night_bounds(observer, night_date)
    midpoint = evening + (morning - evening) / 2
    transit = observer.target_meridian_transit_time(
        midpoint, target, which='nearest', n_grid_points=GRID_POINTS
    )
    if not valid_time(transit) or not evening <= transit <= morning:
        next_transit = observer.target_meridian_transit_time(
            evening, target, which='next', n_grid_points=GRID_POINTS
        )
        if valid_time(next_transit) and next_transit <= morning:
            transit = next_transit
    knots = [evening, morning]
    if valid_time(transit) and evening <= transit <= morning:
        knots.append(transit)
    heights = observer.altaz(Time(knots), target).alt.to_value(u.deg)
    if not np.all(np.isfinite(heights)):
        raise ValueError(f'Could not calculate finite altitudes for {night_date}.')
    if not np.any(heights >= limit):
        return NightVisibility(night_date, [])
    lower_transit = observer.target_meridian_antitransit_time(
        midpoint, target, which='nearest', n_grid_points=GRID_POINTS
    )
    if not valid_time(lower_transit):
        raise ValueError(f'Could not calculate lower transit for {night_date}.')
    if evening < lower_transit < morning:
        knots.append(lower_transit)
    knots.sort(key=lambda moment: moment.jd)
    heights = observer.altaz(Time(knots), target).alt.to_value(u.deg)

    # --- altitude_at ---
    # Supplies altitude in degrees for refinement of a bracketed crossing.
    def altitude_at(moment):
        return observer.altaz(moment, target).alt.to_value(u.deg)

    return NightVisibility(night_date, intervals_above_limit(knots, heights, altitude_at, limit))

# --- clip_night ---
# Keeps only windows at or after a starting instant without changing full-night data.
def clip_night(night, cutoff):
    if cutoff is None:
        return night
    intervals = [(max(start, cutoff), end) for start, end in night.intervals if end >= cutoff]
    return NightVisibility(night.night_date, intervals)

# --- find_visibility_run ---
# Finds the first available run and searches through its first absent night.
# Separate 400-night safeguards bound the wait for a run and the wait for its end.
def find_visibility_run(night_at, start_date, cutoff=None, cycle_nights=SEARCH_CYCLE_NIGHTS):
    previous = night_at(start_date - timedelta(days=1))
    nights, first = [], None
    underway = False
    for offset in range(2 * cycle_nights):
        night = night_at(start_date + timedelta(days=offset))
        if first is None:
            available = clip_night(night, cutoff)
            if available.intervals:
                first = available
                underway = bool(previous.intervals) or (
                    cutoff is not None and night.intervals[0][0] < cutoff
                )
                nights.append(night)
            elif offset + 1 >= cycle_nights:
                return VisibilityRun(None, [], False, False, offset + 1)
        elif not night.intervals:
            return VisibilityRun(first, nights, underway, True, offset + 1)
        else:
            nights.append(night)
        if len(nights) >= cycle_nights:
            return VisibilityRun(first, nights, underway, False, offset + 1)
        previous = night
    raise RuntimeError('Visibility search exceeded its cycle safeguards.')

# --- maximum_night ---
# Chooses the maximum full-night duration, favoring the midpoint of tied dates.
# Two equally central dates are resolved in favor of the earlier date.
def maximum_night(nights):
    maximum = max(night.hours for night in nights)
    tied = [night for night in nights if night.hours == maximum]
    midpoint = (tied[0].night_date.toordinal() + tied[-1].night_date.toordinal()) / 2
    return min(tied, key=lambda night: (abs(night.night_date.toordinal() - midpoint), night.night_date))

# --- format_night ---
# Formats one summary row using full local timestamps and summed observing hours.
# The date's asterisk identifies an already-underway run; an instant stays explicit.
def format_night(label, night, observer, underway=False):
    first = night.intervals[0][0].to_datetime(timezone=observer.timezone)
    last = night.intervals[-1][1].to_datetime(timezone=observer.timezone)
    night_label = f'{night.night_date}{"*" if underway else ""}'
    first_text = first.strftime('%Y-%m-%d %H:%M:%S %Z')
    last_text = last.strftime('%Y-%m-%d %H:%M:%S %Z')
    duration = f'{night.hours:5.2f} h'
    if night.hours == 0:
        duration += ' (instant)'
    return (
        f'{label:<19} {night_label:<11} | '
        f'First {first_text:<25} | '
        f'Last {last_text:<25} | {duration}'
    )

# --- format_summary ---
# Produces exactly three result rows, explicitly labeling absent seasonal events.
def format_summary(result, observer, start_date):
    if result.first is None:
        return [
            f'{"First visibility:":<19} None found in {result.checked} nights from {start_date}.',
            'Maximum visibility: Unavailable; no qualifying visibility found.',
            f'{"Last visibility:":<19} Unavailable; no qualifying run found.',
        ]
    maximum = format_night('Maximum visibility:', maximum_night(result.nights), observer)
    if not result.complete:
        maximum += f' (maximum within {len(result.nights)} checked run nights)'
    return [
        format_night('First visibility:', result.first, observer, result.underway),
        maximum,
        format_night('Last visibility:', result.nights[-1], observer) if result.complete else
        f'{"Last visibility:":<19} No seasonal disappearance found in {len(result.nights)} run nights.',
    ]

# --- run_application ---
# Runs interactive selection, searches silently, and saves the three-line report.
# Expected input/configuration errors are reported without a Python traceback.
def run_application():
    try:
        observer, limit = load_configuration()
        target = select_target()
        start_date, cutoff = select_start(observer)
        print(
            f'\nObject: {target.name} | RA {target.coord.ra.to_string(unit=u.hourangle)} | '
            f'Dec {target.coord.dec.to_string(unit=u.deg, alwayssign=True)}\n'
            f'Site: {observer.name or "Unnamed site"} | '
            f'Lat {observer.location.lat.deg:.4f}, Lon {observer.location.lon.deg:.4f} | '
            f'Elevation {observer.location.height.to_value(u.m):g} m | {observer.timezone.zone}\n'
            f'Visibility criteria: Altitude >= {limit:g} deg; Sun <= -18 deg\n'
            'An asterisk marks a visibility run already underway at the start.'
        )

        # --- night_at ---
        # Calculates a single local observing night for the run-search helper.
        def night_at(night_date):
            return calculate_night(observer, target, limit, night_date)

        result = find_visibility_run(night_at, start_date, cutoff)
        summary = '\n'.join(format_summary(result, observer, start_date))
        print('\n' + summary)
        output_dir = BASE_DIR / 'output'
        output_dir.mkdir(exist_ok=True)
        timestamp = datetime.now(observer.timezone).strftime('%Y%m%d_%H%M%S_%f')
        object_name = ''.join(character if character.isalnum() else '_' for character in target.name)[:60]
        filename = output_dir / f'Single_Object_Visibility_{object_name}_{timestamp}.txt'
        try:
            filename.write_text(summary + '\n', encoding='utf-8')
            print(f'\nResults written to {filename}')
        except OSError as error:
            print(f'\nCould not save the report: {error}', file=sys.stderr)
            return 1
        return 0
    except (EOFError, KeyboardInterrupt):
        print('\nCancelled.')
        return 1
    except Exception as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1

# --- main ---
# Allows approximate seasonal planning with available Earth-rotation predictions
# when cached IERS data are old; Astropy still reports degraded-data warnings.
def main():
    with iers.conf.set_temp('auto_max_age', None):
        return run_application()


if __name__ == '__main__':
    with warnings.catch_warnings():
        # Missing astronomical events are checked explicitly before being used.
        from astroplan.exceptions import TargetAlwaysUpWarning, TargetNeverUpWarning

        warnings.simplefilter('ignore', TargetAlwaysUpWarning)
        warnings.simplefilter('ignore', TargetNeverUpWarning)
        sys.exit(main())
