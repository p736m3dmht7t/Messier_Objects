# Messier_Objects.py
#
# Calculates the visibility of Messier objects during the night between evening
# and morning astronomical twilight. Includes objects that are above 30 degrees
# altitude at evening twilight, observable transit, or morning twilight,
# regardless of their Right Ascension.
#
# Uses astroplan's direct rise/set/transit methods. An object is circumpolar if
# its declination > (90 - latitude) degrees. Rise/set times outside the observation
# window are marked with an asterisk (*). Results are printed to the console and
# saved to a timestamped text file. LST range is included in the output for reference.

import os
import sys
import warnings
from datetime import datetime

# Import third-party libraries
try:
    from dotenv import load_dotenv
except ImportError:
    print("Error: python-dotenv library not found. Please install it: pip install python-dotenv")
    sys.exit(1)

try:
    import numpy as np
    import astropy.units as u
    from astropy.time import Time
    from astropy.coordinates import EarthLocation, Angle
    from astropy.utils.exceptions import AstropyWarning
    warnings.simplefilter('ignore', category=AstropyWarning)
except ImportError:
    print("Error: astropy or numpy library not found. Please install them: pip install astropy numpy")
    sys.exit(1)

try:
    from astroplan import Observer, FixedTarget
except ImportError:
    print("Error: astroplan library not found. Please install it: pip install astroplan")
    sys.exit(1)

# --- Configuration Loading ---

def load_configuration():
    """
    Loads location configuration from the .env file in the current directory.

    Expected variables in .env:
    - LATITUDE: Latitude in decimal degrees.
    - LONGITUDE: Longitude in decimal degrees (West is negative).
    - ALTITUDE: Altitude in meters above sea level.
    - TIMEZONE: IANA Time Zone Name (e.g., 'America/Chicago').

    Returns:
        tuple: (latitude, longitude, altitude, timezone_str) or (None, None, None, None)
               if configuration is invalid or missing.
    """
    if not load_dotenv(override=True):
        print("Warning: .env file not found in the current directory.")

    lat_str = os.getenv('LATITUDE')
    lon_str = os.getenv('LONGITUDE')
    alt_str = os.getenv('ALTITUDE')
    timezone_str = os.getenv('TIMEZONE')

    if not all([lat_str, lon_str, alt_str, timezone_str]):
        print("Error: Missing required variables (LATITUDE, LONGITUDE, ALTITUDE, TIMEZONE) in .env file.")
        return None, None, None, None

    try:
        latitude = float(lat_str.split('#', 1)[0].strip())
        longitude = float(lon_str.split('#', 1)[0].strip())
        altitude = float(alt_str.split('#', 1)[0].strip())
        timezone_str = timezone_str.split('#', 1)[0].strip()

        if not (-90 <= latitude <= 90):
            print(f"Error: Invalid LATITUDE: {latitude}. Must be between -90 and 90.")
            return None, None, None, None
        if not (-180 <= longitude <= 180):
            print(f"Error: Invalid LONGITUDE: {longitude}. Must be between -180 and 180.")
            return None, None, None, None

        print(f"Configuration loaded: Lat={latitude}, Lon={longitude}, Alt={altitude}m, TZ={timezone_str}")
        return latitude, longitude, altitude, timezone_str
    except ValueError as e:
        print(f"Error: Could not convert configuration values to numbers: {e}")
        return None, None, None, None
    except Exception as e:
        print(f"Unexpected error during configuration loading: {e}")
        return None, None, None, None

# --- Formatting Helpers ---

def format_ra(ra_angle):
    """Formats an Astropy Angle (RA) to 'HH MM SS.S' string."""
    return ra_angle.to_string(unit=u.hourangle, sep=' ', precision=1, pad=True)

def format_dec(dec_angle):
    """Formats an Astropy Angle (Dec) to '+DD MM SS' string."""
    return dec_angle.to_string(unit=u.degree, sep=' ', precision=0, pad=True, alwayssign=True)

def format_time_only(astro_time, observer, is_circumpolar=False, is_outside=False):
    """
    Formats an Astropy Time to 'HH:MM:SS TZ', omitting the date. Adds an asterisk (*) if the time is outside the observation window.

    Args:
        astro_time (astropy.time.Time or None): The time to format.
        observer (astroplan.Observer): Observer with timezone info.
        is_circumpolar (bool): If True, returns "Circumpolar".
        is_outside (bool): If True and not circumpolar, appends an asterisk.

    Returns:
        str: Formatted time string, "Circumpolar", or "N/A" if invalid.
    """
    if is_circumpolar:
        return "Circumpolar"
    if astro_time is None or isinstance(astro_time, np.ma.core.MaskedConstant):
        return "N/A"
    try:
        local_dt = astro_time.to_datetime(timezone=observer.timezone)
        time_str = local_dt.strftime('%H:%M:%S %Z')
        return f"{time_str}*" if is_outside else time_str
    except Exception as e:
        print(f"Warning: Could not format time {astro_time}: {e}")
        return "N/A"

# --- Main Calculation Logic ---

def main():
    """
    Main function to calculate Messier object visibility.
    - Loads configuration and sets up observer.
    - Calculates observation window (twilight times) and LST range for reference.
    - Includes Messier objects above 30° at evening twilight, observable transit, or morning twilight.
    - Computes rise, transit, set times, and max altitude using astroplan.
    - Marks rise/set times outside the observation window with an asterisk (*).
    - Outputs results to console and a timestamped file.
    """
    print("Starting Messier Object Visibility Calculation...")

    # 1. Load Configuration
    latitude, longitude, altitude, timezone_str = load_configuration()
    if latitude is None:
        print("Exiting due to configuration errors.")
        sys.exit(1)

    # 2. Set up Observer
    try:
        location = EarthLocation.from_geodetic(
            lon=longitude * u.deg,
            lat=latitude * u.deg,
            height=altitude * u.m
        )
        observer = Observer(location=location, timezone=timezone_str)
        print(f"Observer created for {location.geodetic} in timezone {observer.timezone.zone}")
    except Exception as e:
        print(f"Error creating Observer: {e}")
        print("Ensure LATITUDE, LONGITUDE, ALTITUDE are valid and TIMEZONE is a valid IANA name.")
        sys.exit(1)

    # 3. Calculate Twilight and LST (for output only)
    now = Time.now()
#    now = Time("2025-06-26 01:00:00")
    
    print(f"Current time: {now.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}")

    try:
        t_evening_astro_twil_end = observer.twilight_evening_astronomical(now, which='next')
        if t_evening_astro_twil_end is None or isinstance(t_evening_astro_twil_end, np.ma.core.MaskedConstant):
            raise ValueError("Could not calculate evening astronomical twilight.")
        print(f"Evening twilight ends: {t_evening_astro_twil_end.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}")
        lst_start = observer.local_sidereal_time(t_evening_astro_twil_end)
        print(f"LST at evening twilight: {lst_start.to_string(unit=u.hourangle, sep=':', precision=0)}")

        t_morning_astro_twil_start = observer.twilight_morning_astronomical(t_evening_astro_twil_end, which='next')
        if t_morning_astro_twil_start is None or isinstance(t_morning_astro_twil_start, np.ma.core.MaskedConstant):
            raise ValueError("Could not calculate morning astronomical twilight.")
        print(f"Morning twilight starts: {t_morning_astro_twil_start.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}")
        lst_end = observer.local_sidereal_time(t_morning_astro_twil_start)
        print(f"LST at morning twilight: {lst_end.to_string(unit=u.hourangle, sep=':', precision=0)}")

        t_sun_antimeridian_jd = (t_evening_astro_twil_end.jd + t_morning_astro_twil_start.jd) / 2.0
        if np.isnan(t_sun_antimeridian_jd):
            raise ValueError("Invalid Julian date for Sun's anti-meridian crossing.")
        t_sun_antimeridian = Time(t_sun_antimeridian_jd, format='jd', scale='utc', location=observer.location)
        print(f"Sun's anti-meridian: {t_sun_antimeridian.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}")
        lst_mid = observer.local_sidereal_time(t_sun_antimeridian)
        print(f"LST at Sun anti-meridian: {lst_mid.to_string(unit=u.hourangle, sep=':', precision=0)}")

    except ValueError as e:
        print(f"Error calculating twilight times: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error during twilight calculation: {e}")
        sys.exit(1)

    altitude_limit = 30 * u.deg
    print(f"Altitude threshold: Objects must reach at least {altitude_limit}")

    # 4. Process Messier Objects
    messier_list = [f"M{i}" for i in range(1, 111)]
    visible_objects = []
    print(f"\nChecking {len(messier_list)} Messier objects...")

    for m_name in messier_list:
        try:
            target = FixedTarget.from_name(m_name)
            ra = target.coord.ra
            dec = target.coord.dec

            # Calculate transit, rise, set times
            transit_time = observer.target_meridian_transit_time(
                t_sun_antimeridian, target, which="nearest"
            )
            rise_time = observer.target_rise_time(
                transit_time, target, which="previous", horizon=altitude_limit
            )
            set_time = observer.target_set_time(
                transit_time, target, which="next", horizon=altitude_limit
            )

            # Check altitudes at key times
            altitudes = []
            # Evening twilight
            alt_evening = observer.altaz(t_evening_astro_twil_end, target).alt
            altitudes.append(alt_evening)
            # Transit (if within window)
            alt_transit = 0 * u.deg
            if transit_time is not None and not isinstance(transit_time, np.ma.core.MaskedConstant):
                if t_evening_astro_twil_end <= transit_time <= t_morning_astro_twil_start:
                    alt_transit = observer.altaz(transit_time, target).alt
            altitudes.append(alt_transit)
            # Morning twilight
            alt_morning = observer.altaz(t_morning_astro_twil_start, target).alt
            altitudes.append(alt_morning)

            # Include object if above 30° at any key time
            if not any(alt.to_value(u.deg) > altitude_limit.to_value(u.deg) for alt in altitudes):
                continue

            # Check for circumpolar objects (declination > 90 - latitude)
            is_circumpolar = dec > (90 - latitude) * u.deg

            # Determine if rise/set times are outside the observation window
            rise_outside = rise_time is not None and rise_time < t_evening_astro_twil_end
            set_outside = set_time is not None and set_time > t_morning_astro_twil_start

            # Calculate max altitude at transit (for output)
            max_altitude = alt_transit if alt_transit > 0 * u.deg else observer.altaz(transit_time, target).alt

            # Store results
            visible_objects.append({
                "m_number": m_name,
                "ra": ra,
                "dec": dec,
                "rise_time": rise_time if not is_circumpolar else None,
                "transit_time": transit_time,
                "set_time": set_time if not is_circumpolar else None,
                "max_altitude": max_altitude.to_value(u.deg),
                "is_circumpolar": is_circumpolar,
                "rise_outside": rise_outside,
                "set_outside": set_outside
            })

        except Exception as e:
            print(f"Warning: Could not process {m_name}: {e}")
            continue

    print(f"\nFound {len(visible_objects)} Messier objects meeting the criteria.")

    # 5. Format and Output Results
    if not visible_objects:
        print("No suitable Messier objects found for the upcoming night.")
        return

    # Sort by RA for readability
    visible_objects.sort(key=lambda x: x['ra'].hour)

    col_widths = {'obj': 8, 'ra': 12, 'dec': 12, 'rise': 15, 'transit': 15, 'alt': 8, 'set': 15}
    header = (
        f"{'Object':<{col_widths['obj']}} {'RA (J2000)':<{col_widths['ra']}} {'Dec (J2000)':<{col_widths['dec']}} "
        f"{'Rise > 30°':<{col_widths['rise']}} {'Transit':<{col_widths['transit']}} {'Max Alt':<{col_widths['alt']}} {'Set < 30°':<{col_widths['set']}}"
    )
    separator = "=" * len(header)
    output_lines = [header, separator]

    for obj in visible_objects:
        rise_str = format_time_only(
            obj['rise_time'], observer, obj['is_circumpolar'], obj['rise_outside']
        )
        set_str = format_time_only(
            obj['set_time'], observer, obj['is_circumpolar'], obj['set_outside']
        )
        transit_str = format_time_only(obj['transit_time'], observer)
        max_alt_str = f"{obj['max_altitude']:.1f}°"

        line = (
            f"{obj['m_number']:<{col_widths['obj']}} {format_ra(obj['ra']):<{col_widths['ra']}} {format_dec(obj['dec']):<{col_widths['dec']}} "
            f"{rise_str:<{col_widths['rise']}} {transit_str:<{col_widths['transit']}} {max_alt_str:<{col_widths['alt']}} {set_str:<{col_widths['set']}}"
        )
        output_lines.append(line)

    print("\n" + "\n".join(output_lines))

    # Write to timestamped file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename = f"Messier_Objects_{timestamp}.txt"
    try:
        with open(output_filename, 'w') as f:
            f.write(f"# Messier Object Visibility Report\n")
            f.write(f"# Location: Lat={latitude:.4f}, Lon={longitude:.4f}, Alt={altitude:.0f}m\n")
            f.write(f"# Timezone: {observer.timezone.zone}\n")
            f.write(f"# Observation Window: {t_evening_astro_twil_end.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')} to {t_morning_astro_twil_start.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}\n")
            f.write(f"# LST Range: {lst_start.to_string(unit=u.hourangle, sep=':', precision=0)} to {lst_end.to_string(unit=u.hourangle, sep=':', precision=0)}\n")
            f.write(f"# Sun Anti-Meridian: {t_sun_antimeridian.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}\n")
            f.write(f"# LST at Sun Anti-Meridian: {lst_mid.to_string(unit=u.hourangle, sep=':', precision=0)}\n")
            f.write(f"# Total Observable Objects: {len(visible_objects)}\n")
            f.write(f"# Altitude Threshold: > {altitude_limit}\n")
            f.write("# Note: Objects included if above 30° at evening twilight, observable transit, or morning twilight.\n")
            f.write("# Note: Times are based on direct astroplan calculations.\n")
            f.write("# Note: Times marked with '*' occur before evening twilight (rise) or after morning twilight (set).\n")
            f.write("# Note: Circumpolar objects have declination > (90 - latitude) degrees.\n")
            f.write("#\n")
            f.write("\n".join(output_lines))
        print(f"\nResults written to {output_filename}")
    except IOError as e:
        print(f"\nError writing to {output_filename}: {e}")

    print("\nScript finished.")

if __name__ == "__main__":
    main()