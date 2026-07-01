# Messier_Objects_v3.py
#
# Calculates the visibility of Messier objects and user-defined celestial objects
# during the night between evening and morning astronomical twilight.
# Includes objects that are above 30 degrees altitude at evening twilight,
# observable transit, or morning twilight, regardless of their Right Ascension.
# User-defined objects are specified in the custom_objects.json file.
# Custom variable stars can include T0 in HJD_UTC and P in days.
#
# Uses astroplan's direct rise/set/transit methods. An object is circumpolar if
# its declination > (90 - latitude) degrees. Rise/set times outside the observation
# window are marked with an asterisk (*). Results are printed to the console and
# saved to a timestamped text file. LST range is included in the output for reference.
#
# Sections:
# - load_configuration
# - Formatting Helpers (format_ra, format_dec, format_time_only)
# - main

import os
import sys
import warnings
from datetime import datetime
import json # Added for parsing custom objects

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
    from astropy.coordinates import EarthLocation, Angle, SkyCoord, get_body, get_sun, GeocentricTrueEcliptic
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


# --- load_configuration ---
# Loads location, timezone, and custom astronomical objects from .env.
# See docstring for .env variable details.
def load_configuration():
    """
    Loads location configuration from the .env file and custom astronomical objects from custom_objects.json.

    Expected variables in .env:
    - LATITUDE: Latitude in decimal degrees.
    - LONGITUDE: Longitude in decimal degrees (West is negative).
    - ALTITUDE: Altitude in meters above sea level.
    - TIMEZONE: IANA Time Zone Name (e.g., 'America/Chicago').
    - ALTITUDE_LIMIT: Minimum altitude in degrees for visibility (optional, defaults to 30).
    
    Expected format of custom_objects.json (optional):
    JSON array of objects, each with 'name', 'ra', and 'dec' keys.
    RA and Dec should be strings parsable by astropy's SkyCoord (e.g., "HHhMMmSS.Ss", "DDdMMmSSs").
    Variable stars can include 'T0' as HJD UTC and 'P' as period in days.
    Example: [{"name": "My Star", "ra": "10h30m00s", "dec": "+45d00m00s"}]

    Returns:
        tuple: (latitude, longitude, altitude, timezone_str, altitude_limit, additional_objects_data)
               - latitude (float|None): Latitude in degrees.
               - longitude (float|None): Longitude in degrees.
               - altitude (float|None): Altitude in meters.
               - timezone_str (str|None): IANA timezone string.
               - altitude_limit (float|None): Altitude limit in degrees.
               - additional_objects_data (list): List of dicts for custom objects, empty if none or error.
               Returns (None, None, None, None, None, []) if core configuration is invalid.
    """
    if not load_dotenv(override=True):
        print("Warning: .env file not found in the current directory.")

    lat_str = os.getenv('LATITUDE')
    lon_str = os.getenv('LONGITUDE')
    alt_str = os.getenv('ALTITUDE')
    timezone_str = os.getenv('TIMEZONE')
    altitude_limit_str = os.getenv('ALTITUDE_LIMIT')
    if not all([lat_str, lon_str, alt_str, timezone_str]):
        print("Error: Missing required variables (LATITUDE, LONGITUDE, ALTITUDE, TIMEZONE) in .env file.")
        return None, None, None, None, None, []

    latitude, longitude, altitude, altitude_limit = None, None, None, None
    try:
        latitude = float(lat_str.split('#', 1)[0].strip())
        longitude = float(lon_str.split('#', 1)[0].strip())
        altitude = float(alt_str.split('#', 1)[0].strip())
        timezone_str = timezone_str.split('#', 1)[0].strip()

        # Parse ALTITUDE_LIMIT with default fallback
        if altitude_limit_str:
            altitude_limit = float(altitude_limit_str.split('#', 1)[0].strip())
        else:
            altitude_limit = 30.0  # Default value
            print("Warning: ALTITUDE_LIMIT not found in .env file. Using default value of 30 degrees.")

        if not (-90 <= latitude <= 90):
            print(f"Error: Invalid LATITUDE: {latitude}. Must be between -90 and 90.")
            return None, None, None, None, None, []
        if not (-180 <= longitude <= 180):
            print(f"Error: Invalid LONGITUDE: {longitude}. Must be between -180 and 180.")
            return None, None, None, None, None, []
        if not (0 <= altitude_limit <= 90):
            print(f"Error: Invalid ALTITUDE_LIMIT: {altitude_limit}. Must be between 0 and 90.")
            return None, None, None, None, None, []

        print(f"Configuration loaded: Lat={latitude}, Lon={longitude}, Alt={altitude}m, TZ={timezone_str}, Altitude Limit={altitude_limit}°")
    except ValueError as e:
        print(f"Error: Could not convert core configuration values to numbers: {e}")
        return None, None, None, None, None, []
    except Exception as e:
        print(f"Unexpected error during core configuration loading: {e}")
        return None, None, None, None, None, []

    # Load custom objects from JSON file
    additional_objects_data = []
    custom_objects_file = 'custom_objects.json'
    if os.path.exists(custom_objects_file):
        try:
            with open(custom_objects_file, 'r') as f:
                parsed_objects = json.load(f)
                if not isinstance(parsed_objects, list):
                    print("Error: custom_objects.json must contain a JSON list. Custom objects will not be loaded.")
                else:
                    for i, obj_data in enumerate(parsed_objects):
                        if not isinstance(obj_data, dict):
                            print(f"Warning: Item {i} in custom_objects.json is not a dictionary. Skipping.")
                            continue
                        if not all(key in obj_data for key in ['name', 'ra', 'dec']):
                            print(f"Warning: Item {i} ('{obj_data.get('name', 'N/A')}') in custom_objects.json is missing 'name', 'ra', or 'dec'. Skipping.")
                            continue
                        # Basic validation passed, will be fully validated during FixedTarget creation
                        additional_objects_data.append(obj_data)
                    if additional_objects_data:
                        print(f"Loaded {len(additional_objects_data)} custom object(s) definitions from custom_objects.json.")
        except json.JSONDecodeError as e:
            print(f"Error: Could not parse custom_objects.json: {e}. Custom objects will not be loaded.")
        except Exception as e:
            print(f"Unexpected error loading custom_objects.json: {e}. Custom objects will not be loaded.")
    else:
        print("Info: custom_objects.json file not found. No custom objects will be loaded.")

    return latitude, longitude, altitude, timezone_str, altitude_limit, additional_objects_data

# --- Formatting Helpers ---
# Utility functions for formatting astronomical coordinates and time for display.

# --- format_ra ---
# Formats an Astropy Angle (RA) to 'HH MM SS.S' string.
def format_ra(ra_angle):
    """Formats an Astropy Angle (RA) to 'HH MM SS.S' string."""
    return ra_angle.to_string(unit=u.hourangle, sep=' ', precision=1, pad=True)

# --- format_dec ---
# Formats an Astropy Angle (Dec) to '+DD MM SS' string.
def format_dec(dec_angle):
    """Formats an Astropy Angle (Dec) to '+DD MM SS' string."""
    return dec_angle.to_string(unit=u.degree, sep=' ', precision=0, pad=True, alwayssign=True)

# --- format_lst ---
# Formats an Astropy Angle (LST) to 'HH:MM:SS' string with leading zeros.
def format_lst(lst_angle):
    """Formats an Astropy Angle (LST) to 'HH:MM:SS' string with leading zeros."""
    return lst_angle.to_string(unit=u.hourangle, sep=':', precision=0, pad=True)

# --- is_valid_time ---
# Checks whether an Astropy time value can be used for event-based calculations.
def is_valid_time(astro_time):
    """Returns True when astro_time is available and not masked."""
    return astro_time is not None and not isinstance(astro_time, np.ma.core.MaskedConstant)

# --- get_moon_phase ---
# Calculates the Moon phase name and illuminated percentage for the report header.
def get_moon_phase(time):
    """Returns the Moon phase name and illumination percentage at the given time."""
    sun = get_sun(time).transform_to(GeocentricTrueEcliptic(obstime=time))
    moon = get_body('moon', time).transform_to(GeocentricTrueEcliptic(obstime=time))
    phase_angle = (moon.lon - sun.lon).wrap_at(360 * u.deg).to_value(u.deg)
    illuminated_fraction = (1 - np.cos(np.deg2rad(phase_angle))) / 2

    if phase_angle < 22.5 or phase_angle >= 337.5:
        phase_name = "New Moon"
    elif phase_angle < 67.5:
        phase_name = "Waxing Crescent"
    elif phase_angle < 112.5:
        phase_name = "First Quarter"
    elif phase_angle < 157.5:
        phase_name = "Waxing Gibbous"
    elif phase_angle < 202.5:
        phase_name = "Full Moon"
    elif phase_angle < 247.5:
        phase_name = "Waning Gibbous"
    elif phase_angle < 292.5:
        phase_name = "Last Quarter"
    else:
        phase_name = "Waning Crescent"

    return phase_name, illuminated_fraction * 100

# --- calculate_variable_phase ---
# Calculates variable-star phase using event time converted to HJD UTC.
def calculate_variable_phase(astro_time, coord, observer, t0_hjd_utc, period_days):
    """Returns the phase at astro_time, or None when the time is unavailable."""
    if not is_valid_time(astro_time):
        return None

    light_travel_time = astro_time.light_travel_time(
        coord, kind='heliocentric', location=observer.location
    )
    event_hjd_utc = (astro_time.utc + light_travel_time).jd
    return ((event_hjd_utc - t0_hjd_utc) / period_days) % 1

# --- calculate_moon_separation ---
# Calculates the angular separation between a target and the Moon at a given time.
def calculate_moon_separation(astro_time, coord, observer):
    """Returns Moon separation in degrees, or None when the time is unavailable."""
    if not is_valid_time(astro_time):
        return None

    from astropy.coordinates import AltAz

    altaz_frame = AltAz(obstime=astro_time, location=observer.location)
    target_altaz = coord.transform_to(altaz_frame)
    moon_coord = get_body('moon', astro_time, observer.location)
    moon_altaz = moon_coord.transform_to(altaz_frame)
    return target_altaz.separation(moon_altaz).to_value(u.deg)

# --- format_optional_value ---
# Formats optional phase and Moon separation values for variable-star rows.
def format_optional_value(value, precision, suffix=""):
    """Formats a numeric value or returns N/A when the value is unavailable."""
    if value is None:
        return "N/A"
    return f"{value:.{precision}f}{suffix}"

# --- calculate_horizon_ra ---
# Calculates the Right Ascension of objects at the altitude limit for specific azimuths
# Uses the defined altitude_limit (30°) instead of true horizon (0°)
def calculate_horizon_ra(observer, time, azimuth, altitude):
    """
    Calculate the RA of an object at the specified altitude at a specific azimuth.
    
    Args:
        observer (astroplan.Observer): Observer with location and timezone info.
        time (astropy.time.Time): Time for the calculation.
        azimuth (float): Azimuth in degrees (90° = East, 270° = West).
        altitude (astropy.units.Quantity): Altitude (e.g., 30°).
    
    Returns:
        astropy.coordinates.Angle: Right Ascension at the specified altitude.
    """
    from astropy.coordinates import AltAz, ICRS
    
    # Create horizontal coordinate at specified altitude
    horizon_coord = AltAz(alt=altitude, az=azimuth*u.deg, 
                         obstime=time, location=observer.location)
    
    # Convert to equatorial coordinates (RA/Dec)
    equatorial_coord = horizon_coord.transform_to(ICRS())
    
    return equatorial_coord.ra

# --- format_time_only ---
# Formats an Astropy Time to 'HH:MM:SS TZ', omitting the date, with optional markers.
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

# --- later_valid_time ---
# Selects the later valid time, falling back to the available valid time.
def later_valid_time(first_time, second_time):
    """Returns the later available Astropy time, or None if both are unavailable."""
    first_valid = is_valid_time(first_time)
    second_valid = is_valid_time(second_time)
    if first_valid and second_valid:
        return first_time if first_time >= second_time else second_time
    if first_valid:
        return first_time
    if second_valid:
        return second_time
    return None

# --- earlier_valid_time ---
# Selects the earlier valid time, falling back to the available valid time.
def earlier_valid_time(first_time, second_time):
    """Returns the earlier available Astropy time, or None if both are unavailable."""
    first_valid = is_valid_time(first_time)
    second_valid = is_valid_time(second_time)
    if first_valid and second_valid:
        return first_time if first_time <= second_time else second_time
    if first_valid:
        return first_time
    if second_valid:
        return second_time
    return None

# --- main ---
# Main function to calculate Messier and custom object visibility.
# Handles observer setup, twilight calculation, object processing, and output generation.
def main():
    """
    Main function to calculate Messier and custom object visibility.
    - Loads configuration and sets up observer.
    - Calculates observation window (twilight times) and LST range for reference.
    - Includes objects above 30° at evening twilight, observable transit, or morning twilight.
    - Computes rise, transit, set times, and max altitude using astroplan.
    - Marks rise/set times outside the observation window with an asterisk (*).
    - Outputs results to console and a timestamped file.
    """
    print("Starting Celestial Object Visibility Calculation...")

    # 1. Load Configuration
    latitude, longitude, altitude, timezone_str, altitude_limit, additional_objects_data = load_configuration()
    if latitude is None: # Core configuration failed
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

    # 3. Get observation date from user
    print("\nDate Selection:")
    print("1. Tonight (default)")
    print("2. Specific date")
    
    choice = input("Enter choice (1 or 2, press Enter for tonight): ").strip()
    
    if choice == "2":
        while True:
            date_input = input("Enter date (YYYY-MM-DD) or datetime (YYYY-MM-DD HH:MM:SS): ").strip()
            try:
                # Try to parse as datetime first, then as date
                if len(date_input) > 10:  # Likely includes time
                    now = Time(date_input)
                else:  # Just date, use noon as reference
                    now = Time(f"{date_input} 12:00:00")
                break
            except Exception as e:
                print(f"Invalid date format. Please use YYYY-MM-DD or YYYY-MM-DD HH:MM:SS")
                continue
    else:
        now = Time.now()
    
    print(f"Reference time for observation: {now.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}")

    try:
        t_evening_astro_twil_end = observer.twilight_evening_astronomical(now, which='next', n_grid_points=360)
        if t_evening_astro_twil_end is None or isinstance(t_evening_astro_twil_end, np.ma.core.MaskedConstant):
            raise ValueError("Could not calculate evening astronomical twilight.")
        lst_start = observer.local_sidereal_time(t_evening_astro_twil_end)

        t_morning_astro_twil_start = observer.twilight_morning_astronomical(t_evening_astro_twil_end, which='next', n_grid_points=360)
        if t_morning_astro_twil_start is None or isinstance(t_morning_astro_twil_start, np.ma.core.MaskedConstant):
            raise ValueError("Could not calculate morning astronomical twilight.")
        lst_end = observer.local_sidereal_time(t_morning_astro_twil_start)

        # Calculate time of Sun's anti-meridian passage (approximate local midnight)
        t_sun_antimeridian_jd = (t_evening_astro_twil_end.jd + t_morning_astro_twil_start.jd) / 2.0
        t_sun_antimeridian = Time(t_sun_antimeridian_jd, format='jd', scale='utc', location=observer.location)
        lst_mid = observer.local_sidereal_time(t_sun_antimeridian)
        moon_phase_name, moon_illumination_percent = get_moon_phase(t_sun_antimeridian)
        moon_at_antimeridian = get_body('moon', t_sun_antimeridian)

        print(f"Evening Twilight Ends:    {t_evening_astro_twil_end.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}  (LST: {format_lst(lst_start)})")
        print(f"Sun Anti-Meridian:        {t_sun_antimeridian.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}  (LST: {format_lst(lst_mid)})")
        print(f"Morning Twilight Begins:  {t_morning_astro_twil_start.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}  (LST: {format_lst(lst_end)})")

    except ValueError as e:
        print(f"Error calculating twilight times: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error during twilight calculation: {e}")
        sys.exit(1)

    altitude_limit = altitude_limit * u.deg
    print(f"Altitude threshold: Objects must reach at least {altitude_limit}")

    # Calculate horizon RAs at altitude limit for key times
    ra_east_evening = calculate_horizon_ra(observer, t_evening_astro_twil_end, 90, altitude_limit)
    ra_west_evening = calculate_horizon_ra(observer, t_evening_astro_twil_end, 270, altitude_limit)
    ra_east_midnight = calculate_horizon_ra(observer, t_sun_antimeridian, 90, altitude_limit)
    ra_west_midnight = calculate_horizon_ra(observer, t_sun_antimeridian, 270, altitude_limit)
    ra_east_morning = calculate_horizon_ra(observer, t_morning_astro_twil_start, 90, altitude_limit)
    ra_west_morning = calculate_horizon_ra(observer, t_morning_astro_twil_start, 270, altitude_limit)

    print(f"Horizon RAs at {altitude_limit} altitude:")
    print(f"  Evening:  East RA {format_ra(ra_east_evening)}, West RA {format_ra(ra_west_evening)}")
    print(f"  Midnight: East RA {format_ra(ra_east_midnight)}, West RA {format_ra(ra_west_midnight)}")
    print(f"  Morning:  East RA {format_ra(ra_east_morning)}, West RA {format_ra(ra_west_morning)}")

    # 4. Prepare list of targets (Messier + Custom)
    all_targets_to_process = []

    # Add Messier objects
    messier_names = [f"M{i}" for i in range(1, 111)]
    print(f"\nPreparing {len(messier_names)} Messier objects...")
    for m_name in messier_names:
        try:
            target = FixedTarget.from_name(m_name)
            all_targets_to_process.append({"name": m_name, "target_obj": target, "source": "Messier"})
        except Exception as e:
            # This can happen if object name resolution fails (e.g. network issue, or bad name)
            print(f"Warning: Could not resolve Messier object {m_name}: {e}. Skipping.")

    # Add custom objects from custom_objects.json
    if additional_objects_data:
        print(f"Preparing {len(additional_objects_data)} custom object(s) from custom_objects.json...")
        for custom_obj_info in additional_objects_data:
            name = custom_obj_info['name']
            ra_str = custom_obj_info['ra']
            dec_str = custom_obj_info['dec']
            try:
                # Create SkyCoord from RA/Dec strings. Assumes ICRS frame (common for J2000).
                coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame='icrs')
                target = FixedTarget(coord=coord, name=name)
                variable_info = None
                if 'T0' in custom_obj_info and 'P' in custom_obj_info:
                    t0_hjd_utc = float(custom_obj_info['T0'])
                    period_days = float(custom_obj_info['P'])
                    if period_days <= 0:
                        raise ValueError("P must be greater than zero.")
                    variable_info = {"T0": t0_hjd_utc, "P": period_days}
                all_targets_to_process.append({
                    "name": name,
                    "target_obj": target,
                    "source": "Custom",
                    "variable_info": variable_info
                })
            except ValueError as e: # Catches errors from SkyCoord parsing (e.g., bad format)
                print(f"Warning: Could not create target for custom object '{name}' (RA='{ra_str}', Dec='{dec_str}'): {e}. Skipping.")
            except Exception as e: # Other unexpected errors
                print(f"Warning: Unexpected error processing custom object '{name}': {e}. Skipping.")
    
    if not all_targets_to_process:
        print("No objects (Messier or custom) could be prepared for processing. Exiting.")
        sys.exit(1)

    print(f"\nChecking {len(all_targets_to_process)} total objects for visibility...")
    visible_objects = []

    for item in all_targets_to_process:
        target_name = item["name"]
        target = item["target_obj"]
        
        try:
            ra = target.coord.ra
            dec = target.coord.dec

            transit_time = observer.target_meridian_transit_time(
                t_sun_antimeridian, target, which="nearest", n_grid_points=360
            )
            rise_time = observer.target_rise_time(
                transit_time, target, which="previous", horizon=altitude_limit, n_grid_points=360
            )
            set_time = observer.target_set_time(
                transit_time, target, which="next", horizon=altitude_limit, n_grid_points=360
            )

            altitudes_at_key_times = []
            alt_evening = observer.altaz(t_evening_astro_twil_end, target).alt
            altitudes_at_key_times.append(alt_evening)

            alt_transit_in_window = 0 * u.deg # Default to 0 if not transiting in window
            if transit_time is not None and not isinstance(transit_time, np.ma.core.MaskedConstant):
                if t_evening_astro_twil_end <= transit_time <= t_morning_astro_twil_start:
                    alt_transit_in_window = observer.altaz(transit_time, target).alt
            altitudes_at_key_times.append(alt_transit_in_window)
            
            alt_morning = observer.altaz(t_morning_astro_twil_start, target).alt
            altitudes_at_key_times.append(alt_morning)

            if not any(alt.to_value(u.deg) >= altitude_limit.to_value(u.deg) for alt in altitudes_at_key_times):
                continue # Skip if not above limit at any key time

            is_circumpolar_above_horizon = dec > (90 * u.deg - observer.location.lat)
            is_circumpolar_above_limit = False # Specific check for circumpolar above 30 deg
            if is_circumpolar_above_horizon:
                # For circumpolar objects, check if minimum altitude > altitude_limit
                # Min altitude = dec - (90 - lat) for upper culmination, or (90-lat) - dec for lower (if it passes north)
                # Simpler: check altitude at any point. If always above limit, it's effectively "Circumpolar" for rise/set.
                # A robust check for "always above 30 deg" requires more detailed analysis.
                # The current astroplan rise/set with horizon=30 should handle this by returning masked times.
                # If rise/set times are masked for a non-circumpolar (never reaches 30 deg), it's filtered out by altitude check.
                # If rise/set times are masked for a circumpolar (always above 30 deg), then it's effectively circumpolar *above limit*.
                if isinstance(rise_time, np.ma.core.MaskedConstant) and \
                   isinstance(set_time, np.ma.core.MaskedConstant) and \
                   observer.altaz(transit_time, target).alt > altitude_limit: # Ensure it transits high enough
                    is_circumpolar_above_limit = True


            rise_outside = False
            if rise_time is not None and not isinstance(rise_time, np.ma.core.MaskedConstant):
                 rise_outside = rise_time < t_evening_astro_twil_end

            set_outside = False
            if set_time is not None and not isinstance(set_time, np.ma.core.MaskedConstant):
                set_outside = set_time > t_morning_astro_twil_start
            
            max_alt_at_transit = 0 * u.deg
            if transit_time is not None and not isinstance(transit_time, np.ma.core.MaskedConstant):
                max_alt_at_transit = observer.altaz(transit_time, target).alt
            else: # Should not happen if target_meridian_transit_time is robust
                max_alt_at_transit = observer.target_upper_transit_altitude(target)


            visible_objects.append({
                "name": target_name,
                "ra": ra,
                "dec": dec,
                "target_coord": target.coord,
                "rise_time": rise_time if not is_circumpolar_above_limit else None,
                "transit_time": transit_time,
                "set_time": set_time if not is_circumpolar_above_limit else None,
                "max_altitude": max_alt_at_transit.to_value(u.deg),
                "is_circumpolar_above_limit": is_circumpolar_above_limit, # Renamed for clarity
                "rise_outside": rise_outside,
                "set_outside": set_outside,
                "source": item["source"],
                "variable_info": item.get("variable_info")
            })

        except Exception as e:
            print(f"Warning: Could not process {target_name}: {e}. Skipping.")
            continue

    print(f"\nFound {len(visible_objects)} objects meeting the visibility criteria.")

    # 5. Format and Output Results
    if not visible_objects:
        print("No suitable objects found for the upcoming night.")
        return

    visible_objects.sort(key=lambda x: x['ra'].hour)

    col_widths = {'obj': 20, 'ra': 12, 'dec': 12, 'rise': 15, 'transit': 15, 'alt': 8, 'set': 15}
    header = (
        f"{'Object':<{col_widths['obj']}} {'RA (J2000)':<{col_widths['ra']}} {'Dec (J2000)':<{col_widths['dec']}} "
        f"{'Rise > 30°':<{col_widths['rise']}} {'Transit':<{col_widths['transit']}} {'Max Alt':<{col_widths['alt']}} {'Set < 30°':<{col_widths['set']}}"
    )
    separator = "=" * len(header)
    output_lines = [header, separator]

    for obj in visible_objects:
        rise_str = format_time_only(
            obj['rise_time'], observer, obj['is_circumpolar_above_limit'], obj['rise_outside']
        )
        set_str = format_time_only(
            obj['set_time'], observer, obj['is_circumpolar_above_limit'], obj['set_outside']
        )
        transit_str = format_time_only(obj['transit_time'], observer) # No outside marker for transit
        max_alt_str = f"{obj['max_altitude']:.1f}°"

        line = (
            f"{obj['name']:<{col_widths['obj']}} {format_ra(obj['ra']):<{col_widths['ra']}} {format_dec(obj['dec']):<{col_widths['dec']}} "
            f"{rise_str:<{col_widths['rise']}} {transit_str:<{col_widths['transit']}} {max_alt_str:<{col_widths['alt']}} {set_str:<{col_widths['set']}}"
        )
        output_lines.append(line)

    variable_objects = [obj for obj in visible_objects if obj.get("variable_info")]
    if variable_objects:
        detail_label_width = col_widths['obj'] + col_widths['ra'] + col_widths['dec'] + 3
        variable_header = (
            f"{'Object':<{col_widths['obj']}} {'RA (J2000)':<{col_widths['ra']}} {'Dec (J2000)':<{col_widths['dec']}} "
            f"{'First':<{col_widths['rise']}} {'Transit':<{col_widths['transit']}} {'Max Alt':<{col_widths['alt']}} {'Last':<{col_widths['set']}}"
        )
        variable_separator = "=" * len(variable_header)
        output_lines.extend(["", "Variable Stars", "=" * len("Variable Stars"), variable_header, variable_separator])
        for obj in variable_objects:
            first_time = later_valid_time(obj['rise_time'], t_evening_astro_twil_end)
            last_time = earlier_valid_time(obj['set_time'], t_morning_astro_twil_start)
            first_str = format_time_only(first_time, observer)
            last_str = format_time_only(last_time, observer)
            transit_str = format_time_only(obj['transit_time'], observer)
            max_alt_str = f"{obj['max_altitude']:.1f}°"

            line = (
                f"{obj['name']:<{col_widths['obj']}} {format_ra(obj['ra']):<{col_widths['ra']}} {format_dec(obj['dec']):<{col_widths['dec']}} "
                f"{first_str:<{col_widths['rise']}} {transit_str:<{col_widths['transit']}} {max_alt_str:<{col_widths['alt']}} {last_str:<{col_widths['set']}}"
            )
            output_lines.append(line)

            t0_hjd_utc = obj["variable_info"]["T0"]
            period_days = obj["variable_info"]["P"]
            phase_first = calculate_variable_phase(
                first_time, obj['target_coord'], observer, t0_hjd_utc, period_days
            )
            phase_transit = calculate_variable_phase(
                obj['transit_time'], obj['target_coord'], observer, t0_hjd_utc, period_days
            )
            phase_last = calculate_variable_phase(
                last_time, obj['target_coord'], observer, t0_hjd_utc, period_days
            )
            moon_sep_first = calculate_moon_separation(first_time, obj['target_coord'], observer)
            moon_sep_transit = calculate_moon_separation(obj['transit_time'], obj['target_coord'], observer)
            moon_sep_last = calculate_moon_separation(last_time, obj['target_coord'], observer)
            phase_first_str = format_optional_value(phase_first, 2)
            phase_transit_str = format_optional_value(phase_transit, 2)
            phase_last_str = format_optional_value(phase_last, 2)
            moon_sep_first_str = format_optional_value(moon_sep_first, 0, ' deg')
            moon_sep_transit_str = format_optional_value(moon_sep_transit, 0, ' deg')
            moon_sep_last_str = format_optional_value(moon_sep_last, 0, ' deg')

            output_lines.append(
                f"{'  Phase':<{detail_label_width}}"
                f"{phase_first_str:<{col_widths['rise']}} "
                f"{phase_transit_str:<{col_widths['transit']}} "
                f"{'':<{col_widths['alt']}} "
                f"{phase_last_str:<{col_widths['set']}}"
            )
            output_lines.append(
                f"{'  Moon Sep':<{detail_label_width}}"
                f"{moon_sep_first_str:<{col_widths['rise']}} "
                f"{moon_sep_transit_str:<{col_widths['transit']}} "
                f"{'':<{col_widths['alt']}} "
                f"{moon_sep_last_str:<{col_widths['set']}}"
            )

    print("\n" + "\n".join(output_lines))

    # Write to timestamped file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    output_filename = os.path.join(output_dir, f"Celestial_Objects_Visibility_{timestamp}.txt")
    try:
        with open(output_filename, 'w') as f:
            f.write(f"# Celestial Object Visibility Report\n")
            f.write(f"# Location: Lat={latitude:.4f}, Lon={longitude:.4f}, Alt={altitude:.0f}m\n")
            f.write(f"# Timezone: {observer.timezone.zone}\n")
            f.write(f"# Moon Phase: {moon_phase_name} ({moon_illumination_percent:.1f}% illuminated)\n")
            f.write(f"# Moon at Sun Anti-Meridian: RA {format_ra(moon_at_antimeridian.ra)}  Dec {format_dec(moon_at_antimeridian.dec)}\n")
            f.write(f"# Evening Twilight Ends:    {t_evening_astro_twil_end.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}  (LST: {format_lst(lst_start)})  East RA: {format_ra(ra_east_evening)}  West RA: {format_ra(ra_west_evening)}\n")
            f.write(f"# Sun Anti-Meridian:        {t_sun_antimeridian.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}  (LST: {format_lst(lst_mid)})  East RA: {format_ra(ra_east_midnight)}  West RA: {format_ra(ra_west_midnight)}\n")
            f.write(f"# Morning Twilight Begins:  {t_morning_astro_twil_start.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}  (LST: {format_lst(lst_end)})  East RA: {format_ra(ra_east_morning)}  West RA: {format_ra(ra_west_morning)}\n")
            f.write(f"# Total Observable Objects in Report: {len(visible_objects)}\n")
            f.write(f"# Altitude Threshold: > {altitude_limit}\n")
            f.write("# Objects included if: above 30° at evening twilight, OR transit is observable and above 30°, OR above 30° at morning twilight.\n")
            f.write("# Includes standard Messier objects and any custom objects defined in custom_objects.json.\n")
            f.write("# Times are based on direct astroplan calculations for the specified horizon.\n")
            f.write("# Times marked with '*' occur before evening twilight (rise) or after morning twilight (set).\n")
            f.write(f"# 'Circumpolar' indicates object remains above {altitude_limit.value}° for the entire night period.\n")
            f.write("#\n")
            f.write("\n".join(output_lines))
        print(f"\nResults written to {output_filename}")
    except IOError as e:
        print(f"\nError writing to {output_filename}: {e}")

    print("\nScript finished.")

if __name__ == "__main__":
    main()
