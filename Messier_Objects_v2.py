# Messier_Objects.py
#
# Calculates the visibility of Messier objects and user-defined celestial objects
# during the night between evening and morning astronomical twilight.
# Includes objects that are above 30 degrees altitude at evening twilight,
# observable transit, or morning twilight, regardless of their Right Ascension.
# User-defined objects are specified in the .env file via the ADDITIONAL_OBJECTS variable.
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
    from astropy.coordinates import EarthLocation, Angle, SkyCoord # Added SkyCoord
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
    Loads location configuration and custom astronomical objects from the .env file.

    Expected variables in .env:
    - LATITUDE: Latitude in decimal degrees.
    - LONGITUDE: Longitude in decimal degrees (West is negative).
    - ALTITUDE: Altitude in meters above sea level.
    - TIMEZONE: IANA Time Zone Name (e.g., 'America/Chicago').
    - ADDITIONAL_OBJECTS (optional): JSON string representing a list of custom objects.
      Each object in the list must be a dictionary with 'name', 'ra', and 'dec' keys.
      RA and Dec should be strings parsable by astropy's SkyCoord (e.g., "HHhMMmSS.Ss", "DDdMMmSSs").
      Example: ADDITIONAL_OBJECTS=[{"name": "My Star", "ra": "10h30m00s", "dec": "+45d00m00s"}]

    Returns:
        tuple: (latitude, longitude, altitude, timezone_str, additional_objects_data)
               - latitude (float|None): Latitude in degrees.
               - longitude (float|None): Longitude in degrees.
               - altitude (float|None): Altitude in meters.
               - timezone_str (str|None): IANA timezone string.
               - additional_objects_data (list): List of dicts for custom objects, empty if none or error.
               Returns (None, None, None, None, []) if core configuration is invalid.
    """
    if not load_dotenv(override=True):
        print("Warning: .env file not found in the current directory.")

    lat_str = os.getenv('LATITUDE')
    lon_str = os.getenv('LONGITUDE')
    alt_str = os.getenv('ALTITUDE')
    timezone_str = os.getenv('TIMEZONE')
    additional_objects_str = os.getenv('ADDITIONAL_OBJECTS')

    if not all([lat_str, lon_str, alt_str, timezone_str]):
        print("Error: Missing required variables (LATITUDE, LONGITUDE, ALTITUDE, TIMEZONE) in .env file.")
        return None, None, None, None, []

    latitude, longitude, altitude = None, None, None
    try:
        latitude = float(lat_str.split('#', 1)[0].strip())
        longitude = float(lon_str.split('#', 1)[0].strip())
        altitude = float(alt_str.split('#', 1)[0].strip())
        timezone_str = timezone_str.split('#', 1)[0].strip()

        if not (-90 <= latitude <= 90):
            print(f"Error: Invalid LATITUDE: {latitude}. Must be between -90 and 90.")
            return None, None, None, None, []
        if not (-180 <= longitude <= 180):
            print(f"Error: Invalid LONGITUDE: {longitude}. Must be between -180 and 180.")
            return None, None, None, None, []

        print(f"Configuration loaded: Lat={latitude}, Lon={longitude}, Alt={altitude}m, TZ={timezone_str}")
    except ValueError as e:
        print(f"Error: Could not convert core configuration values to numbers: {e}")
        return None, None, None, None, []
    except Exception as e:
        print(f"Unexpected error during core configuration loading: {e}")
        return None, None, None, None, []

    # Process additional objects
    additional_objects_data = []
    if additional_objects_str:
        cleaned_additional_objects_str = additional_objects_str.split('#', 1)[0].strip()
        if cleaned_additional_objects_str: # Check if string is not empty after stripping comments
            try:
                parsed_objects = json.loads(cleaned_additional_objects_str)
                if not isinstance(parsed_objects, list):
                    print("Error: ADDITIONAL_OBJECTS in .env must be a JSON list. Custom objects will not be loaded.")
                else:
                    for i, obj_data in enumerate(parsed_objects):
                        if not isinstance(obj_data, dict):
                            print(f"Warning: Item {i} in ADDITIONAL_OBJECTS is not a dictionary. Skipping.")
                            continue
                        if not all(key in obj_data for key in ['name', 'ra', 'dec']):
                            print(f"Warning: Item {i} ('{obj_data.get('name', 'N/A')}') in ADDITIONAL_OBJECTS is missing 'name', 'ra', or 'dec'. Skipping.")
                            continue
                        # Basic validation passed, will be fully validated during FixedTarget creation
                        additional_objects_data.append(obj_data)
                    if additional_objects_data:
                        print(f"Loaded {len(additional_objects_data)} custom object(s) definitions from .env.")
            except json.JSONDecodeError as e:
                print(f"Error: Could not parse ADDITIONAL_OBJECTS JSON: {e}. Custom objects will not be loaded.")
            except Exception as e:
                print(f"Unexpected error loading ADDITIONAL_OBJECTS: {e}. Custom objects will not be loaded.")
        else:
            print("Info: ADDITIONAL_OBJECTS variable found in .env but is empty or fully commented out.")

    return latitude, longitude, altitude, timezone_str, additional_objects_data

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
    latitude, longitude, altitude, timezone_str, additional_objects_data = load_configuration()
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

    # 3. Calculate Twilight and LST (for output only)
    now = Time.now()
    # Fixed time for consistent testing/demonstration; remove or comment out for live runs
    now = Time("2025-06-26 01:00:00")
    
    print(f"Reference time for 'tonight': {now.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}")

    try:
        t_evening_astro_twil_end = observer.twilight_evening_astronomical(now, which='next')
        if t_evening_astro_twil_end is None or isinstance(t_evening_astro_twil_end, np.ma.core.MaskedConstant):
            raise ValueError("Could not calculate evening astronomical twilight.")
        lst_start = observer.local_sidereal_time(t_evening_astro_twil_end)

        t_morning_astro_twil_start = observer.twilight_morning_astronomical(t_evening_astro_twil_end, which='next')
        if t_morning_astro_twil_start is None or isinstance(t_morning_astro_twil_start, np.ma.core.MaskedConstant):
            raise ValueError("Could not calculate morning astronomical twilight.")
        lst_end = observer.local_sidereal_time(t_morning_astro_twil_start)

        # Calculate time of Sun's anti-meridian passage (approximate local midnight)
        # Ensure that twilight times are valid before attempting to use their JD values
        if np.isnan(t_evening_astro_twil_end.jd) or np.isnan(t_morning_astro_twil_start.jd):
             raise ValueError("Invalid Julian dates for twilight times.")
        t_sun_antimeridian_jd = (t_evening_astro_twil_end.jd + t_morning_astro_twil_start.jd) / 2.0
        if np.isnan(t_sun_antimeridian_jd): # Should not happen if inputs are not NaN
            raise ValueError("Invalid Julian date for Sun's anti-meridian crossing.")
        t_sun_antimeridian = Time(t_sun_antimeridian_jd, format='jd', scale='utc', location=observer.location)
        lst_mid = observer.local_sidereal_time(t_sun_antimeridian)

        print(f"Evening twilight ends: {t_evening_astro_twil_end.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')} (LST: {lst_start.to_string(unit=u.hourangle, sep=':', precision=0)})")
        print(f"Morning twilight starts: {t_morning_astro_twil_start.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')} (LST: {lst_end.to_string(unit=u.hourangle, sep=':', precision=0)})")
        print(f"Sun's anti-meridian: {t_sun_antimeridian.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')} (LST: {lst_mid.to_string(unit=u.hourangle, sep=':', precision=0)})")

    except ValueError as e:
        print(f"Error calculating twilight times: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error during twilight calculation: {e}")
        sys.exit(1)

    altitude_limit = 30 * u.deg
    print(f"Altitude threshold: Objects must reach at least {altitude_limit}")

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

    # Add custom objects from .env
    if additional_objects_data:
        print(f"Preparing {len(additional_objects_data)} custom object(s) from .env...")
        for custom_obj_info in additional_objects_data:
            name = custom_obj_info['name']
            ra_str = custom_obj_info['ra']
            dec_str = custom_obj_info['dec']
            try:
                # Create SkyCoord from RA/Dec strings. Assumes ICRS frame (common for J2000).
                coord = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg), frame='icrs')
                target = FixedTarget(coord=coord, name=name)
                all_targets_to_process.append({"name": name, "target_obj": target, "source": "Custom"})
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
                t_sun_antimeridian, target, which="nearest"
            )
            rise_time = observer.target_rise_time(
                transit_time, target, which="previous", horizon=altitude_limit
            )
            set_time = observer.target_set_time(
                transit_time, target, which="next", horizon=altitude_limit
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
                "rise_time": rise_time if not is_circumpolar_above_limit else None,
                "transit_time": transit_time,
                "set_time": set_time if not is_circumpolar_above_limit else None,
                "max_altitude": max_alt_at_transit.to_value(u.deg),
                "is_circumpolar_above_limit": is_circumpolar_above_limit, # Renamed for clarity
                "rise_outside": rise_outside,
                "set_outside": set_outside,
                "source": item["source"]
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

    print("\n" + "\n".join(output_lines))

    # Write to timestamped file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename = f"Celestial_Objects_Visibility_{timestamp}.txt" # Renamed for generality
    try:
        with open(output_filename, 'w') as f:
            f.write(f"# Celestial Object Visibility Report\n")
            f.write(f"# Location: Lat={latitude:.4f}, Lon={longitude:.4f}, Alt={altitude:.0f}m\n")
            f.write(f"# Timezone: {observer.timezone.zone}\n")
            f.write(f"# Observation Window (Astronomical Twilight): {t_evening_astro_twil_end.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')} to {t_morning_astro_twil_start.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')}\n")
            f.write(f"# LST Range: {lst_start.to_string(unit=u.hourangle, sep=':', precision=0)} to {lst_end.to_string(unit=u.hourangle, sep=':', precision=0)}\n")
            f.write(f"# Sun Anti-Meridian (Approx. Local Midnight): {t_sun_antimeridian.to_datetime(timezone=observer.timezone).strftime('%Y-%m-%d %H:%M:%S %Z')} (LST: {lst_mid.to_string(unit=u.hourangle, sep=':', precision=0)})\n")
            f.write(f"# Total Observable Objects in Report: {len(visible_objects)}\n")
            f.write(f"# Altitude Threshold: > {altitude_limit}\n")
            f.write("# Objects included if: above 30° at evening twilight, OR transit is observable and above 30°, OR above 30° at morning twilight.\n")
            f.write("# Includes standard Messier objects and any custom objects defined in .env (ADDITIONAL_OBJECTS).\n")
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