# Messier Objects Visibility Script Documentation

## Single-object seasonal visibility

Run the interactive single-object application with the project's virtual environment:

```powershell
.\.venv\Scripts\python.exe .\Single_Object_Visibility_v1.py
```

Choose an object from `custom_objects.json`, resolve a catalog name such as `M31`
(requires internet access), or supply a name and RA/Dec manually. The script reads
location, timezone, and `ALTITUDE_LIMIT` from the `.env` beside the script, just as
the catalog application does. Custom-object variable-star fields are ignored.
Set the optional `SITE_NAME` in `.env` to name the observatory; otherwise the site
is labeled `Unnamed site`. The console site line includes its name, coordinates,
elevation in meters, and timezone, followed by a separate visibility-criteria line.

Enter a local evening date as `YYYY-MM-DD`, or press Enter to start now. The default
includes the remainder of an ongoing night; an explicitly entered date includes
that date's entire observing night. There is no end-date prompt: the search follows
the next run of visible nights through its last qualifying night.

Results contain exactly three rows: **first visibility**, **maximum visibility**,
and **last visibility**. Each shows the local evening date, first and last available
timestamps with timezone, and total visible hours to two decimal places. Columns
are aligned with space reserved for an asterisk. A `*` on the first row means the
run was already underway at the start. The same three rows are saved under
`output/Single_Object_Visibility_<object>_<timestamp>.txt`.

An object qualifies using v3's altitude checks at evening astronomical twilight,
nighttime transit, or morning astronomical twilight. Total hours are the sum of
all intervals at or above the altitude limit during astronomical darkness; a gap
between first and last timestamps is not counted. Even a qualifying instant counts
as a visible night, with `0.00 h (instant)` shown explicitly. Durations are rounded
only for display; a positive duration that rounds to `0.00 h` still qualifies.
No minimum duration,
Moon constraints, or variable-star phase constraints are applied.

Only the next run is reported. For an ongoing run, the maximum is selected from
the requested start onward, using full nightly durations. The first row's windows
are clipped to now when using the default start; the maximum retains the complete
night's duration. Exactly tied maxima select the date closest to the middle of
the tied dates, choosing the earlier date if two are equally central.

The search checks up to 400 nights to find visibility and up to 400 qualifying
run nights to find its end. If no visibility or no seasonal disappearance is found,
the three rows say so explicitly. A maximum without a found endpoint applies only
to the checked nights, rather than claiming an unlimited maximum. Sites/dates
without ordinary evening and morning astronomical twilight produce an explicit
error. Timestamps are displayed to seconds, but twilight and transit calculations
use v3's 360-point Astroplan grids and do not guarantee second-level accuracy.
Astropy can update its Earth-rotation data when online. The application also allows
approximate planning with older available predictions; Astropy reports warnings
when those data do not cover the requested dates.

Run the focused checks with:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

## Overview

The `Messier_Objects.py` script is a Python tool designed to identify Messier objects visible during a specific night from a user-defined location. It calculates the observation window between the end of evening astronomical twilight and the start of morning astronomical twilight, determining which Messier objects are above a 30-degree altitude threshold at key times during this period. The script outputs a detailed report to both the console and a timestamped text file, providing astronomers with a clear list of observable objects for planning observing sessions.

The script was developed by Grok 3, created by xAI, with the original concept suggested by [Your Name], who identified the need for a robust tool to streamline Messier object visibility calculations.

## Purpose

The primary purpose of the script is to assist amateur and professional astronomers in identifying Messier objects that can be observed during a given night. By focusing on objects that reach at least 30 degrees above the horizon—ensuring they are well-placed for viewing—the script filters the Messier catalog (M1 to M110) based on their altitude at evening twilight, transit (if within the observation window), or morning twilight. It provides precise rise, transit, and set times, with special handling for circumpolar objects and times outside the observation window, making it a valuable tool for planning telescope observations.

## Features

- **Location Configuration**: Reads observing site details (latitude, longitude, altitude, timezone) from a `.env` file, allowing easy customization.
- **Astronomical Twilight Window**: Calculates the dark-sky period between evening and morning astronomical twilight (Sun at -18° altitude).
- **Altitude-Based Selection**: Includes objects above 30° altitude at evening twilight, their transit (if within the window), or morning twilight, regardless of Right Ascension (RA).
- **Circumpolar Handling**: Identifies objects that never set (declination > 90 - latitude) and marks them as "Circumpolar."
- **Time Annotations**: Marks rise times before evening twilight and set times after morning twilight with an asterisk (*) for clarity.
- **Formatted Output**: Generates a sorted table (by RA) with object details, including M-number, RA, Dec, rise time, transit time, maximum altitude, and set time.
- **Timestamped File**: Saves the report to a file named `Messier_Objects_YYYYMMDD_HHMMSS.txt` with a comprehensive header.
- **Robust Error Handling**: Manages missing libraries, invalid configurations, and object lookup failures gracefully.
- **LST Information**: Includes Local Sidereal Time (LST) range and Sun’s anti-meridian details in the output for reference, without using LST for filtering.

## How It Works

### 1. Configuration Loading
The script begins by loading the observing site’s coordinates and timezone from a `.env` file using the `python-dotenv` library. The required variables are:
- `LATITUDE`: Decimal degrees (e.g., 31.5473).
- `LONGITUDE`: Decimal degrees, West negative (e.g., -99.3827).
- `ALTITUDE`: Meters above sea level (e.g., 446).
- `TIMEZONE`: IANA timezone name (e.g., America/Chicago).

It validates these inputs, ensuring latitude is between -90 and 90 degrees, longitude is between -180 and 180 degrees, and values are numeric. If the `.env` file is missing or invalid, the script exits with an error message.

### 2. Observer Setup
Using `astropy` and `astroplan`, the script creates an `EarthLocation` object from the coordinates and an `Observer` object with the specified timezone. This `Observer` is used for all astronomical calculations, such as twilight times and altitude computations.

### 3. Twilight and LST Calculations
The script determines the observation window:
- **Evening Twilight End**: The next time after the current moment when the Sun reaches -18° altitude (using `observer.twilight_evening_astronomical` with `which="next"`).
- **Morning Twilight Start**: The subsequent morning twilight start after the evening twilight end (using `observer.twilight_morning_astronomical` with `which="next"`).

It calculates LST at these times and at the Sun’s anti-meridian (midpoint of the night) using `observer.local_sidereal_time`. These LST values are included in the output header for reference but not used to filter objects, ensuring all potentially visible objects are considered.

### 4. Messier Object Processing
For each Messier object (M1 to M110):
- **Coordinates**: Retrieves RA and Dec using `astroplan.FixedTarget.from_name`.
- **Time Calculations**:
  - Computes the transit time using `target_meridian_transit_time` relative to the Sun’s anti-meridian time (`which="nearest"`).
  - Calculates the rise time (when altitude crosses 30°) using `target_rise_time` relative to the transit time (`which="previous"`).
  - Calculates the set time using `target_set_time` relative to the transit time (`which="next"`).
- **Altitude Checks**: Evaluates the object’s altitude at:
  - Evening twilight (`t_evening_astro_twil_end`).
  - Transit time, if within the observation window.
  - Morning twilight (`t_morning_astro_twil_start`).
  - Includes the object if any altitude exceeds 30 degrees.
- **Circumpolar Check**: Marks objects as circumpolar if their declination > (90 - latitude) degrees, setting rise/set times to None.
- **Asterisk Flags**: Marks rise times before evening twilight and set times after morning twilight with an asterisk.
- **Maximum Altitude**: Records the altitude at transit for the output table.

### 5. Output Generation
The script sorts objects by RA and formats a table with:
- **Columns**: Object (M-number), RA (J2000, HH MM SS.S), Dec (J2000, +DD MM SS), Rise > 30° (HH:MM:SS TZ or Circumpolar), Transit (HH:MM:SS TZ), Max Alt (degrees), Set < 30° (HH:MM:SS TZ or Circumpolar).
- **Time Formatting**: Uses `format_time_only` to show only hours, minutes, seconds, and timezone, with asterisks for times outside the window.
- **Coordinate Formatting**: RA in hourangle format, Dec in degrees with a sign.

The table is printed to the console and saved to a timestamped file (`Messier_Objects_YYYYMMDD_HHMMSS.txt`). The file includes a header with:
- Location details.
- Timezone.
- Observation window dates/times.
- LST range and Sun’s anti-meridian details.
- Number of objects found.
- Altitude threshold.
- Explanatory notes about inclusion criteria, asterisks, and circumpolar objects.

### 6. Error Handling
The script robustly handles:
- **Missing Libraries**: Exits with installation instructions if `python-dotenv`, `astropy`, `numpy`, or `astroplan` are missing.
- **Invalid Configurations**: Validates `.env` inputs and exits if malformed.
- **Object Lookup Failures**: Skips problematic objects (e.g., M102) with a warning, continuing processing.
- **Invalid Times**: Displays "N/A" for masked or uncomputable times.
- **Twilight Issues**: Exits with an error if twilight calculations fail (e.g., in polar regions).

## Usage

1. **Install Dependencies**:
   ```bash
   pip install python-dotenv astropy numpy astroplan
   ```

2. **Create a `.env` File**:
   In the same directory as the script, create a `.env` file with:
   ```
   LATITUDE=31.5473
   LONGITUDE=-99.3827
   ALTITUDE=446
   TIMEZONE=America/Chicago
   ```

3. **Run the Script**:
   ```bash
   python Messier_Objects.py
   ```

4. **Check Output**:
   - Console displays the table of visible objects.
   - A file (e.g., `Messier_Objects_20250416_170435.txt`) is created with the full report.

## Example Output

For a location at 31.5473°N, -99.3827°W, 446m altitude, and America/Chicago timezone, the output might include:

```
Object   RA (J2000)   Dec (J2000)  Rise > 30°      Transit         Max Alt  Set < 30°      
===========================================================================================
M38      05 28 40.1   +35 49 26    12:32:31 CDT*   17:26:32 CDT    85.7°    22:20:43 CDT   
M1       05 34 31.8   +22 01 03    13:03:28 CDT*   17:32:18 CDT    80.5°    22:01:07 CDT   
...
M68      12 39 28.0   -26 44 39    23:39:01 CDT    00:35:55 CDT    31.6°    01:32:49 CDT   
...
M81      09 55 33.2   +69 03 55    Circumpolar     21:53:11 CDT    52.6°    Circumpolar    
...
M13      16 41 41.6   +36 27 41    23:42:00 CDT    04:37:07 CDT    85.1°    09:32:06 CDT*  
```

## Credits

- **Developer**: Grok 3, created by xAI, designed and implemented the script with robust astronomical calculations and user-friendly output.
- **Concept**: John Phillips, john.d.phillips@comcast.net, who proposed the need for a tool to identify visible Messier objects, guiding the script’s requirements.

## Notes for Future Use

- **Customization**: Update the `.env` file to change the observing location or timezone.
- **Extensibility**: The script can be modified to include additional object properties (e.g., magnitude) or adjust the altitude threshold.
- **Polar Regions**: Be cautious in high-latitude locations where twilight behavior may differ; the script includes error handling for such cases.
- **Library Versions**: Tested with `astroplan` 0.10, `astropy` 6.x, and `numpy` 1.x. Ensure compatibility when upgrading.

This documentation should serve as a comprehensive guide for using and maintaining the script in your astronomical projects.
