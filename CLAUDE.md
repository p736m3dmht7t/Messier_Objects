# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Workflow Instructions

1. First think through the problem, read the codebase for relevant files, and write a plan to tasks/todo.md.
2. The plan should have a list of todo items that you can check off as you complete them
3. Before you begin working, check in with me and I will verify the plan.
4. Then, begin working on the todo items, marking them as complete as you go.
5. Please every step of the way just give me a high level explanation of what changes you made
6. Make every task and code change you do as simple as possible. We want to avoid making any massive or complex changes. Every change should impact as little code as possible. Everything is about simplicity.
7. Finally, add a review section to the [todo.md](http://todo.md/) file with a summary of the changes you made and any other relevant information.
8. Follow my preferred commenting style in ./comment_style.md

## Running the Script

```bash
# Activate the virtual environment first
source .venv/bin/activate

# Run the main script
python Messier_Objects_v2.py
```

## Environment Setup

The script requires a `.env` file in the project root (not committed to git):
```
LATITUDE=31.5473
LONGITUDE=-99.3827
ALTITUDE=446
TIMEZONE=America/Chicago
ALTITUDE_LIMIT=30   # optional, defaults to 30
```

Dependencies are in `.venv/` (Python 3.13). Key packages: `astropy`, `astroplan`, `numpy`, `python-dotenv`.

## Architecture

There are two script versions:
- `Messier_Objects.py` — original, Messier catalog only (M1–M110)
- `Messier_Objects_v2.py` — current version, adds support for user-defined objects from `custom_objects.json`

Both scripts follow the same flow:
1. `load_configuration()` — reads `.env` and (v2 only) `custom_objects.json`
2. Build `astroplan.Observer` from the location
3. Compute evening/morning astronomical twilight times (Sun at -18°) for the current night
4. For each object, compute transit, rise (above altitude limit), and set times using `astroplan` methods
5. Filter to objects above the altitude threshold at evening twilight, transit, or morning twilight
6. Sort by RA and print a formatted table to console + timestamped `.txt` file

**Key design decisions:**
- Circumpolar threshold: declination > (90 − latitude) degrees
- Rise/set times outside the observation window are flagged with `*`
- LST range is shown in the header for reference but not used for filtering
- `target_meridian_transit_time` uses `which="nearest"` relative to the Sun's anti-meridian (midpoint of night)

**`custom_objects.json`** — JSON array of objects with `name`, `ra`, and `dec` fields. RA/Dec must be strings parsable by `astropy.coordinates.SkyCoord` (e.g., `"15h08m09.14s"`, `"+39d58m12.9s"`).

Output files are named `Celestial_Objects_Visibility_YYYYMMDD_HHMMSS.txt` and saved in the project root.
