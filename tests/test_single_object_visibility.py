###
# test_single_object_visibility.py
# John Phillips, john.d.phillips@comcast.net
# Creation date: 2026-10-29
###

# Checks observing-run boundaries, summed intervals, interactive selection,
# and astronomical agreement with v3 using bundled, offline Earth-rotation data.

import contextlib
import io
import json
import tempfile
import unittest
from datetime import date, timedelta
from pathlib import Path
from unittest.mock import patch

import astropy.units as u
import numpy as np
from astroplan import Observer
from astroplan.exceptions import TargetAlwaysUpWarning
from astropy.coordinates import EarthLocation
from astropy.time import Time
from astropy.utils import iers

import Single_Object_Visibility_v1 as visibility


# --- make_night ---
# Creates a synthetic observing night of a specified duration in hours.
def make_night(night_date, hours):
    if hours is None:
        return visibility.NightVisibility(night_date, [])
    start = Time(f'{night_date}T20:00:00')
    return visibility.NightVisibility(night_date, [(start, start + hours * u.hour)])


class VisibilityTests(unittest.TestCase):
    # --- setUpClass ---
    # Uses bundled IERS data so validation does not require network access.
    @classmethod
    def setUpClass(cls):
        cls.iers_context = iers.conf.set_temp('auto_download', False)
        cls.iers_context.__enter__()
        cls.age_context = iers.conf.set_temp('auto_max_age', None)
        cls.age_context.__enter__()
        cls.observer = Observer(
            location=EarthLocation.from_geodetic(-99.3827 * u.deg, 31.5473 * u.deg, 446 * u.m),
            timezone='America/Chicago',
        )
        cls.origin = date(2026, 1, 1)

    # --- tearDownClass ---
    # Restores the caller's IERS settings after offline validation.
    @classmethod
    def tearDownClass(cls):
        cls.age_context.__exit__(None, None, None)
        cls.iers_context.__exit__(None, None, None)

    # --- test_next_run_ends_before_second_season ---
    # Verifies that an instant starts a run and later seasons are excluded.
    def test_next_run_ends_before_second_season(self):
        values = {-1: None, 0: None, 1: 0, 2: 3, 3: 0, 4: None, 5: 9}

        # --- night_at ---
        # Exposes synthetic nights and fails if the search reaches the second run.
        def night_at(day):
            offset = (day - self.origin).days
            self.assertLess(offset, 5)
            return make_night(day, values[offset])

        result = visibility.find_visibility_run(night_at, self.origin)
        self.assertTrue(result.complete)
        self.assertFalse(result.underway)
        self.assertEqual([night.hours for night in result.nights], [0, 3, 0])
        self.assertEqual(result.first.night_date, self.origin + timedelta(days=1))

    # --- test_ongoing_run_and_midpoint_tie ---
    # Retains the current season and selects the middle of its maximum plateau.
    def test_ongoing_run_and_midpoint_tie(self):
        values = {-1: 2, 0: 4, 1: 5, 2: 5, 3: 5, 4: 1, 5: None}
        result = visibility.find_visibility_run(
            lambda day: make_night(day, values[(day - self.origin).days]), self.origin
        )
        self.assertTrue(result.underway)
        self.assertEqual(visibility.maximum_night(result.nights).night_date, self.origin + timedelta(days=2))
        self.assertIn('2026-01-01*', visibility.format_summary(result, self.observer, self.origin)[0])

    # --- test_even_tie_uses_earlier_middle ---
    # Selects the earlier central date rather than the first of four tied nights.
    def test_even_tie_uses_earlier_middle(self):
        nights = [make_night(self.origin + timedelta(days=offset), 4) for offset in range(4)]
        self.assertEqual(visibility.maximum_night(nights).night_date, self.origin + timedelta(days=1))

    # --- test_ongoing_night_clipping_and_full_night_maximum ---
    # Clips first availability while preserving full-night duration for the maximum.
    def test_ongoing_night_clipping_and_full_night_maximum(self):
        cutoff = Time('2026-01-01T22:00:00')
        values = {-1: None, 0: 4, 1: None}
        result = visibility.find_visibility_run(
            lambda day: make_night(day, values[(day - self.origin).days]), self.origin, cutoff
        )
        self.assertTrue(result.underway)
        self.assertEqual(result.first.intervals[0][0], cutoff)
        self.assertAlmostEqual(result.first.hours, 2)
        self.assertAlmostEqual(visibility.maximum_night(result.nights).hours, 4)

    # --- test_exhausted_first_night_moves_to_next_available_window ---
    # Skips elapsed windows without skipping the ongoing run's next night.
    def test_exhausted_first_night_moves_to_next_available_window(self):
        values = {-1: 2, 0: 2, 1: 2, 2: None}
        result = visibility.find_visibility_run(
            lambda day: make_night(day, values[(day - self.origin).days]),
            self.origin, Time('2026-01-01T23:00:00'),
        )
        self.assertEqual(result.first.night_date, self.origin + timedelta(days=1))
        self.assertTrue(result.underway)

    # --- test_cycle_safeguards ---
    # Distinguishes never-found visibility from a run without a found seasonal end.
    def test_cycle_safeguards(self):
        for hours in (None, 8):
            with self.subTest(hours=hours):
                result = visibility.find_visibility_run(
                    lambda day: make_night(day, hours), self.origin, cycle_nights=5
                )
                lines = visibility.format_summary(result, self.observer, self.origin)
                self.assertEqual(result.checked, 5)
                self.assertFalse(result.complete)
                self.assertEqual(len(lines), 3)
                self.assertIn('None found' if hours is None else 'No seasonal disappearance', '\n'.join(lines))

    # --- test_search_limits_are_separate ---
    # Allows a late-starting run its own complete cycle of endpoint searching.
    def test_search_limits_are_separate(self):
        result = visibility.find_visibility_run(
            lambda day: make_night(day, 2 if 4 <= (day - self.origin).days <= 7 else None),
            self.origin, cycle_nights=5,
        )
        self.assertTrue(result.complete)
        self.assertEqual(len(result.nights), 4)

    # --- test_separate_windows_exclude_the_gap ---
    # Checks two dark-window intersections against a curve with known crossings.
    def test_separate_windows_exclude_the_gap(self):
        start = Time('2026-01-01T00:00:00')
        knots = [start + hours * u.hour for hours in (0, 5, 10)]

        # --- altitude_at ---
        # Places crossings at hours three and seven, leaving a four-hour gap.
        def altitude_at(moment):
            return ((moment - start).to_value(u.hour) - 5) ** 2 - 4

        intervals = visibility.intervals_above_limit(knots, [altitude_at(t) for t in knots], altitude_at, 0)
        night = visibility.NightVisibility(self.origin, intervals)
        self.assertEqual(len(intervals), 2)
        self.assertAlmostEqual(night.hours, 6, delta=0.00006)

    # --- test_real_night_has_two_separate_windows ---
    # Confirms both sides of a real lower-transit dip are counted independently.
    def test_real_night_has_two_separate_windows(self):
        target = visibility.make_target('Split northern field', '18h00m00s', '+80d00m00s')
        night = visibility.calculate_night(self.observer, target, 25, date(2026, 1, 15))
        self.assertEqual(len(night.intervals), 2)
        start, end = night.intervals[0][0], night.intervals[-1][1]
        self.assertLess(night.hours, (end - start).to_value(u.hour) - 4)
        for first, last in night.intervals:
            self.assertGreaterEqual(self.observer.altaz(first + (last - first) / 2, target).alt.deg, 25)

    # --- test_isolated_contact_is_not_discarded ---
    # Preserves a zero-duration contact at a transit or twilight endpoint.
    def test_isolated_contact_is_not_discarded(self):
        start = Time('2026-01-01T00:00:00')
        knots = [start + hours * u.hour for hours in (0, 1, 2)]
        for heights in ([-1, 0, -1], [0, -1, -2], [-2, -1, 0]):
            with self.subTest(heights=heights):
                intervals = visibility.intervals_above_limit(knots, heights, lambda t: -1, 0)
                night = visibility.NightVisibility(self.origin, intervals)
                self.assertEqual(len(intervals), 1)
                self.assertEqual(night.hours, 0)
                self.assertIn('(instant)', visibility.format_night('First visibility:', night, self.observer))

    # --- test_v3_qualification_and_independent_duration ---
    # Compares v3 event checks with the new result and independently integrates
    # one-minute altitudes to detect missing windows or incorrectly counted gaps.
    def test_v3_qualification_and_independent_duration(self):
        targets = [
            ('M31', '00h42m44.3s', '+41d16m09s'),
            ('TZ Boo', '15h08m09.14s', '+39d58m12.9s'),
            ('Northern field', '06h00m00s', '+70d00m00s'),
            ('Southern field', '06h00m00s', '-80d00m00s'),
        ]
        for target_data in targets:
            target = visibility.make_target(*target_data)
            for night_date in (date(2026, 1, 15), date(2026, 7, 15)):
                with self.subTest(target=target.name, night=night_date):
                    evening, morning = visibility.night_bounds(self.observer, night_date)
                    midpoint = Time((evening.jd + morning.jd) / 2, format='jd')
                    transit = self.observer.target_meridian_transit_time(
                        midpoint, target, which='nearest', n_grid_points=360
                    )
                    if not evening <= transit <= morning:
                        next_transit = self.observer.target_meridian_transit_time(
                            evening, target, which='next', n_grid_points=360
                        )
                        if next_transit <= morning:
                            transit = next_transit
                    heights = [self.observer.altaz(t, target).alt.deg for t in (evening, morning)]
                    if evening <= transit <= morning:
                        heights.append(self.observer.altaz(transit, target).alt.deg)
                    night = visibility.calculate_night(self.observer, target, 30, night_date)
                    self.assertEqual(bool(night.intervals), max(heights) >= 30)
                    samples = evening + np.linspace(0, (morning - evening).to_value(u.day), 1441) * u.day
                    above = self.observer.altaz(samples, target).alt.deg >= 30
                    integrated_hours = np.trapezoid(above.astype(float), (samples - evening).to_value(u.hour))
                    self.assertAlmostEqual(night.hours, integrated_hours, delta=0.025)

    # --- test_real_twilight_contact ---
    # Sets the limit exactly at evening altitude when the target is descending.
    def test_real_twilight_contact(self):
        target = visibility.make_target('Evening contact', '00h42m44.3s', '+41d16m09s')
        night_date = date(2026, 2, 15)
        evening, _ = visibility.night_bounds(self.observer, night_date)
        limit = self.observer.altaz(evening, target).alt.deg
        night = visibility.calculate_night(self.observer, target, limit, night_date)
        self.assertTrue(night.intervals)
        self.assertEqual(night.intervals[0][0], evening)
        self.assertEqual(night.hours, 0)

    # --- test_default_start_includes_after_midnight_remainder ---
    # Labels an ongoing night by yesterday's evening date while retaining now.
    def test_default_start_includes_after_midnight_remainder(self):
        now = Time('2026-01-16T08:00:00')
        with patch('builtins.input', return_value=''):
            start_date, cutoff = visibility.select_start(self.observer, now)
        self.assertEqual(start_date, date(2026, 1, 15))
        self.assertEqual(cutoff, now)

    # --- test_dst_nights_use_elapsed_time ---
    # Checks physical observing hours rather than differences in wall-clock labels.
    def test_dst_nights_use_elapsed_time(self):
        target = visibility.make_target('Northern field', '06h00m00s', '+85d00m00s')
        for night_date in (date(2026, 3, 7), date(2026, 10, 31)):
            with self.subTest(night=night_date):
                night = visibility.calculate_night(self.observer, target, 20, night_date)
                evening, morning = visibility.night_bounds(self.observer, night_date)
                self.assertAlmostEqual(night.hours, (morning - evening).to_value(u.hour))

    # --- test_missing_twilight_is_explicit ---
    # Does not treat absent high-latitude twilight as a seasonal endpoint.
    def test_missing_twilight_is_explicit(self):
        polar_observer = Observer(location=EarthLocation.from_geodetic(0 * u.deg, 80 * u.deg), timezone='UTC')
        with self.assertWarns(TargetAlwaysUpWarning):
            with self.assertRaisesRegex(ValueError, 'No evening astronomical twilight'):
                visibility.night_bounds(polar_observer, date(2026, 6, 21))

    # --- test_all_object_input_modes ---
    # Exercises real custom/manual parsing and a mocked online resolver.
    def test_all_object_input_modes(self):
        with tempfile.TemporaryDirectory(dir=visibility.BASE_DIR) as directory:
            Path(directory, 'custom_objects.json').write_text(json.dumps([
                {'name': 'TZ Boo', 'ra': '15h08m09.14s', 'dec': '+39d58m12.9s'}
            ]), encoding='utf-8')
            with contextlib.redirect_stdout(io.StringIO()), patch.object(visibility, 'BASE_DIR', Path(directory)):
                with patch('builtins.input', side_effect=['1', 'TZ Boo']):
                    custom = visibility.select_target()
        with contextlib.redirect_stdout(io.StringIO()):
            with patch('builtins.input', side_effect=['3', 'Manual', '15h08m09.14s', '+39d58m12.9s']):
                manual = visibility.select_target()
            with patch('builtins.input', side_effect=['2', 'M31']):
                with patch.object(visibility.FixedTarget, 'from_name', return_value=manual) as resolver:
                    self.assertIs(visibility.select_target(), manual)
                    resolver.assert_called_once_with('M31')
        self.assertLess(custom.coord.separation(manual.coord).arcsec, 0.0001)

    # --- test_planning_allows_old_predictions_and_restores_settings ---
    # Prevents offline seasonal planning from failing on the default IERS age limit.
    def test_planning_allows_old_predictions_and_restores_settings(self):
        with iers.conf.set_temp('auto_max_age', 30):
            with patch.object(visibility, 'run_application', side_effect=lambda: iers.conf.auto_max_age):
                self.assertIsNone(visibility.main())
            self.assertEqual(iers.conf.auto_max_age, 30)

    # --- test_main_saves_exact_three_line_summary ---
    # Runs the real calculation/output path with a small synthetic seasonal run.
    def test_main_saves_exact_three_line_summary(self):
        values = {-1: 1, 0: 2, 1: 3, 2: None}
        target = visibility.make_target('Test target', '12h00m00s', '+20d00m00s')
        with tempfile.TemporaryDirectory(dir=visibility.BASE_DIR) as directory:
            with contextlib.redirect_stdout(io.StringIO()) as output:
                with patch.object(visibility, 'BASE_DIR', Path(directory)), \
                     patch.object(visibility, 'load_configuration', return_value=(self.observer, 30)), \
                     patch.object(visibility, 'select_target', return_value=target), \
                     patch.object(visibility, 'select_start', return_value=(self.origin, None)), \
                     patch.object(visibility, 'calculate_night', side_effect=lambda o, t, l, day: make_night(day, values[(day - self.origin).days])):
                    self.assertEqual(visibility.main(), 0)
            reports = list((Path(directory) / 'output').glob('*.txt'))
            self.assertEqual(len(reports), 1)
            rows = reports[0].read_text(encoding='utf-8').splitlines()
            self.assertEqual(len(rows), 3)
            self.assertIn('\n'.join(rows), output.getvalue())


if __name__ == '__main__':
    unittest.main()
