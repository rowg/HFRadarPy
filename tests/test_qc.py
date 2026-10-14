"""Tests for QC method correctness — validates actual flag values, not just column presence."""
import shutil
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

from hfradarpy.radials import Radial

data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()
seab_path = data_path / "radials" / "ruv" / "SEAB"


def _load(name="RDLi_SEAB_2019_01_01_0100.ruv"):
    rad = Radial(seab_path / name, replace_invalid=False)
    rad.initialize_qc()
    return rad


def _with_rows(rad, **columns):
    """Keep the first n rows of a radial and overwrite the given columns, for hand-built QC inputs."""
    n = len(next(iter(columns.values())))
    rad.data = rad.data.iloc[:n].copy()
    for col, values in columns.items():
        rad.data[col] = np.asarray(values, dtype=float)
    return rad


@pytest.fixture
def codar_radial():
    """Load a CODAR radial for QC testing."""
    return _load()


@pytest.fixture
def codar_radial_previous():
    """Load the previous timestep radial for temporal tests."""
    return seab_path / "RDLi_SEAB_2019_01_01_0000.ruv"


class TestQCFlagValues:
    """Verify that QC flags contain only valid QARTOD values: 1 (pass), 2 (not evaluated), 3 (suspect), 4 (fail)."""

    VALID_FLAGS = {1, 2, 3, 4}

    def test_syntax_flag_values(self, codar_radial):
        codar_radial.qc_qartod_syntax()
        flags = set(codar_radial.data["Q201"].unique())
        assert flags.issubset(self.VALID_FLAGS), f"Invalid flags found: {flags - self.VALID_FLAGS}"

    def test_max_velocity_flag_values(self, codar_radial):
        codar_radial.qc_qartod_maximum_velocity()
        flags = set(codar_radial.data["Q202"].unique())
        assert flags.issubset(self.VALID_FLAGS)

    def test_valid_location_flag_values(self, codar_radial):
        codar_radial.qc_qartod_valid_location()
        flags = set(codar_radial.data["Q203"].unique())
        assert flags.issubset(self.VALID_FLAGS)

    def test_radial_count_flag_values(self, codar_radial):
        codar_radial.qc_qartod_radial_count()
        flags = set(codar_radial.data["Q204"].unique())
        assert flags.issubset(self.VALID_FLAGS)

    def test_spatial_median_flag_values(self, codar_radial):
        codar_radial.qc_qartod_spatial_median()
        flags = set(codar_radial.data["Q205"].unique())
        assert flags.issubset(self.VALID_FLAGS)

    def test_temporal_gradient_flag_values(self, codar_radial, codar_radial_previous):
        codar_radial.qc_qartod_temporal_gradient(codar_radial_previous)
        flags = set(codar_radial.data["Q206"].unique())
        assert flags.issubset(self.VALID_FLAGS)

    def test_avg_bearing_flag_values(self, codar_radial):
        codar_radial.qc_qartod_avg_radial_bearing(reference_bearing=180)
        flags = set(codar_radial.data["Q207"].unique())
        assert flags.issubset(self.VALID_FLAGS)

    def test_primary_flag_values(self, codar_radial):
        codar_radial.qc_qartod_syntax()
        codar_radial.qc_qartod_maximum_velocity()
        codar_radial.qc_qartod_primary_flag()
        flags = set(codar_radial.data["PRIM"].unique())
        assert flags.issubset(self.VALID_FLAGS)


class TestSyntax:
    def test_valid_file_passes(self, codar_radial):
        codar_radial.qc_qartod_syntax()
        assert (codar_radial.data["Q201"] == 1).all()

    def test_timestamp_mismatch_fails(self, codar_radial):
        """The filename says 01:00, so a header timestamp of 02:00 must fail."""
        codar_radial.metadata["TimeStamp"] = "2019 01 01  02 00 00"
        codar_radial.qc_qartod_syntax()
        assert (codar_radial.data["Q201"] == 4).all()

    def test_origin_out_of_range_fails(self, codar_radial):
        codar_radial.metadata["Origin"] = "40.3668167  -273.9735333"
        codar_radial.qc_qartod_syntax()
        assert (codar_radial.data["Q201"] == 4).all()

    def test_filename_without_timestamp_fails(self, codar_radial):
        codar_radial.file_name = "radial_without_timestamp.ruv"
        codar_radial.qc_qartod_syntax()
        assert (codar_radial.data["Q201"] == 4).all()

    @pytest.mark.parametrize("key", ["PatternType", "TimeStamp", "Origin"])
    def test_missing_metadata_fails(self, codar_radial, key):
        del codar_radial.metadata[key]
        codar_radial.qc_qartod_syntax()
        assert (codar_radial.data["Q201"] == 4).all()

    @pytest.mark.parametrize("name", ["RDL_csw_2019_10_24_162300.ruv", "RDL_UMiami_STF_2019_06_01_0000.hfrweralluv1.0"])
    def test_wera_file_passes_without_pattern_type(self, name):
        """PatternType only applies to direction-finding systems, so WERA files don't need it."""
        rad = Radial(data_path / "radials" / "ruv" / "WERA" / name)
        rad.initialize_qc()
        rad.qc_qartod_syntax()
        assert (rad.data["Q201"] == 1).all()


class TestMaxVelocityLogic:
    """Test that maximum velocity QC correctly classifies velocities."""

    def test_all_pass_with_high_thresholds(self, codar_radial):
        """With very high thresholds, all data should pass."""
        codar_radial.qc_qartod_maximum_velocity(max_speed=99999, high_speed=99998)
        assert (codar_radial.data["Q202"] == 1).all()

    def test_all_fail_with_zero_threshold(self, codar_radial):
        """With zero threshold, all non-zero velocities should fail."""
        non_zero = codar_radial.data["VELO"].abs() > 0
        assert non_zero.sum() > 100
        codar_radial.qc_qartod_maximum_velocity(max_speed=0, high_speed=0)
        assert (codar_radial.data.loc[non_zero, "Q202"] == 4).all()

    def test_exact_flags(self, codar_radial):
        """|VELO| <= high passes, high < |VELO| < max is suspect, |VELO| > max fails, in either direction."""
        rad = _with_rows(codar_radial, VELO=[10, -150, 151, -249, 251, -251])
        rad.qc_qartod_maximum_velocity(max_speed=250, high_speed=150)
        assert rad.data["Q202"].tolist() == [1, 1, 3, 3, 4, 4]

    def test_speed_equal_to_max_is_suspect(self, codar_radial):
        """Only speeds beyond max_speed fail; max_speed itself is the top of the suspect band."""
        rad = _with_rows(codar_radial, VELO=[250, -250])
        rad.qc_qartod_maximum_velocity(max_speed=250, high_speed=150)
        assert (rad.data["Q202"] == 3).all()


class TestRadialCountLogic:
    def test_pass_with_low_thresholds(self, codar_radial):
        """With thresholds well below data count, should pass."""
        codar_radial.qc_qartod_radial_count(min_count=1, low_count=2)
        assert (codar_radial.data["Q204"] == 1).all()

    def test_fail_with_high_min_count(self, codar_radial):
        """With min_count above data count, should fail."""
        codar_radial.qc_qartod_radial_count(min_count=999999, low_count=999999)
        assert (codar_radial.data["Q204"] == 4).all()

    def test_warn_between_thresholds(self, codar_radial):
        """With data count between min and low, should warn."""
        count = (codar_radial.data["VFLG"] != 128).sum()
        codar_radial.qc_qartod_radial_count(min_count=count - 1, low_count=count + 1)
        assert (codar_radial.data["Q204"] == 3).all()

    @pytest.mark.parametrize("min_count, low_count, expected", [
        (397, 500, 3),  # count == min_count is suspect, not fail
        (398, 500, 4),
        (100, 397, 3),  # count == low_count is still suspect
        (100, 396, 1),
    ])
    def test_vflg_128_rows_are_not_counted(self, codar_radial, min_count, low_count, expected):
        """The file has 733 rows, 336 with VFLG=128, so the count used is 397."""
        assert len(codar_radial.data) == 733
        assert (codar_radial.data["VFLG"] != 128).sum() == 397
        codar_radial.qc_qartod_radial_count(min_count=min_count, low_count=low_count)
        assert (codar_radial.data["Q204"] == expected).all()


class TestValidLocation:
    def test_vflg_128_fails_everything_else_passes(self, codar_radial):
        codar_radial.qc_qartod_valid_location(use_mask=False)
        npt.assert_array_equal(codar_radial.data["Q203"].to_numpy(),
                               np.where(codar_radial.data["VFLG"] == 128, 4, 1))

    def test_angseg_wrapping_through_north(self, codar_radial):
        """A segment from 350 to 20 degrees covers bearings 350-360 and 0-20, within its range limits only."""
        codar_radial.data["VFLG"] = 0
        seg = dict(range_min=10, range_max=30, bearingTrue_start=350, bearingTrue_stop=20)
        codar_radial.qc_qartod_valid_location(use_mask=False, angseg=[seg])
        bear, rnge = codar_radial.data["BEAR"], codar_radial.data["RNGE"]
        in_segment = rnge.between(10, 30) & ((bear >= 350) | (bear <= 20))
        assert in_segment.sum() > 5
        npt.assert_array_equal(codar_radial.data["Q203"].to_numpy(), np.where(in_segment, 4, 1))

    def test_not_evaluated_without_any_method(self, codar_radial):
        codar_radial.data = codar_radial.data.drop(columns="VFLG")
        codar_radial.qc_qartod_valid_location(use_mask=False)
        assert (codar_radial.data["Q203"] == 2).all()


def spatial_median_reference(data, range_cell_limit=2, angular_limit=10, current_difference=30):
    """Brute-force QARTOD spatial median: compare each radial with the median of all radials (itself included)
    within range_cell_limit range cells and angular_limit degrees of bearing (wrapping through north)."""
    bear, cell, velo = (data[c].to_numpy(dtype=float) for c in ("BEAR", "SPRC", "VELO"))
    dbear = np.abs(bear[:, None] - bear[None, :])
    dbear = np.minimum(dbear, 360 - dbear)
    neighbors = (dbear <= angular_limit) & (np.abs(cell[:, None] - cell[None, :]) <= range_cell_limit)
    median = np.array([np.median(velo[n]) for n in neighbors])
    return np.where(np.abs(velo - median) > current_difference, 4, 1)


class TestSpatialMedian:
    @pytest.mark.parametrize("name", ["RDLi_SEAB_2019_01_01_0000.ruv", "RDLi_SEAB_2019_01_01_0100.ruv"])
    def test_matches_reference(self, name):
        rad = _load(name)
        rad.qc_qartod_spatial_median()
        expected = spatial_median_reference(rad.data)
        assert (expected == 4).sum() > 0
        npt.assert_array_equal(rad.data["Q205"].to_numpy(), expected)

    def test_injected_spike_fails(self, codar_radial):
        baseline = spatial_median_reference(codar_radial.data)
        spike = int(np.flatnonzero(baseline == 1)[len(codar_radial.data) // 2])
        codar_radial.data.loc[spike, "VELO"] += 100
        codar_radial.qc_qartod_spatial_median()
        assert codar_radial.data.loc[spike, "Q205"] == 4
        npt.assert_array_equal(codar_radial.data["Q205"].to_numpy(), spatial_median_reference(codar_radial.data))

    def test_uniform_field_passes(self, codar_radial):
        codar_radial.data["VELO"] = 12.5
        codar_radial.qc_qartod_spatial_median()
        assert (codar_radial.data["Q205"] == 1).all()


class TestTemporalGradient:
    def _differences(self, rad, previous_file):
        previous = Radial(previous_file)
        merged = rad.data[["LOND", "LATD", "VELO"]].merge(previous.data[["LOND", "LATD", "VELO"]],
                                                          on=["LOND", "LATD"], how="left", suffixes=("", "_prev"))
        return (merged["VELO"] - merged["VELO_prev"]).abs()

    def test_matches_reference(self, codar_radial, codar_radial_previous):
        """Flags follow the documented thresholds: missing previous value 2, >= fail 4, >= warn 3, otherwise 1."""
        diff = self._differences(codar_radial, codar_radial_previous)
        codar_radial.qc_qartod_temporal_gradient(codar_radial_previous, gradient_temp_fail=54, gradient_temp_warn=36)
        expected = np.select([diff.isna(), diff >= 54, diff >= 36], [2, 4, 3], 1)
        assert set(expected) == {1, 2, 3, 4}  # the data exercises every branch
        npt.assert_array_equal(codar_radial.data["Q206"].to_numpy(), expected)

    def test_missing_previous_file_not_evaluated(self, codar_radial, tmp_path):
        codar_radial.qc_qartod_temporal_gradient(tmp_path / "RDLi_SEAB_2019_01_01_0000.ruv")
        assert (codar_radial.data["Q206"] == 2).all()

    def test_difference_equal_to_fail_threshold_fails(self, codar_radial, codar_radial_previous):
        diff = self._differences(codar_radial, codar_radial_previous)
        k = diff.idxmax()
        codar_radial.qc_qartod_temporal_gradient(codar_radial_previous, gradient_temp_fail=diff[k],
                                                 gradient_temp_warn=diff[k] - 1)
        assert codar_radial.data.loc[k, "Q206"] == 4

    def test_difference_equal_to_warn_threshold_is_suspect(self, codar_radial, codar_radial_previous):
        diff = self._differences(codar_radial, codar_radial_previous)
        k = diff.idxmax()
        codar_radial.qc_qartod_temporal_gradient(codar_radial_previous, gradient_temp_fail=diff[k] + 1,
                                                 gradient_temp_warn=diff[k])
        assert codar_radial.data.loc[k, "Q206"] == 3


class TestStuckValue:
    def test_repeated_file_is_stuck(self, tmp_path):
        """The same radial three hours in a row is stuck everywhere."""
        source = seab_path / "RDLi_SEAB_2019_01_01_0100.ruv"
        for name in ["RDLi_SEAB_2019_01_01_0100.ruv", "RDLi_SEAB_2019_01_01_0000.ruv",
                     "RDLi_SEAB_2018_12_31_2300.ruv"]:
            shutil.copy(source, tmp_path / name)
        rad = Radial(tmp_path / "RDLi_SEAB_2019_01_01_0100.ruv", replace_invalid=False)
        rad.initialize_qc()
        rad.qc_qartod_stuck_value(resolution=0.01, N=3)
        assert (rad.data["Q209"] == 4).all()

    def test_not_evaluated_without_previous_files(self, codar_radial):
        """There is no 2018-12-31 23:00 file, so N=3 cannot be evaluated."""
        codar_radial.qc_qartod_stuck_value(N=3)
        assert (codar_radial.data["Q209"] == 2).all()

    def test_matches_reference(self):
        rad = _load("RDLi_SEAB_2019_01_01_0200.ruv")
        ts = rad.data[["LOND", "LATD", "VELO"]]
        for i, hour in enumerate(["0100", "0000"], start=1):
            previous = Radial(seab_path / f"RDLi_SEAB_2019_01_01_{hour}.ruv").data[["LOND", "LATD", "VELO"]]
            ts = ts.merge(previous, on=["LOND", "LATD"], how="left", suffixes=(None, f"_{i}"))
        velo = ts[["VELO", "VELO_1", "VELO_2"]].to_numpy()
        stuck = np.abs(np.diff(velo, axis=1)).max(axis=1) < 0.01
        expected = np.select([np.isnan(velo).any(axis=1), stuck], [2, 4], 1)

        rad.qc_qartod_stuck_value(resolution=0.01, N=3)
        assert set(expected) >= {1, 2}
        npt.assert_array_equal(rad.data["Q209"].to_numpy(), expected)


class TestAvgBearingLogic:
    BEARINGS = [100, 110, 130, 180]  # mean 130

    def test_pass_with_matching_bearing(self, codar_radial):
        """Reference bearing equal to the mean should pass."""
        rad = _with_rows(codar_radial, BEAR=self.BEARINGS)
        rad.qc_qartod_avg_radial_bearing(reference_bearing=130, warning_threshold=15, failure_threshold=30)
        assert (rad.data["Q207"] == 1).all()

    @pytest.mark.parametrize("offset, expected", [(14.9, 1), (15, 3), (29.9, 3), (30, 4), (-30, 4), (180, 4)])
    def test_exact_flags(self, codar_radial, offset, expected):
        rad = _with_rows(codar_radial, BEAR=self.BEARINGS)
        rad.qc_qartod_avg_radial_bearing(reference_bearing=130 + offset, warning_threshold=15, failure_threshold=30)
        assert (rad.data["Q207"] == expected).all()

    @pytest.mark.parametrize("rotation", [0, 90, 200, 260, 300])
    def test_mean_is_rotation_invariant(self, codar_radial, rotation):
        """Rotating every bearing (including through north) rotates the mean, so the flag is unchanged."""
        rad = _with_rows(codar_radial, BEAR=np.mod(np.array(self.BEARINGS) + rotation, 360))
        rad.qc_qartod_avg_radial_bearing(reference_bearing=(130 + rotation + 20) % 360,
                                         warning_threshold=15, failure_threshold=30)
        assert (rad.data["Q207"] == 3).all()

    @pytest.mark.parametrize("reference, expected", [(0, 1), (360, 1), (20, 3), (340, 3), (180, 4)])
    def test_coverage_through_north(self, codar_radial, reference, expected):
        """Bearings 350-10 average to ~2.5, not ~180, and the reference difference wraps too."""
        rad = _with_rows(codar_radial, BEAR=[350, 355, 5, 10])
        rad.qc_qartod_avg_radial_bearing(reference_bearing=reference, warning_threshold=15, failure_threshold=30)
        assert (rad.data["Q207"] == expected).all()


def primary_flag_reference(data, columns):
    """Worst flag across the given columns, with "not evaluated" (2) ignored."""
    flags = data[columns]
    return np.select([flags.eq(4).any(axis=1), flags.eq(3).any(axis=1)], [4, 3], 1)


class TestPrimaryFlagLogic:
    @pytest.fixture
    def mixed_flags(self, codar_radial, codar_radial_previous):
        """A radial whose QC columns contain every flag value: 1-4."""
        codar_radial.qc_qartod_syntax()                                          # all 1
        codar_radial.qc_qartod_maximum_velocity(max_speed=30, high_speed=20)     # 1, 3 and 4
        codar_radial.qc_qartod_temporal_gradient(codar_radial_previous)          # includes 2
        codar_radial.qc_qartod_radial_count(min_count=100, low_count=500)        # all 3
        return codar_radial

    def test_primary_flag_worst_case(self, mixed_flags):
        """Primary flag should be the worst flag across all tests, ignoring "not evaluated"."""
        columns = ["Q201", "Q202", "Q206", "Q204"]
        assert set(np.unique(mixed_flags.data[columns])) == {1, 2, 3, 4}
        mixed_flags.qc_qartod_primary_flag()
        npt.assert_array_equal(mixed_flags.data["PRIM"].to_numpy(), primary_flag_reference(mixed_flags.data, columns))

    def test_primary_flag_with_subset(self, mixed_flags):
        """`include` takes QC method names; excluded tests must not affect the primary flag."""
        mixed_flags.qc_qartod_primary_flag(include=["qc_qartod_syntax", "qc_qartod_temporal_gradient"])
        expected = primary_flag_reference(mixed_flags.data, ["Q201", "Q206"])
        assert not (expected == 4).all()  # Q202 failures were excluded
        npt.assert_array_equal(mixed_flags.data["PRIM"].to_numpy(), expected)


class TestQCMetadata:
    def test_qc_metadata_stored(self, codar_radial):
        """QC test metadata should be recorded."""
        codar_radial.qc_qartod_syntax()
        codar_radial.qc_qartod_maximum_velocity()
        assert "Q201" in codar_radial.metadata["QCTest"]
        assert "Q202" in codar_radial.metadata["QCTest"]
        assert "qc_qartod_syntax" in codar_radial.metadata["QCTest"]["Q201"]
        assert "qc_qartod_maximum_velocity" in codar_radial.metadata["QCTest"]["Q202"]

    def test_initialize_qc_creates_empty_dict(self):
        radial_file = seab_path / "RDLi_SEAB_2019_01_01_0100.ruv"
        rad = Radial(radial_file)
        rad.initialize_qc()
        assert isinstance(rad.metadata["QCTest"], dict)
        assert len(rad.metadata["QCTest"]) == 0
