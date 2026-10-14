"""Tests for QC method correctness — validates actual flag values, not just column presence."""
from pathlib import Path

import numpy as np
import pytest

from hfradarpy.radials import Radial

data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()


@pytest.fixture
def codar_radial():
    """Load a CODAR radial for QC testing."""
    radial_file = data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0100.ruv"
    rad = Radial(radial_file, replace_invalid=False)
    rad.initialize_qc()
    return rad


@pytest.fixture
def codar_radial_previous():
    """Load the previous timestep radial for temporal tests."""
    return data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0000.ruv"


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


class TestMaxVelocityLogic:
    """Test that maximum velocity QC correctly classifies velocities."""

    def test_all_pass_with_high_thresholds(self, codar_radial):
        """With very high thresholds, all data should pass."""
        codar_radial.qc_qartod_maximum_velocity(max_speed=99999, high_speed=99998)
        assert (codar_radial.data["Q202"] == 1).all()

    def test_all_fail_with_zero_threshold(self, codar_radial):
        """With zero threshold, all non-zero velocities should fail."""
        codar_radial.qc_qartod_maximum_velocity(max_speed=0, high_speed=0)
        non_zero = codar_radial.data["VELO"].abs() > 0
        if non_zero.any():
            assert (codar_radial.data.loc[non_zero, "Q202"] == 4).all()


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
        # Account for VFLG filtering: num_radials excludes VFLG==128
        if "VFLG" in codar_radial.data:
            count = len(codar_radial.data[codar_radial.data["VFLG"] != 128])
        else:
            count = len(codar_radial.data)
        codar_radial.qc_qartod_radial_count(min_count=count - 1, low_count=count + 1)
        assert (codar_radial.data["Q204"] == 3).all()


class TestAvgBearingLogic:
    def test_pass_with_matching_bearing(self, codar_radial):
        """Reference bearing close to actual mean should pass."""
        actual_mean = codar_radial.data["BEAR"].mean()
        codar_radial.qc_qartod_avg_radial_bearing(
            reference_bearing=actual_mean,
            warning_threshold=15,
            failure_threshold=30
        )
        assert (codar_radial.data["Q207"] == 1).all()

    def test_fail_with_opposite_bearing(self, codar_radial):
        """Reference bearing 180 degrees off should fail."""
        actual_mean = codar_radial.data["BEAR"].mean()
        opposite = (actual_mean + 180) % 360
        codar_radial.qc_qartod_avg_radial_bearing(
            reference_bearing=opposite,
            warning_threshold=15,
            failure_threshold=30
        )
        flag = codar_radial.data["Q207"].iloc[0]
        assert flag in (3, 4)  # should be suspect or fail


class TestPrimaryFlagLogic:
    def test_primary_flag_worst_case(self, codar_radial):
        """Primary flag should be the worst (highest) flag value across all tests."""
        codar_radial.qc_qartod_syntax()
        codar_radial.qc_qartod_maximum_velocity()
        codar_radial.qc_qartod_radial_count()
        codar_radial.qc_qartod_primary_flag()

        qc_cols = ["Q201", "Q202", "Q204"]
        for _, row in codar_radial.data.iterrows():
            worst = max(row[col] for col in qc_cols)
            assert row["PRIM"] == worst

    def test_primary_flag_with_subset(self, codar_radial):
        """Primary flag computed from a subset of tests should only consider those tests."""
        codar_radial.qc_qartod_syntax()
        codar_radial.qc_qartod_maximum_velocity(max_speed=0, high_speed=0)  # force failures
        codar_radial.qc_qartod_radial_count()

        # Only include Q201 (should pass), excluding Q202 (which we forced to fail)
        codar_radial.qc_qartod_primary_flag(include=["Q201", "Q204"])
        # Primary flag should NOT reflect Q202 failures
        assert not (codar_radial.data["PRIM"] == 4).all() or len(codar_radial.data) == 0


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
        radial_file = data_path / "radials" / "ruv" / "SEAB" / "RDLi_SEAB_2019_01_01_0100.ruv"
        rad = Radial(radial_file)
        rad.initialize_qc()
        assert isinstance(rad.metadata["QCTest"], dict)
        assert len(rad.metadata["QCTest"]) == 0
