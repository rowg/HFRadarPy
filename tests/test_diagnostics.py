"""Tests for hfradarpy.diagnostics module."""
import os
import shutil
from glob import glob
from pathlib import Path

import pandas as pd
import pytest

from hfradarpy.diagnostics import Diags, concat

data_path = (Path(__file__).parent.with_name("examples") / "data").resolve()


class TestDiagsLoading:
    def test_load_hdt_file(self):
        hdt_file = data_path / "diagnostics" / "SEAB" / "STAT_SEAB_2018_01_01.hdt"
        d = Diags(hdt_file)
        assert not d.data.empty
        assert d.file_name == "STAT_SEAB_2018_01_01.hdt"

    def test_load_rdt_file(self):
        rdt_file = data_path / "diagnostics" / "SEAB" / "STAT_SEAB_2018_01_01.rdt"
        d = Diags(rdt_file)
        assert not d.data.empty
        assert d.file_name == "STAT_SEAB_2018_01_01.rdt"

    def test_load_xdt_file(self):
        xdt_file = data_path / "diagnostics" / "SEAB" / "STAT_SEAB_2018_01_01.xdt"
        d = Diags(xdt_file)
        assert not d.data.empty
        assert d.file_name == "STAT_SEAB_2018_01_01.xdt"

    def test_repr(self):
        hdt_file = data_path / "diagnostics" / "SEAB" / "STAT_SEAB_2018_01_01.hdt"
        d = Diags(hdt_file)
        assert "Diagnostic" in repr(d)
        assert "STAT_SEAB_2018_01_01.hdt" in repr(d)

    def test_data_has_time_index(self):
        hdt_file = data_path / "diagnostics" / "SEAB" / "STAT_SEAB_2018_01_01.hdt"
        d = Diags(hdt_file)
        assert d.data.index.name == "time"

    def test_metadata_populated(self):
        hdt_file = data_path / "diagnostics" / "SEAB" / "STAT_SEAB_2018_01_01.hdt"
        d = Diags(hdt_file)
        assert "FileType" in d.metadata
        assert "TimeStamp" in d.metadata

    def test_is_valid(self):
        hdt_file = data_path / "diagnostics" / "SEAB" / "STAT_SEAB_2018_01_01.hdt"
        d = Diags(hdt_file)
        assert d.is_valid()


class TestDiagsConcat:
    def test_concat_hdt_files(self):
        files = sorted(glob(str(data_path / "diagnostics" / "SEAB" / "*.hdt")))
        assert len(files) == 5
        result = concat(files)
        assert isinstance(result, pd.DataFrame)
        assert not result.empty
        # Concatenated result should have more rows than a single file
        single = Diags(files[0])
        assert len(result) > len(single.data)

    def test_concat_rdt_files(self):
        files = sorted(glob(str(data_path / "diagnostics" / "SEAB" / "*.rdt")))
        result = concat(files)
        assert isinstance(result, pd.DataFrame)
        assert not result.empty


class TestDiagsExport:
    @pytest.mark.xfail(reason="Known bug: _TableHeader column count mismatch in to_file()")
    def test_to_file(self, tmp_path):
        hdt_file = data_path / "diagnostics" / "SEAB" / "STAT_SEAB_2018_01_01.hdt"
        d = Diags(hdt_file)
        out_file = tmp_path / "STAT_SEAB_2018_01_01_test.hdt"
        d.to_file(out_file)
        assert os.path.isfile(out_file)

    @pytest.mark.xfail(reason="Known bug: _TableHeader column count mismatch in to_file()")
    def test_to_file_no_overwrite(self, tmp_path):
        """Exporting to the same filename should create a .mod version."""
        # Work on a copy so a failing export never leaves files in examples/data
        hdt_file = tmp_path / "STAT_SEAB_2018_01_01.hdt"
        shutil.copy(data_path / "diagnostics" / "SEAB" / hdt_file.name, hdt_file)
        d = Diags(hdt_file)
        d.to_file(hdt_file, overwrite=False)
        mod_file = hdt_file.with_suffix(".hdt.mod")
        assert os.path.isfile(mod_file)
