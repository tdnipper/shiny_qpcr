import numpy as np
import pandas as pd
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from shared import (
    filter_targets,
    group_ct_data,
    welch_ttest,
    one_sample_ttest,
    correct_pvalues,
    p_to_star,
    export_to_csv,
    export_to_excel,
    _split_sample_name,
    _transform_quantstudio_df,
    _is_quantstudio_export,
    import_file,
)


@pytest.fixture
def sample_df():
    return pd.DataFrame({
        "Group": ["Grp1", "Grp1", "Grp1", "Grp1"],
        "Condition": ["A", "A", "B", "B"],
        "Target Name": ["Gene1", "Gene1", "Gene1", "Gene1"],
        "CT": [20.0, 21.0, 22.0, 23.0],
    })


class TestFilterTargets:
    def test_exact_match(self, sample_df):
        result = filter_targets(sample_df, "Gene1")
        assert len(result) == 4

    def test_no_match(self, sample_df):
        result = filter_targets(sample_df, "GeneX")
        assert len(result) == 0


class TestGroupCtData:
    def test_groups_correctly(self, sample_df):
        result = group_ct_data(sample_df)
        assert "Gene1_Grp1_A" in result
        assert "Gene1_Grp1_B" in result
        assert len(result["Gene1_Grp1_A"]) == 2
        assert len(result["Gene1_Grp1_B"]) == 2


class TestWelchTtest:
    def test_identical_groups(self):
        g1 = np.array([20.0, 20.0, 20.0])
        g2 = np.array([20.0, 20.0, 20.0])
        t, p = welch_ttest(g1, g2)
        assert np.isnan(t)  # 0/0 with no variance

    def test_different_groups(self):
        g1 = np.array([10.0, 11.0, 12.0])
        g2 = np.array([20.0, 21.0, 22.0])
        t, p = welch_ttest(g1, g2)
        assert p < 0.05

    def test_single_replicate_returns_nan(self):
        t, p = welch_ttest(np.array([10.0]), np.array([20.0]))
        assert np.isnan(p)


class TestOneSampleTtest:
    def test_mean_zero(self):
        vals = np.array([0.0, 0.0, 0.0])
        t, p = one_sample_ttest(vals, popmean=0)
        assert np.isnan(t)  # 0/0

    def test_mean_nonzero(self):
        vals = np.array([5.0, 5.1, 4.9])
        t, p = one_sample_ttest(vals, popmean=0)
        assert p < 0.05


class TestCorrectPvalues:
    def test_none_correction(self):
        p = np.array([0.01, 0.04, 0.1])
        result = correct_pvalues(p, method="none")
        np.testing.assert_array_equal(result, p)

    def test_bonferroni(self):
        p = np.array([0.01, 0.04, 0.1])
        result = correct_pvalues(p, method="bonferroni")
        np.testing.assert_allclose(result, [0.03, 0.12, 0.3])

    def test_bonferroni_caps_at_one(self):
        p = np.array([0.5, 0.8])
        result = correct_pvalues(p, method="bonferroni")
        assert all(r <= 1.0 for r in result)

    def test_bh_correction(self):
        p = np.array([0.01, 0.04, 0.1])
        result = correct_pvalues(p, method="bh")
        # BH should be less aggressive than Bonferroni
        bonf = correct_pvalues(p, method="bonferroni")
        assert all(result <= bonf + 1e-10)


class TestPToStar:
    def test_significant(self):
        assert p_to_star(0.0001) == "***"
        assert p_to_star(0.005) == "**"
        assert p_to_star(0.03) == "*"
        assert p_to_star(0.1) == "ns"

    def test_nan(self):
        assert p_to_star(np.nan) == ""


class TestExport:
    def test_csv_export(self, sample_df):
        result = export_to_csv(sample_df)
        assert isinstance(result, bytes)
        assert b"Target Name" in result

    def test_excel_export(self, sample_df):
        result = export_to_excel(sample_df)
        assert isinstance(result, bytes)
        assert len(result) > 0


class TestSplitSampleName:
    def test_simple_split(self):
        group, condition = _split_sample_name("WT_treated")
        assert group == "WT"
        assert condition == "treated"

    def test_underscore_in_group(self):
        group, condition = _split_sample_name("anti_Flag_enrich")
        assert group == "anti_Flag"
        assert condition == "enrich"

    def test_no_underscore_raises(self):
        with pytest.raises(ValueError, match="no underscore"):
            _split_sample_name("nosplit")

    def test_single_underscore(self):
        group, condition = _split_sample_name("A_B")
        assert group == "A"
        assert condition == "B"


class TestTransformQuantstudioDf:
    @pytest.fixture
    def raw_df(self):
        return pd.DataFrame({
            "Sample Name": ["WT_ctrl", "WT_ctrl", "KO_treated"],
            "Target Name": ["Gene1", "Gene2", "Gene1"],
            "CT": [20.0, 21.0, 22.0],
            "Task": ["Unknown", "Unknown", "Unknown"],
        })

    def test_columns_present(self, raw_df):
        result = _transform_quantstudio_df(raw_df)
        assert set(["Group", "Condition", "Target Name", "CT", "Biological Replicate"]).issubset(result.columns)

    def test_sample_name_removed(self, raw_df):
        result = _transform_quantstudio_df(raw_df)
        assert "Sample Name" not in result.columns

    def test_group_condition_values(self, raw_df):
        result = _transform_quantstudio_df(raw_df)
        assert list(result["Group"]) == ["WT", "WT", "KO"]
        assert list(result["Condition"]) == ["ctrl", "ctrl", "treated"]

    def test_task_renamed_to_biological_replicate(self, raw_df):
        result = _transform_quantstudio_df(raw_df)
        assert "Biological Replicate" in result.columns
        assert "Task" not in result.columns

    def test_input_not_mutated(self, raw_df):
        original_cols = list(raw_df.columns)
        _transform_quantstudio_df(raw_df)
        assert list(raw_df.columns) == original_cols


class TestIsQuantstudioExport:
    def test_csv_returns_false(self, tmp_path):
        f = tmp_path / "data.csv"
        f.write_text("a,b\n1,2\n")
        assert _is_quantstudio_export(str(f)) is False

    def test_xlsx_without_results_sheet(self, tmp_path):
        import openpyxl
        wb = openpyxl.Workbook()
        wb.active.title = "Sheet1"
        path = str(tmp_path / "no_results.xlsx")
        wb.save(path)
        assert _is_quantstudio_export(path) is False

    def test_xlsx_with_results_sheet(self, tmp_path):
        import openpyxl
        wb = openpyxl.Workbook()
        wb.active.title = "Results"
        path = str(tmp_path / "with_results.xlsx")
        wb.save(path)
        assert _is_quantstudio_export(path) is True


def _build_quantstudio_xlsx(tmp_path, rows):
    """Write a minimal QuantStudio xlsx: 46 blank header rows then data rows."""
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Results"
    # 46 blank header rows
    for _ in range(46):
        ws.append([])
    # Header row (row 47 = skiprows=46 reads from here)
    ws.append(["Sample Name", "Target Name", "CT", "Task"])
    for row in rows:
        ws.append(row)
    path = str(tmp_path / "instrument.xlsx")
    wb.save(path)
    return path


class TestImportFileQuantstudioFormat:
    def test_expected_columns(self, tmp_path):
        path = _build_quantstudio_xlsx(tmp_path, [
            ["WT_ctrl", "Gene1", 20.0, "Unknown"],
        ])
        df = import_file(path)
        assert set(["Group", "Condition", "Target Name", "CT", "Biological Replicate"]).issubset(df.columns)

    def test_group_condition_parsed(self, tmp_path):
        path = _build_quantstudio_xlsx(tmp_path, [
            ["anti_Flag_enrich", "Gene1", 20.0, "Unknown"],
        ])
        df = import_file(path)
        assert df["Group"].iloc[0] == "anti_Flag"
        assert df["Condition"].iloc[0] == "enrich"

    def test_undetermined_row_dropped(self, tmp_path):
        path = _build_quantstudio_xlsx(tmp_path, [
            ["WT_ctrl", "Gene1", 20.0, "Unknown"],
            ["WT_ctrl", "Gene1", "Undetermined", "Unknown"],
        ])
        df = import_file(path)
        assert len(df) == 1

    def test_biological_replicate_populated(self, tmp_path):
        path = _build_quantstudio_xlsx(tmp_path, [
            ["WT_ctrl", "Gene1", 20.0, "Unknown"],
        ])
        df = import_file(path)
        assert df["Biological Replicate"].iloc[0] == "Unknown"
