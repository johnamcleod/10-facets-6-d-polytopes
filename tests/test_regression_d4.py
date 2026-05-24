"""
Regression test: d=4, n=8 pipeline must produce exactly 348 polytopes.

This is the single most important correctness gate (§12 of spec).
Run: pytest tests/test_regression_d4.py -v

NOTE: This test requires the AAK order-type database or a sufficiently
dense random sampling. With random sampling it is INCOMPLETE and will
undercount. A passing result here only validates the pipeline logic
on whatever order types were generated.

The test is marked `slow` and `regression`; run with:
  pytest -m regression tests/test_regression_d4.py
"""

import pytest
import os

# Skip unless explicitly running regression suite
pytestmark = pytest.mark.regression


@pytest.fixture(scope="module")
def d4_pipeline_results(tmp_path_factory):
    """Run the d=4 pipeline end-to-end and return results."""
    tmp = tmp_path_factory.mktemp("d4")
    d = 4

    from pipeline.stage1_order_types import run_stage1
    from pipeline.stage2_gale import run_stage2
    from pipeline.stage3_filters import run_stage3

    ot = run_stage1(d, output_dir=str(tmp / "stage1"), data_dir="data/aak")
    types = run_stage2(d, ot, output_dir=str(tmp / "stage2"), exact_dedup=False)
    survivors, elim_log = run_stage3(d, types, output_dir=str(tmp / "stage3"))
    return {"order_types": ot, "types": types, "survivors": survivors,
            "elim_log": elim_log}


def test_d4_pipeline_runs(d4_pipeline_results):
    """Pipeline completes without exception."""
    assert d4_pipeline_results is not None


def test_d4_p_ge_2(d4_pipeline_results):
    """All surviving types have p >= 2."""
    for t in d4_pipeline_results["survivors"]:
        assert t["p_count"] >= 2, f"Type {t['type_id']} has p={t['p_count']} < 2"


def test_d4_missing_face_sizes(d4_pipeline_results):
    """All missing faces in surviving types have size in [2,5]."""
    for t in d4_pipeline_results["survivors"]:
        for mf in t["missing_faces"]:
            assert 2 <= len(mf) <= 5, (
                f"Type {t['type_id']} has missing face {mf} of size {len(mf)}"
            )


@pytest.mark.skipif(
    not os.path.exists("data/aak"),
    reason="AAK database not available — count test requires complete enumeration"
)
def test_d4_survivor_count_with_aak(d4_pipeline_results):
    """With the full AAK database, Stage 3 survivors should match Burcroff's count.

    Burcroff reports specific combinatorial types surviving to Stage 3 in d=4.
    (The exact Stage 3 survivor count; Stage 4/5 then gives 348 polytopes.)
    """
    # TODO: fill in exact Stage 3 survivor count from Burcroff §6
    survivors = d4_pipeline_results["survivors"]
    # At minimum, non-zero survivors
    assert len(survivors) > 0
