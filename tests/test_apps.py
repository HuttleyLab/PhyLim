import pathlib
import sys

import pytest
from cogent3 import get_app, load_aligned_seqs
from cogent3.app.composable import NotCompleted
from cogent3.app.result import model_result
from cogent3.core.table import Table
from cogent3.util.deserialise import deserialise_object
from numpy import allclose

from phylim.apps import (
    PhyloLimitRec,
    _get_lf,
    check_fit_boundary,
    classify_model_psubs,
    load_param_values,
    load_psubs,
    phylim,
    phylim_filter,
    phylim_to_model_result,
)
from phylim.check_boundary import BoundsViolation, ParamRules
from phylim.classify_matrix import ModelMatrixCategories, ModelPsubs

DATADIR = pathlib.Path(__file__).parent / "data"

# set alignment for computing likelihood
_algn = load_aligned_seqs(f"{DATADIR}/piqtree/four_otu.fasta", moltype="dna")


_model_res = deserialise_object(
    f"{DATADIR}/eval_identifiability/unid_model_result.json"
)

_model_res_split = deserialise_object(
    f"{DATADIR}/split_codon_model_result/split_codon_model_result.json"
)


def test_get_lf():
    assert _get_lf(_model_res) is _model_res.lf


def test_lf_get_lf():
    assert _get_lf(_model_res.lf) is _model_res.lf


def test_get_lf_multiple():
    with pytest.raises(
        ValueError, match="Model result must contain exactly one likelihood function."
    ):
        _get_lf(_model_res_split)


def test_load_param_values():
    result = load_param_values(_model_res.lf)
    assert isinstance(result, ParamRules)


def test_load_psubs():
    result = load_psubs(_model_res.lf)
    assert isinstance(result, ModelPsubs)


def test_generate_record():
    # two I
    rec_app = phylim()  # default `strict` == F
    record = rec_app(_model_res)
    assert isinstance(record, PhyloLimitRec)
    assert not record.check.strict


def test_generate_record_strict_control():
    # two I
    rec_app = phylim(strict=True)  # set `strict`
    record = rec_app(_model_res)
    assert isinstance(record, PhyloLimitRec)
    assert record.check.strict
    assert not record.is_identifiable


def test_to_rich_dict_phylolimitrec():
    rec_app = phylim(strict=True)
    record = rec_app(_model_res)
    result = record.to_rich_dict()
    assert isinstance(result, dict)
    assert all(
        k in result
        for k in [
            "model_name",
            "boundary_values",
            "nondlc_and_identity",
            "source",
            "names",
            "violation_type",
            "strict",
            "version",
        ]
    )


def test_violation_type_phylolimitrec():
    rec_app = phylim()
    record = rec_app(_model_res)
    assert record.violation_type is None


def test_has_bv_phylolimitrec():
    rec_app = phylim()
    record = rec_app(_model_res)
    assert record.has_BV


def test_to_table_phylolimitrec():
    rec_app = phylim()
    record = rec_app(_model_res)
    result = record.to_table()
    assert isinstance(result, Table)


def test_repr_html_phylolimitrec():
    rec_app = phylim()
    record = rec_app(_model_res)
    result = record._repr_html_()
    assert isinstance(result, str)


def test_classify_model_psubs():
    check_app = classify_model_psubs()
    res = check_app(_model_res)
    assert isinstance(res, ModelMatrixCategories)


def test_check_fit_boundary():
    check_app = check_fit_boundary()
    res = check_app(_model_res)
    assert isinstance(res, BoundsViolation)


@pytest.mark.xfail(reason="need new version of piqtree")
@pytest.mark.parametrize("model_name", ["HKY", "GTR", "UNREST"])
def test_convert_piq_build_treeto_model_result(model_name):
    fit = get_app("piq_fit_tree", tree=_algn.quick_tree(), model=model_name)
    tree = fit(_algn)
    converter = phylim_to_model_result()
    res = converter(tree)
    res.lf.set_alignment(_algn)
    assert allclose(res.lf.lnL, tree.params["lnL"])


@pytest.mark.xfail(reason="need new version of piqtree")
@pytest.mark.skipif(
    sys.platform.startswith("win"), reason="Test not supported on Windows"
)
def test_piqtree_app():
    phylo = get_app("piq_build_tree", "GTR")
    tree = phylo(_algn)
    lf_from = get_app("phylim_to_model_result")
    result = lf_from(tree)
    checker = get_app("phylim")
    checked = checker(result)
    assert isinstance(checked, PhyloLimitRec)


def test_phylim_filter_app_pass():
    filter_app = phylim_filter(strict=True)
    result1 = filter_app(_model_res)
    assert isinstance(result1, NotCompleted)


def test_phylim_filter_app_fail():
    filter_app = phylim_filter(strict=False)
    result2 = filter_app(_model_res)
    assert isinstance(result2, model_result)


def test_delta_col():
    rec_app = phylim()
    record = rec_app(_model_res)
    delta_col_dict = record.delta_col
    assert isinstance(delta_col_dict, dict)
    assert all(isinstance(v, float) for v in delta_col_dict.values())


def test_calc_delta_col_tree_coverage():
    tree = _model_res.lf.tree
    expected_edge_count = len(list(tree.preorder(include_self=False)))
    rec_app = phylim()
    record = rec_app(_model_res)
    delta_col_dict = record.delta_col

    assert len(delta_col_dict) == expected_edge_count
