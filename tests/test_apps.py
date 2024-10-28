import pathlib

import pytest

from cogent3 import load_aligned_seqs
from cogent3.util.deserialise import deserialise_object
from cogent3.util.table import Table
from numpy import allclose

from phylo_limits.apps import (
    PhyloLimitRec,
    check_fit_boundary,
    classify_model_psubs,
    load_param_values,
    load_psubs,
    phylo_limits,
    phylo_limits_tree_to_likelihoodfunction,
)
from phylo_limits.check_boundary import BoundsViolation, ParamRules
from phylo_limits.classify_matrix import ModelMatrixCategories, ModelPsubs


DATADIR = pathlib.Path(__file__).parent / "data"

# set alignment for computing likelihood
_algn = load_aligned_seqs(f"{DATADIR}/piqtree2/four_otu.fasta", moltype="dna")


def test_load_param_values():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    result = load_param_values(model_res)
    assert isinstance(result, ParamRules) == True


def test_load_psubs():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    result = load_psubs(model_res)
    assert isinstance(result, ModelPsubs) == True


def test_generate_record():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )  # two I
    rec_app = phylo_limits()  # default `strict` == F
    record = rec_app(model_res)
    assert isinstance(record, PhyloLimitRec) == True
    assert record.check.strict == False


def test_generate_record_strict_control():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )  # two I
    rec_app = phylo_limits(strict=True)  # set `strict`
    record = rec_app(model_res)
    assert isinstance(record, PhyloLimitRec) == True
    assert record.check.strict == True
    assert record.is_identifiable == False


def test_to_rich_dict_phylolimitrec():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    rec_app = phylo_limits(strict=True)
    record = rec_app(model_res)
    result = record.to_rich_dict()
    assert isinstance(result, dict) == True
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
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    rec_app = phylo_limits()
    record = rec_app(model_res)
    assert record.violation_type == None


def test_has_bv_phylolimitrec():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    rec_app = phylo_limits()
    record = rec_app(model_res)
    assert record.has_BV == True


def test_to_table_phylolimitrec():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    rec_app = phylo_limits()
    record = rec_app(model_res)
    result = record.to_table()
    assert isinstance(result, Table) == True


def test_repr_html_phylolimitrec():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    rec_app = phylo_limits()
    record = rec_app(model_res)
    result = record._repr_html_()
    assert isinstance(result, str) == True


def test_classify_model_psubs():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    check_app = classify_model_psubs()
    res = check_app(model_res)
    assert isinstance(res, ModelMatrixCategories) == True


def test_check_fit_boundary():
    model_res = deserialise_object(
        f"{DATADIR}/eval_identifiability/unid_model_result.json"
    )
    check_app = check_fit_boundary()
    res = check_app(model_res)
    assert isinstance(res, BoundsViolation) == True


@pytest.mark.parametrize("tree_name", ["hky_tree", "gtr_tree"])
def test_convert_piqtree_to_lf(tree_name):
    tree = deserialise_object(f"{DATADIR}/piqtree2/{tree_name}.json")
    converter = phylo_limits_tree_to_likelihoodfunction()
    lf = converter(tree)
    lf.set_alignment(_algn)
    assert allclose(lf.lnL, tree.params["lnL"])
