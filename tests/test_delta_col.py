import numpy
import pytest
from cogent3.util.dict_array import DictArray

from phylim.classify_matrix import ModelPsubs
from phylim.delta_col import calc_delta_col


@pytest.mark.parametrize(
    "matrix,expected_delta_col",
    [
        (
            numpy.array(
                [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 0.9, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.1, 0.0, 1.0],
                ]
            ),
            0.8,
        ),
        (
            numpy.array(
                [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 0.1, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.9, 0.0, 1.0],
                ]
            ),
            -0.8,
        ),
        (
            numpy.array(
                [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 0.5, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.5, 0.0, 1.0],
                ]
            ),
            0.0,
        ),
    ],
)
def test_calc_delta_col_hand_constructed(matrix, expected_delta_col):
    dict_array = DictArray(matrix)
    psubs = ModelPsubs(psubs={("test_edge",): dict_array}, source="test")
    result = calc_delta_col(psubs)
    assert numpy.isclose(result[("test_edge",)], expected_delta_col)
