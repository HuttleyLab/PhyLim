import numpy

from phylim.classify_matrix import ModelPsubs


def min_diff_from_diag(
    m: numpy.ndarray, diag_indices: tuple, off_diag_indices: numpy.ndarray
) -> numpy.ndarray:
    """compute difference for each column between diagonal and
    the largest off-diagonal element

    Args:
        m (numpy.ndarray): a matrix
        diag_indices (tuple): diagonal indices
        off_diag_indices (numpy.ndarray): off-diagonal indices in boolean
    """
    return m[diag_indices] - m.max(axis=0, where=off_diag_indices, initial=m.min())


def min_col_diff(
    m: numpy.ndarray,
    diag_indices: tuple,  # = DIAG,
    off_diag_indices: numpy.ndarray,  # = OFFDIAG,
) -> float:
    return min_diff_from_diag(m, diag_indices, off_diag_indices).min()


def calc_delta_col(psubs: ModelPsubs) -> dict:
    """calculate delta_col for given psubs"""
    psubs_dict = {}
    for key, value in psubs.items():
        p = value.to_array()
        shape = p.shape
        mask = numpy.ones(shape, bool)
        mask[numpy.diag_indices(shape[0])] = False
        diag_indices = numpy.diag_indices(shape[0])
        offdiag_indices = mask
        psubs_dict[key] = min_col_diff(p, diag_indices, offdiag_indices)
    return psubs_dict
