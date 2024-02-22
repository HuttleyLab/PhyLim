#TODO: this file will be filled with functions for projecting the fit: for given likelihood function, report the ens where DLC is broken.
#       (it should project only matrix who is DLC)
import numpy
import scipy

from cogent3 import maths
from cogent3.app.typing import SerialisableType
from cogent3.util.dict_array import DictArray

from phylo_limits.delta_col import (
    EstimateMinColDelta,
    get_stat_pi_via_eigen,
    invalid_dlc,
)
from phylo_limits.p_classifier import check_chainsaw


expm = scipy.linalg.expm



def num_saw_srch(t, q, p0, tau, model_name):
    start = float(tau) * t
    end = 1000  # change the search upper as 1000, step as 0.1
    step = 0.1

    interval_l = 0.0
    intervals = []
    for i in numpy.arange(start, end, step):
        p = scipy.linalg.expm(q.to_array() * i)
        pi = get_stat_pi_via_eigen(p) if model_name in ["ssGN", "GN"] else p0
        if check_chainsaw(p_matrix=p, p_limit=numpy.array([pi, pi, pi, pi])):
            if interval_l == 0.0:
                interval_l = i
            if i + step >= end:
                interval = [
                    maths.matrix_exponential_integration.expected_number_subs(
                        p0=p0, Q=q, t=interval_l
                    ),
                    maths.matrix_exponential_integration.expected_number_subs(
                        p0=p0, Q=q, t=end
                    ),
                ]

                intervals.append(interval)
        elif interval_l != 0.0:
            interval = [
                maths.matrix_exponential_integration.expected_number_subs(
                    p0=p0, Q=q, t=interval_l
                ),
                maths.matrix_exponential_integration.expected_number_subs(
                    p0=p0, Q=q, t=i - step
                ),
            ]
            interval_l = 0.0
            intervals.append(interval)
    return intervals




def get_node_parent(tree: ..., node: str) -> str:
    """given a node, seek its parent in tree"""
    return tree.get_node_matching_name(node).parent.name




def get_sub_num(q:DictArray, p0, model_name, edge_name, qc, lengths) -> SerialisableType:
    """predict the dlc violating point before the limit for a rate matrix run on the edge"""
    # load each Q, no calibration
    estimate = EstimateMinColDelta(q)
    estimate.optimise()

    tau = estimate.stat[0]
    delta_col = estimate.fit

    sub_num = maths.matrix_exponential_integration.expected_number_subs(
        p0=p0, Q=q, t=tau
    )  # trans tau to sub_num

    p = expm(estimate.calc.Q * tau)

    if model_name in ["ssGN", "GN"]:
        pi = get_stat_pi_via_eigen(p)
    else:
        pi=p0

    if invalid_dlc(
        delta_col,
        limit=numpy.array([pi, pi, pi, pi]),
        to_check=p,
    ):# check if the tau_ estimated render the p into a sympathetic matrix not in limit
        chainsaw_intervals = num_saw_srch(t=lengths(par_name="length", edge=edge_name), q=qc, p0=p0, tau=tau, model_name=model_name)
        return {"edge_name":edge_name, "tau_":tau, "ENS":sub_num, "delta_col":delta_col, "chainsaws":chainsaw_intervals}




def project(qsub_dict:dict, motif_probs:DictArray, tree, model_name, qcsub_dict, lengths) -> ...:
    """
    assumption: the model is identifiable.
    des: it will calc every Q in a lf, and return a table consist of the tau_ point plus the chainsaws
    """
    valid_transits=[]
    for k, v in qsub_dict.items():
        edge_name = k[0]
        p0 = motif_probs[get_node_parent(tree=tree, node=edge_name)]
        if transit:= get_sub_num(q=v,p0=p0,model_name=model_name, edge_name=edge_name, qc=qcsub_dict[k], lengths=lengths):
            valid_transits.append(transit)
    return valid_transits

