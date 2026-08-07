"""Module for calculating switch error rate and phase error rate of haplotype files."""
import argparse
import os

import numpy as np

from src.accuracy_assessment import get_hap_switch_error_metrics
from src.accuracy_core import get_params


def calc_switch_error_rate():
    """
    Run this function at root directory and input the path to the
    assessed haplotype file.
    The SER and PER_intra calculation follows the definition from:
    https://www.cell.com/hgg-advances/fulltext/S2666-2477(25)00082-X
    """

    parser = argparse.ArgumentParser(
        prog="calc_switch_error_rate",
        description="Calculate switch error rate, phase error rate",
        epilog="The SER and PER_intra calculation follows the definition from: "
        "https://www.cell.com/hgg-advances/fulltext/S2666-2477(25)00082-X",
    )

    params = get_params()

    n_loci_all = int(params["nLociAll"])
    n_ind = int(params["nInd"])

    parser.add_argument(
        "-true_path",
        type=str,
        required=False,
        default=os.path.join(
            "tests", "accuracy_tests", "sim_for_alphapeel_accu_test", "true-hap_0.5.txt"
        ),
        help="Enter the path of the true haplotype file, default is the simulation path",
    )
    parser.add_argument(
        "-called_path",
        type=str,
        required=True,
        help="Enter the path of the assessed haplotype file",
    )

    args = parser.parse_args()
    true_path = args.true_path
    called_path = args.called_path

    called_file = np.loadtxt(called_path, usecols=np.arange(1, n_loci_all + 1))
    true_file = np.loadtxt(true_path, usecols=np.arange(1, n_loci_all + 1))

    metrics = dict(
        get_hap_switch_error_metrics(called_file, true_file, n_ind, n_loci_all)
    )

    print(f"Switch error rate: {metrics['switch_error_rate']}")
    print(f"Phase error (intra) rate: {metrics['phase_error_rate']}")
    print(f"Uncalled rate: {metrics['uncalled_rate']}")
    print(
        "Proportion of genotypes wrongly called as homozygote: "
        f"{metrics['wrong_homozygote_rate']}"
    )
    print(
        "Proportion of genotypes correctly called as heterozygote: "
        f"{metrics['correct_heterozygote_rate']}"
    )
    print(f"Homozygote count in true genotype: {metrics['homozygote_count']}")
    print(f"Heterozygote count in true genotype: {metrics['heterozygote_count']}")
    print(f"Homo to hetero ratio: {metrics['homo_to_hetero_ratio']}")
