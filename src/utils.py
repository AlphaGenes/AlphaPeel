import numpy as np
import os
import argparse


def get_params():
    param_file = os.path.join("tests", "accuracy_tests", "simulation_parameters.txt")
    with open(param_file, "r") as file:
        sim_params = [line.strip().split() for line in file]

    params = {}
    for param_name, param_value in sim_params:
        params[param_name] = float(param_value)

    return params


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
        epilog="The SER and PER_intra calculation follows the definition from: https://www.cell.com/hgg-advances/fulltext/S2666-2477(25)00082-X",
    )

    params = get_params()

    nLociAll = int(params["nLociAll"])
    nInd = int(params["nInd"])

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

    called_file = np.loadtxt(called_path, usecols=np.arange(1, nLociAll + 1))
    true_file = np.loadtxt(true_path, usecols=np.arange(1, nLociAll + 1))

    switch_error_count = 0
    phase_error_count = 0
    uncalled_count = 0
    wrong_homo_count = 0
    true_hetero_count = 0
    homo_count = 0
    hetero_count = 0

    for ind in range(nInd):
        hap_p_new = called_file[ind * 2]
        hap_m_new = called_file[ind * 2 + 1]
        hap_p_true = true_file[ind * 2]
        hap_m_true = true_file[ind * 2 + 1]

        switched = False
        for loci in range(nLociAll):
            if hap_p_true[loci] == hap_m_true[loci]:
                homo_count += 1
            else:
                hetero_count += 1
            if hap_p_new[loci] == 9 or hap_m_new[loci] == 9:
                uncalled_count += 1
                continue
            if hap_p_true[loci] + hap_m_true[loci] == 1:
                if hap_p_new[loci] == hap_m_new[loci]:
                    wrong_homo_count += 1
                    continue
                true_hetero_count += 1
                if hap_p_new[loci] != hap_p_true[loci]:
                    phase_error_count += 1
                if (hap_p_new[loci] != hap_p_true[loci] and switched is False) or (
                    hap_p_new[loci] == hap_p_true[loci] and switched is True
                ):
                    switched = not switched
                    switch_error_count += 1

    print(f"Switch error rate: {switch_error_count / (nInd * (nLociAll - 1))}")
    print(f"Phase error (intra) rate: {phase_error_count / (nInd * nLociAll)}")
    print(f"Uncalled rate: {uncalled_count / (nInd * nLociAll)}")
    print(
        f"Proportion of genotypes wrongly called as homozygote: {wrong_homo_count / (nInd * nLociAll)}"
    )
    print(
        f"Proportion of genotypes correctly called as heterozygote: {true_hetero_count / (nInd * nLociAll)}"
    )
    print(f"Homozygote count in true genotype: {homo_count}")
    print(f"Heterozygote count in true genotype: {hetero_count}")
    print(f"Homo to hetero ratio: {homo_count / hetero_count}")


def main():
    calc_switch_error_rate()


if __name__ == "__main__":
    main()
