import os
import shutil
import numpy as np
import warnings
import pytest


@pytest.fixture(scope="session")
def sim_path():
    return os.path.join("tests", "accuracy_tests", "sim_for_alphapeel_accu_test")


def prepare_path(output_path):
    """
    Prepare an empty output folder
    """
    if os.path.exists(output_path):
        shutil.rmtree(output_path)

    os.mkdir(output_path)


def generate_output_path(name):
    return os.path.join(
        "tests",
        "accuracy_tests",
        "outputs",
        name,
    )


def generate_command(
    sim_path,
    method,
    est_alt_allele_prob,
    est_geno_error_prob,
    est_seq_error_prob,
    seq_file,
    alt_allele_prob_file,
    update_alt_allele_prob,
    metafounder,
    sex_chrom,
    output_path,
):
    command = "AlphaPeel "
    input_file = ["ped_file"]
    arguments = {
        "method": method,
        "n_cycle": "5",
        "n_thread": "6",
        "geno_threshold": ".1",
        "hap_threshold": ".1",
        "geno": None,
        "hap": None,
        "seg_prob": None,
        "geno_prob": None,
        "phased_geno_prob": None,
    }

    if est_alt_allele_prob:
        arguments["est_alt_allele_prob"] = None
    if est_geno_error_prob and est_seq_error_prob:
        arguments["est_geno_error_prob"] = None
        arguments["est_seq_error_prob"] = None
    if update_alt_allele_prob:
        arguments["update_alt_allele_prob"] = None
    if seq_file:
        input_file.append("seq_file")
    else:
        input_file.append("geno_file")
    if alt_allele_prob_file:
        input_file.append("alt_allele_prob_file")
    if method == "hybrid":
        input_file.append("map_file")
        input_file.append("seg_map_file")
        input_file.append("seg_file")
    if metafounder:
        for file in input_file:
            command += f"-{file} {os.path.join(sim_path, f'metafounder_{file}.txt')} "
    elif sex_chrom:
        for file in input_file:
            command +=  f"-{file} {os.path.join(sim_path, f'X_chr_{file}.txt')} "
        command += "-sex_chrom "
    else:
        for file in input_file:
            command += f"-{file} {os.path.join(sim_path, f'{file}.txt')} "

    for key, value in arguments.items():
        if value is not None:
            command += f"-{key} {value} "
        else:
            command += f"-{key} "

    command += f"-out_file {output_path}{os.sep}"

    return command


def make_directory(path):
    """
    Prepare a empty folder at the input path
    """
    if os.path.exists(path):
        shutil.rmtree(path)

    os.mkdir(path)


def get_marker_accu(output, real):
    """
    Get marker accuracy between the output and the real data
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        accus = np.array(
            [np.corrcoef(real[:, i], output[:, i])[0, 1] for i in range(real.shape[1])]
        )
        return round(np.nanmean(accus), 3)


def get_ind_accu(output, real, nIndPerGen, n_row_per_ind, gen=None):
    """
    Get individual accuracy between the output and the real data
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        accus = np.array(
            [np.corrcoef(real[i, :], output[i, :])[0, 1] for i in range(real.shape[0])]
        )
        if type(gen) == int:
            accus = accus[
                gen
                * (nIndPerGen * n_row_per_ind) : (gen + 1)
                * (nIndPerGen * n_row_per_ind)
            ]
        return round(np.nanmean(accus), 3)


def assess_peeling(
    sim_path, get_params, output_path, name, method, metafounder, sex_chrom
):
    """
    Assess the performance of the peeling
    """
    if metafounder:
        file_to_check = [
            "dosage",
            "geno_prob",
            "phased_geno_prob",
        ]
    else:
        file_to_check = [
            "dosage",
            "geno_0.3333333333333333",
            "hap_0.5",
            "geno_prob",
            "phased_geno_prob",
        ]
    if method == "multi":
        file_to_check.append("seg_prob")

    nGen = int(get_params["nGen"])
    nIndPerGen = int(get_params["nInd"] / nGen)
    nLociAll = int(get_params["nLociAll"])

    print(" ")
    print(f"Test: {name}")

    for file in file_to_check:
        if file in ["dosage", "geno_0.3333333333333333"]:
            n_row_per_ind = 1
        elif file in ["phased_geno_prob", "seg_prob"]:
            n_row_per_ind = 4
        elif file == "geno_prob":
            n_row_per_ind = 3
        elif file == "hap_0.5":
            n_row_per_ind = 2

        file_path = os.path.join(output_path, f".{file}.txt")
        if metafounder:
            true_path = os.path.join(sim_path, f"true-metafounder_{file}.txt")
        elif sex_chrom:
            true_path = os.path.join(sim_path, f"true-X_chr_{file}.txt")
        else:
            true_path = os.path.join(sim_path, f"true-{file}.txt")

        new_file = np.loadtxt(file_path, usecols=np.arange(1, nLociAll + 1))
        true_file = np.loadtxt(true_path, usecols=np.arange(1, nLociAll + 1))

        if metafounder:
            print(f"File: metafounder_{file}")
        elif sex_chrom:
            print(f"File: sex_chrom_{file}")
        else:
            print(f"File: {file}")

        Marker_accu = [str(get_marker_accu(new_file[:, 1:], true_file[:, 1:]))]
        for gen in range(nGen):
            Marker_accu.append(
                str(
                    get_marker_accu(
                        new_file[
                            gen
                            * (nIndPerGen * n_row_per_ind) : (gen + 1)
                            * (nIndPerGen * n_row_per_ind)
                        ],
                        true_file[
                            gen
                            * (nIndPerGen * n_row_per_ind) : (gen + 1)
                            * (nIndPerGen * n_row_per_ind)
                        ],
                    )
                )
            )

        print("Marker_accuracies", " ".join(Marker_accu))

        Ind_accu = [
            str(get_ind_accu(new_file[:, 1:], true_file[:, 1:], nIndPerGen, None))
        ]
        for gen in range(nGen):
            Ind_accu.append(
                str(
                    get_ind_accu(
                        new_file[:, 1:],
                        true_file[:, 1:],
                        nIndPerGen,
                        n_row_per_ind,
                        gen,
                    )
                )
            )

        print("Individual_accuracies", " ".join(Ind_accu))


@pytest.mark.parametrize(
    "method, est_alt_allele_prob, est_geno_error_prob, est_seq_error_prob, seq_file, alt_allele_prob_file, update_alt_allele_prob, metafounder, sex_chrom",
    [
        ("single", None, None, None, None, None, None, None, None),
        ("single", "est_alt_allele_prob", None, None, None, None, None, None, None),
        ("multi", None, None, None, None, None, None, None, None),
        ("multi", "est_alt_allele_prob", None, None, None, None, None, None, None),
        (
            "multi",
            "est_alt_allele_prob",
            "est_geno_error_prob",
            "est_seq_error_prob",
            None,
            None,
            None,
            None,
            None,
        ),
        ("multi", None, None, None, "seq_file", None, None, None, None),
        (
            "multi",
            "est_alt_allele_prob",
            None,
            None,
            "seq_file",
            None,
            None,
            None,
            None,
        ),
        (
            "multi",
            "est_alt_allele_prob",
            "est_geno_error_prob",
            "est_seq_error_prob",
            "seq_file",
            None,
            None,
            None,
            None,
        ),
        ("hybrid", None, None, None, None, None, None, None, None),
        ("hybrid", None, None, None, "seq_file", None, None, None, None),
        (
            "single",
            None,
            None,
            None,
            None,
            "alt_allele_prob_file",
            None,
            "metafounder",
            None,
        ),
        (
            "single",
            "est_alt_allele_prob",
            None,
            None,
            None,
            None,
            None,
            "metafounder",
            None,
        ),
        (
            "single",
            None,
            None,
            None,
            None,
            None,
            "update_alt_allele_prob",
            "metafounder",
            None,
        ),
        (
            "single",
            "est_alt_allele_prob",
            None,
            None,
            None,
            None,
            "update_alt_allele_prob",
            "metafounder",
            None,
        ),
        (
            "single",
            None,
            None,
            None,
            None,
            "alt_allele_prob_file",
            "update_alt_allele_prob",
            "metafounder",
            None,
        ),
        ("single", None, None, None, None, None, None, None, "sex_chrom"),
        ("multi", None, None, None, None, None, None, None, "sex_chrom"),
        ("multi", None, None, None, "seq_file", None, None, None, "sex_chrom"),
        ("hybrid", None, None, None, None, None, None, None, "sex_chrom"),
        ("hybrid", None, None, None, "seq_file", None, None, None, "sex_chrom"),
    ],
)
def test_accu(
    get_params,
    method,
    est_alt_allele_prob,
    est_geno_error_prob,
    est_seq_error_prob,
    seq_file,
    alt_allele_prob_file,
    update_alt_allele_prob,
    metafounder,
    sex_chrom,
    sim_path,
    benchmark,
):
    name = "_".join(
        [
            param
            for param in filter(
                lambda param: True if param else False,
                [
                    method,
                    est_alt_allele_prob,
                    est_geno_error_prob,
                    est_seq_error_prob,
                    seq_file,
                    alt_allele_prob_file,
                    update_alt_allele_prob,
                    metafounder,
                    sex_chrom,
                ],
            )
        ]
    )
    output_path = generate_output_path(name)
    prepare_path(output_path)

    if method == "hybrid":
        # create subset of segregation file
        # make sure run a multi test with the same arguments in advance
        multi_name = "_".join(
            [
                param
                for param in filter(
                    lambda param: True if param else False,
                    [
                        "multi",
                        est_alt_allele_prob,
                        est_geno_error_prob,
                        est_seq_error_prob,
                        seq_file,
                        alt_allele_prob_file,
                        update_alt_allele_prob,
                        metafounder,
                        sex_chrom,
                    ],
                )
            ]
        )
        multi_path = generate_output_path(multi_name)

        nSegMap = int(get_params["nSegMap"])
        nLociAll = int(get_params["nLociAll"])

        subset = np.floor(np.linspace(1, nLociAll, num=nSegMap)).astype(dtype=int)
        subset = np.concatenate(([0], subset))
        seg_path = os.path.join(multi_path, ".seg_prob.txt")
        seg_file_path = os.path.join(sim_path, "seg_file.txt")

        seg = np.loadtxt(seg_path)
        np.savetxt(seg_file_path, seg[:, subset])

    command = generate_command(
        sim_path,
        method,
        est_alt_allele_prob,
        est_geno_error_prob,
        est_seq_error_prob,
        seq_file,
        alt_allele_prob_file,
        update_alt_allele_prob,
        metafounder,
        sex_chrom,
        output_path,
    )

    benchmark(os.system, command)

    assess_peeling(
        sim_path, get_params, output_path, name, method, metafounder, sex_chrom
    )
