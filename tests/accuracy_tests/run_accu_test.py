import os
import shutil
import numpy as np
import warnings


def make_directory(path):
    """
    Prepare a empty folder at the input path
    """
    if os.path.exists(path):
        shutil.rmtree(path)

    os.mkdir(path)


def read_file(file_path):
    """
    INPUT:
    file_path: str, the path of the file to be read
    OUTPUT:
    values: 2d list of str, store the values of the records
    """
    with open(file_path, "r") as file:
        values = [line.strip().split() for line in file]

    return np.array(values)


def write_file(file_path, list_of_data):
    """
    INPUT:
    file_path: str, the path of the file to write
    list_of_data: list of str, the data to be written
    OUTPUT:
    values: 2d list of str, store the values of the records
    """
    with open(file_path, "w") as file:
        for row in list_of_data:
            file.write(" ".join(row) + "\n")


def generate_file_path(file_name):
    """
    returns the corresponding path of the input files
    """
    return os.path.join("tests", "accuracy_tests", "baseData", f"{file_name}")


def standard_input_command(seq):
    """
    generates the general input command
    seq(Boolean) indicates whether the input file is seqfile or genotypes
    """
    if seq:
        command = "AlphaPeel -seqfile " + generate_file_path("sequence.2")
    else:
        command = "AlphaPeel -genotypes " + generate_file_path("genotypes.txt")
    command += (
        " -pedigree "
        + generate_file_path("pedigree.txt")
        + " -nCycles 5"
        + " -maxthreads 6"
    )
    return command


def output_path_command(file_prefix):
    """
    returns the output path command
    """
    path = os.path.join("tests", "accuracy_tests", "outputs", file_prefix)
    command = f"-out {path}"
    return command


def assess_peeling(file_prefix):
    """
    assess the performance of the peeling
    """

    output_path = os.path.join("tests", "accuracy_tests", "outputs", file_prefix)
    true_path = generate_file_path("trueGenotypes.txt")

    new_file = np.loadtxt(output_path)
    trueGenotypes = np.loadtxt(true_path)

    print(" ")
    print("Assessing peeling file: " + file_prefix)

    Marker_accu = [str(get_marker_accu(new_file[:, 1:], trueGenotypes[:, 1:]))]
    for gen in range(5):
        Marker_accu.append(
            str(
                get_marker_accu(
                    new_file[gen * 200 : (gen + 1) * 200],
                    trueGenotypes[gen * 200 : (gen + 1) * 200],
                )
            )
        )

    print("Marker_accuracies", " ".join(Marker_accu))

    Ind_accu = [str(get_ind_accu(new_file[:, 1:], trueGenotypes[:, 1:], None))]
    for gen in range(5):
        Ind_accu.append(str(get_ind_accu(new_file[:, 1:], trueGenotypes[:, 1:], gen)))

    print("Individual_accuracies", " ".join(Ind_accu))


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


def get_ind_accu(output, real, gen=None):
    """
    Get individual accuracy between the output and the real data
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        accus = np.array(
            [np.corrcoef(real[i, :], output[i, :])[0, 1] for i in range(real.shape[0])]
        )
    if type(gen) == int:
        accus = accus[gen * 200 : (gen + 1) * 200]
    return round(np.nanmean(accus), 3)


def test_single(benchmark):
    if os.path.exists("tests/accuracy_tests/accu_report.txt"):
        os.remove("tests/accuracy_tests/accu_report.txt")
    path = os.path.join("tests", "accuracy_tests", "outputs")
    make_directory(path)
    command = (
        standard_input_command(seq=False)
        + " -runType single "
        + output_path_command("peeling.single")
    )
    benchmark(os.system, command)

    assess_peeling("peeling.single.dosages")


def test_single_estmaf(benchmark):
    command = (
        standard_input_command(seq=False)
        + " -runType single"
        + " -estmaf "
        + output_path_command("peeling.single.estmaf")
    )
    benchmark(os.system, command)

    assess_peeling("peeling.single.estmaf.dosages")


def test_multi(benchmark):
    path = os.path.join("tests", "accuracy_tests", "outputs")
    make_directory(path)
    command = (
        standard_input_command(seq=False)
        + " -runType multi "
        + output_path_command("peeling.multi")
    )
    benchmark(os.system, command)

    assess_peeling("peeling.multi.dosages")


def test_multi_estmaf(benchmark):
    command = (
        standard_input_command(seq=False)
        + " -runType multi"
        + " -estmaf "
        + output_path_command("peeling.multi.estmaf")
    )
    benchmark(os.system, command)

    assess_peeling("peeling.multi.estmaf.dosages")


def test_multi_estmaf_estsrrors(benchmark):
    command = (
        standard_input_command(seq=False)
        + " -runType multi"
        + " -esterrors"
        + " -estmaf "
        + output_path_command("peeling.estErrorsAndMaf")
    )
    benchmark(os.system, command)

    assess_peeling("peeling.estErrorsAndMaf.dosages")


def test_multi_seq(benchmark):
    command = (
        standard_input_command(seq=True)
        + " -runType multi "
        + output_path_command("peeling.multi.seq")
    )
    benchmark(os.system, command)

    assess_peeling("peeling.multi.seq.dosages")


def test_hybrid_single(benchmark):
    # generate seg file for first hybrid test
    path = os.path.join("tests", "accuracy_tests", "outputs")

    subset = np.floor(np.linspace(1, 1000, num=200)).astype(dtype=int)
    subset = np.concatenate(([0], subset))

    file_r = os.path.join(path, "peeling.multi.seg")
    seg = read_file(file_r)
    file_w = generate_file_path("seg.subset.txt")
    write_file(file_w, seg[:, subset])

    command = (
        standard_input_command(seq=False)
        + " -runType single"
        + " -mapfile "
        + generate_file_path("map.txt")
        + " -segmapfile "
        + generate_file_path("segmap.txt")
        + " -segfile "
        + generate_file_path("seg.subset.txt")
        + " "
        + output_path_command("peeling.hybrid")
    )
    benchmark(os.system, command)

    assess_peeling("peeling.hybrid.dosages")


def test_hybrid_single_seq(benchmark):
    command = (
        standard_input_command(seq=True)
        + " -runType single"
        + " -mapfile "
        + generate_file_path("map.txt")
        + " -segmapfile "
        + generate_file_path("segmap.txt")
        + " -segfile "
        + generate_file_path("seg.subset.txt")
        + " "
        + output_path_command("peeling.hybrid.seq")
    )
    benchmark(os.system, command)

    assess_peeling("peeling.hybrid.seq.dosages")
