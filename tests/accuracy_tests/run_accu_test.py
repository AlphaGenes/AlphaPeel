import pytest
import os
import shutil
import numpy as np


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
    linesep = os.linesep
    with open(file_path, "w") as file:
        for row in list_of_data:
            file.write(" ".join(row) + linesep)


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
    newAcc = get_gwas_cor(new_file[:, 1:], trueGenotypes[:, 1:])
    print("Accuracy: ", round(newAcc, 3))


def get_gwas_cor(output, real):
    """
    Get correlation between the output and the real data
    """
    # print(f"real.shape[1]: {real.shape[1]}")
    # for i in range(real.shape[1]):
    #     real[:, i]
    #     output[:, i]

    cors = np.array(
        [np.corrcoef(real[:, i], output[:, i])[0, 1] for i in range(real.shape[1])]
    )
    return np.mean(cors)


@pytest.fixture
def commands():
    """
    Return a list of commands
    """
    command_1 = (
        standard_input_command(seq=False)
        + " -runType multi "
        + output_path_command("peeling")
    )
    command_2 = (
        standard_input_command(seq=False)
        + " -runType multi"
        + " -esterrors"
        + " -estmaf "
        + output_path_command("peeling.estErrorsAndMaf")
    )
    # Multi locus with sequence
    command_3 = (
        standard_input_command(seq=True)
        + " -runType multi "
        + output_path_command("peeling.multi.seq")
    )
    # Single locus
    command_4 = (
        standard_input_command(seq=False)
        + " -runType single"
        + " -estmaf "
        + output_path_command("peeling.single")
    )
    command_5 = (
        standard_input_command(seq=False)
        + " -runType single"
        + " -mapfile "
        + generate_file_path("map.txt")
        + " -segmapfile "
        + generate_file_path("segmap.txt")
        + " -segfile "
        + os.path.join("tests", "accuracy_tests", "outputs", "seg.subset.txt")
        + " "
        + output_path_command("peeling.hybrid")
    )
    command_6 = (
        standard_input_command(seq=True)
        + " -runType single"
        + " -mapfile "
        + generate_file_path("map.txt")
        + " -segmapfile "
        + generate_file_path("segmap.txt")
        + " -segfile "
        + os.path.join("tests", "accuracy_tests", "outputs", "seg.subset.txt")
        + " "
        + output_path_command("peeling.hybrid.seq")
    )
    return [command_1, command_2, command_3, command_4, command_5, command_6]


def test_cases(commands):
    """
    Run the tests
    """
    path = os.path.join("tests", "accuracy_tests", "outputs")
    make_directory(path)
    for i in range(4):
        os.system(commands[i])

    subset = np.floor(np.linspace(1, 1000, num=200))
    subset = np.concatenate(([0], subset), dtype=int, casting="unsafe")

    file_r = os.path.join(path, "peeling.seg")
    seg = read_file(file_r)
    file_w = os.path.join(path, "seg.subset.txt")
    write_file(file_w, seg[:, subset])

    os.system(commands[4])
    os.system(commands[5])

    assess_peeling("peeling.dosages")
    assess_peeling("peeling.multi.seq.dosages")
    assess_peeling("peeling.single.dosages")
    assess_peeling("peeling.hybrid.dosages")
    assess_peeling("peeling.hybrid.seq.dosages")
