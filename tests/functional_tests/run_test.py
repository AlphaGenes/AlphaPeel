import os
import shutil
import subprocess
import pytest


def read_and_sort_file(file_path, id_list=None, **kwargs):
    """
    INPUT:
    file_path: str, the path of the file to be read
    id_list: None or list of int, the ids of the records to be read
    OUTPUT:
    values: 2d list of str, store the sorted values of the (selected) records
    """
    with open(file_path, "r") as file:
        values = [line.strip().split() for line in file]

    if id_list is not None:
        values = [row for row in values if row[0] in id_list]

    values.sort(key=lambda row: row[0])

    return values


@pytest.fixture
def commands_and_paths():
    """
    Return a dictionary of commands and a dictionary of paths
    """
    commands = {}
    paths = {}

    command_1 = "AlphaPeel -genotypes test1/genotypes.txt \
                          -phasefile test1/phasefile.txt \
                          -penetrance test1/penetrance.txt \
                          -seqfile test1/seqfile.txt \
                          -pedigree test1/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 \
                          -esterrors \
                          -out test1/outputs/output"

    command_2 = "AlphaPeel -genotypes test2/genotypes.txt \
                          -phasefile test2/phasefile.txt \
                          -penetrance test2/penetrance.txt \
                          -seqfile test2/seqfile.txt \
                          -pedigree test2/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 \
                          -startsnp 2 \
                          -stopsnp 4 \
                          -out test2/outputs/output"

    command_3 = "AlphaPeel -genotypes test3/genotypes.txt \
                          -phasefile test3/phasefile.txt \
                          -penetrance test3/penetrance.txt \
                          -seqfile test3/seqfile.txt \
                          -pedigree test3/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 0.99\
                          -startsnp 2 \
                          -stopsnp 4 \
                          -binary_call_files \
                          -out test3/outputs/output"

    command_3b = ""

    command_3c = ""

    command_4 = ""

    command_5 = ""

    command_6 = ""

    command_7 = "AlphaPeel -genotypes test7/genotypes.txt \
                          -phasefile test7/phasefile.txt \
                          -penetrance test7/penetrance.txt \
                          -seqfile test7/seqfile.txt \
                          -pedigree test7/pedigree.txt \
                          -runType multi \
                          -esterrors \
                          -calling_threshold .1 \
                          -out test7/outputs/output"

    command_7b = ""

    command_7c = ""

    command_8 = ""

    local_variables = locals()

    for n in range(1, 9):
        test_n = str(n)
        paths[test_n] = f"test{n}/outputs"
        commands[test_n] = local_variables["command_" + test_n]
        if n in [3, 7]:
            paths[test_n + "b"] = f"test{n}b/outputs"
            commands[test_n + "b"] = local_variables["command_" + test_n + "b"]
            paths[test_n + "c"] = f"test{n}c/outputs"
            commands[test_n + "c"] = local_variables["command_" + test_n + "c"]

    return commands, paths


def make_directory(path):
    """
    Prepare a empty folder at the input path
    """
    if os.path.exists(path):
        shutil.rmtree(path)

    os.mkdir(path)


def delete_columns(two_d_list, col_del):
    """
    Delete the columns in col_del of two_d_list
    """
    for n in range(len(col_del)):
        for row in two_d_list:
            del row[col_del[n] - n - 1]


def test_cases(commands_and_paths):
    """
    Run the tests
    """
    tests = ["1", "2", "7"]

    for test_number in tests:
        command, path = (
            commands_and_paths[0][test_number],
            commands_and_paths[1][test_number],
        )

        make_directory(path)

        subprocess.run(command, shell=True, capture_output=True, text=True)

        output_file_path = path + "/output.called.0.1"
        expected_file_path = path[:-7] + "/trueGenotypes.txt"

        output = read_and_sort_file(output_file_path)
        expected = read_and_sort_file(expected_file_path)

        if test_number == "2":
            delete_columns(expected, [2, 6])

        if test_number == "3":
            delete_columns(expected, [2, 6])
            delete_columns(output, [1, 3, 4, 5, 6])

        assert output == expected
