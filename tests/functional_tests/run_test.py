import os
import shutil
import subprocess
import pytest


def read_and_sort_file(file_name, id_list=None, **kwargs):
    with open(file_name, "r") as file:
        values = [line.strip().split() for line in file]

    if id_list is not None:
        values = [row for row in values if row[0] in id_list]

    values.sort(key=lambda row: row[0])

    return values


@pytest.fixture
def commands_and_paths():
    commands = []
    paths = []

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

    local_variables = locals()

    for variables in local_variables:
        if variables[:8] == "command_":
            commands.append(local_variables[variables])

    for n in range(1, 9):
        paths.append(f"test{n}/outputs")
        if n in [3, 7]:
            paths.append(f"test{n}b/outputs")
            paths.append(f"test{n}c/outputs")

    return commands, paths


def make_direction(path):
    if os.path.exists(path):
        shutil.rmtree(path)

    os.mkdir(path)


def test_cases(commands_and_paths):
    for test_number in [1, 2]:
        command, path = (
            commands_and_paths[0][test_number - 1],
            commands_and_paths[1][test_number - 1],
        )

        make_direction(path)

        subprocess.run(command, shell=True, capture_output=True, text=True)

        output_file_path = path + "/output.called.0.1"
        expected_file_path = path[:-7] + "/trueGenotypes.txt"

        output = read_and_sort_file(output_file_path)
        expected = read_and_sort_file(expected_file_path)

        if test_number == 2:
            for row in expected:
                del row[1]
                del row[4]

        assert output == expected
