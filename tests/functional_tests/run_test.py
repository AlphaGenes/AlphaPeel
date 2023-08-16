import os
import shutil
import subprocess
import platform
import pytest


def read_file(file_path):
    """
    INPUT:
    file_path: str, the path of the file to be read
    OUTPUT:
    values: 2d list of str, store the values of the records
    """
    with open(file_path, "r") as file:
        values = [line.strip().split() for line in file]

    return values


def read_and_sort_file(file_path, id_list=None, **kwargs):
    """
    INPUT:
    file_path: str, the path of the file to be read
    id_list: None or list of int, the ids of the records to be read
    OUTPUT:
    values: 2d list of str, store the sorted values of the (selected) records
    """
    values = read_file(file_path)

    if id_list is not None:
        values = [row for row in values if row[0] in id_list]

    values.sort(key=lambda row: row[0])

    return values


def generate_file_path(file_type, test_number):
    """
    returns the corresponding path of the input files
    """
    return os.path.join(f"test{test_number}", f"{file_type}.txt")


def standard_input_command(test_number):
    """
    returns the standard input command
    """
    #     command = f"""
    # AlphaPeel -genotypes {generate_file_path('genotypes', test_number)} \
    #     -phasefile {generate_file_path('phasefile', test_number)} \
    #     -penetrance {generate_file_path('penetrance', test_number)} \
    #     -seqfile {generate_file_path('seqfile', test_number)} \
    #     -pedigree {generate_file_path('pedigree', test_number)} \
    # """
    command = (
        "AlphaPeel -genotypes "
        + generate_file_path("genotypes", test_number)
        + " -phasefile "
        + generate_file_path("phasefile", test_number)
        + " -penetrance "
        + generate_file_path("penetrance", test_number)
        + " -seqfile "
        + generate_file_path("seqfile", test_number)
        + " -pedigree "
        + generate_file_path("pedigree", test_number)
    )

    return command


def output_path_command(test_number, file_prefix):
    """
    returns the output path command
    """
    path = os.path.join(f"test{test_number}", "outputs", file_prefix)
    command = f"-out {path}"
    return command


@pytest.fixture
def commands_and_paths():
    """
    Return a dictionary of commands and a dictionary of paths
    """
    commands = {}
    paths = {}

    linesep = os.linesep

    # Test 1: Can we read in unrelated individuals from multiple file formats and
    # output the values to a normal dosage file
    test_number = "1"
    command_1 = (
        standard_input_command(test_number)
        + " -runType multi"
        + " -calling_threshold .1"
        + " -esterrors "
        + output_path_command(test_number, "output")
    )

    # Test 2: Can we read in a subset of values as in Test 1 output them and
    # make sure it's the same chunk?
    test_number = 2

    command_2 = (
        standard_input_command(test_number)
        + " -runType multi"
        + " -calling_threshold .1"
        + " -startsnp 2"
        + " -stopsnp 4 "
        + output_path_command(test_number, "output")
    )

    # Test 3: Can we read in values, call the values and output them as binary?
    command_3 = """
AlphaPeel -genotypes test3/genotypes.txt \
                          -phasefile test3/phasefile.txt \
                          -penetrance test3/penetrance.txt \
                          -seqfile test3/seqfile.txt \
                          -pedigree test3/pedigree.txt \
                          -runType multi \
                          -calling_threshold .1 0.99\
                          -startsnp 2 \
                          -stopsnp 4 \
                          -binary_call_files \
                          -out test3/outputs/output

plink --bfile test3/outputs/output.called.0.1 --real-ref-alleles --recode A --out test3/outputs/output.called.0.1
plink --bfile test3/outputs/output.called.0.99 --real-ref-alleles --recode A --out test3/outputs/output.called.0.99
"""

    # Test 3b: Can we read in a binary file, run the algorithm call the values and check that the output is the same.
    command_3b = """
for nind in 1 2 3 4; do

  AlphaPeel -genotypes test3b/genotypes-$nind.txt \
                            -runType multi \
                            -calling_threshold .1 \
                            -binary_call_files \
                            -out test3b/outputs/output.$nind
  plink --bfile test3b/outputs/output.$nind.called.0.1 --real-ref-alleles --recode A --out test3b/outputs/output.$nind.called.0.1

  AlphaPeel -bfile test3b/outputs/output.$nind.called.0.1 \
                            -runType multi \
                            -calling_threshold .1 \
                            -binary_call_files \
                            -out test3b/outputs/round2.$nind

  plink --bfile test3b/outputs/round2.$nind.called.0.1 --real-ref-alleles --recode A --out test3b/outputs/round2.$nind.called.0.1

done
"""

    # Test 3c: Will the pedigree file be correctly read from the bed file?
    # Create the binary file. Re-run.
    command_3c = """
AlphaPeel -genotypes test3c/genotypes.txt \
                          -runType multi \
                          -calling_threshold .9 \
                          -binary_call_files \
                          -out test3c/outputs/output
plink --bfile test3c/outputs/output.called.0.9 --real-ref-alleles --recode A --out test3c/outputs/output.called.0.9

AlphaPeel -bfile test3c/outputs/output.called.0.9 \
                          -runType multi \
                          -calling_threshold .9 \
                          -out test3c/outputs/noFamNoPedigree

AlphaPeel -bfile test3c/outputs/output.called.0.9 \
                          -pedigree test3c/pedigree.txt \
                          -runType multi \
                          -calling_threshold .9 \
                          -out test3c/outputs/noFamPedigree

cp test3c/fake.fam test3c/outputs/output.called.0.9.fam
AlphaPeel -bfile test3c/outputs/output.called.0.9 \
                          -runType multi \
                          -calling_threshold .9 \
                          -out test3c/outputs/famNoPedigree
"""

    # Test 4a-d: Can we read in values and return them in the correct order.
    # Check id, pedigree, genotypes, sequence, segregation. Also check onlykeyed.
    test_number = "4"
    methods = ["id", "pedigree", "genotypes", "sequence"]
    command_4 = ""
    for method in methods:
        file_prefix = f"output.{method}"
        command_4 += (
            standard_input_command(test_number)
            + " -runType multi"
            + " -calling_threshold .1 "
            + output_path_command(test_number, file_prefix)
            + " -writekey "
            + method
            + linesep
        )
    command_4 += (
        standard_input_command(test_number)
        + " -runType multi"
        + " -calling_threshold .1 "
        + output_path_command(test_number, "output.only")
        + " -writekey sequence"
        + " -onlykeyed"
    )

    # Test 5: Read in an error rate and output genotypes.
    # Should use some silly values, and some not-so-silly values.
    command_5 = ""

    # Test 6: Sex Chromosome
    command_6 = ""

    # Test 7: Check -esterrors just to make sure it runs.
    test_number = "7"
    command_7 = (
        standard_input_command(test_number)
        + " -runType multi"
        + " -esterrors"
        + " -calling_threshold .1 "
        + output_path_command(test_number, "output")
    )

    # Test 7b: Check -estmaf just to make sure it runs.
    test_number = "7b"
    command_7b = (
        standard_input_command(test_number)
        + " -runType multi"
        + " -estmaf"
        + " -calling_threshold .1 "
        + output_path_command(test_number, "output")
    )

    # Test 7c: Check -estmaf just to make sure it runs.
    test_number = "7c"
    command_7c = (
        standard_input_command(test_number)
        + " -runType multi"
        + " -length 1.0"
        + " -calling_threshold .1 "
        + output_path_command(test_number, "output")
    )

    # Test 8: Check to make sure the no_dosages, no_seg, no_params
    # flags work, and the haps file works.
    test_number = "8"
    command_8 = (
        standard_input_command("8")
        + " -runType multi"
        + " -no_dosages "
        + output_path_command(test_number, "no_dosages")
        + linesep
        + standard_input_command("8")
        + " -runType multi"
        + " -no_seg "
        + output_path_command(test_number, "no_seg")
        + linesep
        + standard_input_command("8")
        + " -runType multi"
        + " -no_params "
        + output_path_command(test_number, "no_params")
        + linesep
        + standard_input_command("8")
        + " -runType multi"
        + " -haps "
        + output_path_command(test_number, "haps")
    )

    local_variables = locals()

    for n in range(1, 9):
        test_n = str(n)
        paths[test_n] = os.path.join(f"test{n}", "outputs")
        commands[test_n] = local_variables["command_" + test_n]
        if n in [3, 7]:
            paths[test_n + "b"] = os.path.join(f"test{n}b", "outputs")
            commands[test_n + "b"] = local_variables["command_" + test_n + "b"]
            paths[test_n + "c"] = os.path.join(f"test{n}c", "outputs")
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


def check_for_files(file_prefix):
    """
    Check whether the output files exist
    """
    output = [
        os.path.exists(f"{file_prefix}.dosages"),
        os.path.exists(f"{file_prefix}.seg"),
        os.path.exists(f"{file_prefix}.maf"),
        os.path.exists(f"{file_prefix}.genoError"),
        os.path.exists(f"{file_prefix}.seqError"),
        os.path.exists(f"{file_prefix}.haps"),
    ]
    return output


def test_cases(commands_and_paths):
    """
    Run the tests
    """
    # the numbers of the tests to be run
    tests = ["1", "2", "4", "7", "7b", "7c", "8"]
    os = platform.system()

    for test_number in tests:
        command, path = (
            commands_and_paths[0][test_number],
            commands_and_paths[1][test_number],
        )

        make_directory(path)

        # run the command
        if os == "Windows":
            commands = command.split(os.linesep)
            for one_line_command in commands:
                subprocess.run(command, shell=True, capture_output=True, text=True)
        else:
            subprocess.run(command, shell=True, capture_output=True, text=True)

        if test_number == "4":
            methods = ["id", "pedigree", "genotypes", "sequence"]
            answer = ["genotypes", "penetrance", "genotypes", "seq"]
            for i in range(len(methods)):
                output_file_path = os.path.join(path, f"output.{methods[i]}.called.0.1")
                output = read_file(output_file_path)

                assert len(output) == 4
                assert output[0][0] == answer[i]

            output_file_path = os.path.join(path, "output.only.called.0.1")
            output = read_file(output_file_path)

            assert len(output) == 1
            assert output[0][0] == "seq"

        elif test_number == "8":
            assert check_for_files(os.path.join(path, "no_dosages")) == [
                False,
                True,
                True,
                True,
                True,
                False,
            ]
            assert check_for_files(os.path.join(path, "no_seg")) == [
                True,
                False,
                True,
                True,
                True,
                False,
            ]
            assert check_for_files(os.path.join(path, "no_params")) == [
                True,
                True,
                False,
                False,
                False,
                False,
            ]
            assert check_for_files(os.path.join(path, "haps")) == [
                True,
                True,
                True,
                True,
                True,
                True,
            ]

        else:
            output_file_path = os.path.join(path, "output.called.0.1")
            expected_file_path = os.path.join(path[:-7], "trueGenotypes.txt")

            if test_number == "1":
                assert sorted(os.listdir("test1/outputs/")) == [
                    "output.called.0.1",
                    "output.dosages",
                    "output.genoError",
                    "output.maf",
                    "output.seg",
                    "output.seqError",
                ]
            output = read_and_sort_file(output_file_path)
            expected = read_and_sort_file(expected_file_path)

            if test_number == "2":
                delete_columns(expected, [2, 6])

            if test_number == "3":
                delete_columns(expected, [2, 6])
                delete_columns(output, [1, 3, 4, 5, 6])

            assert output == expected
