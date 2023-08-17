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


@pytest.fixture
def commands():
    """
    Return a list of commands
    """
    command_1 = """
    AlphaPeel -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling \
                        -nCycles 5 \
                        -runType multi \
                        -maxthreads 6
"""
    command_2 = """
    AlphaPeel -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.estErrorsAndMaf \
                        -nCycles 5 \
                        -runType multi \
                        -maxthreads 6 \
                        -esterrors
"""
    # Multi locus with sequence
    command_3 = """
    AlphaPeel -seqfile baseData/sequence.2 \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.multi.seq \
                        -nCycles 5 \
                        -runType multi \
                        -maxthreads 6
"""
    # Single locus
    command_4 = """
    AlphaPeel -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.single \
                        -nCycles 5 \
                        -runType single \
                        -maxthreads 6 \
                        -estmaf
"""
    command_5 = """
    AlphaPeel -genotypes baseData/genotypes.txt \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.hybrid \
                        -nCycles 5 \
                        -runType single \
                        -maxthreads 6 \
                        -mapfile baseData/map.txt\
                        -segmapfile baseData/segmap.txt\
                        -segfile outputs/seg.subset.txt
"""
    command_6 = """
    AlphaPeel -seqfile baseData/sequence.2 \
                        -pedigree baseData/pedigree.txt \
                        -out outputs/peeling.hybrid.seq \
                        -nCycles 5 \
                        -runType single \
                        -maxthreads 6 \
                        -mapfile baseData/map.txt\
                        -segmapfile baseData/segmap.txt\
                        -segfile outputs/seg.subset.txt
"""
    return [command_1, command_2, command_3, command_4, command_5, command_6]


def test_cases(commands):
    """
    Run the tests
    """
    make_directory("outputs")
    for i in range(4):
        os.system(commands[i])

    subset = np.floor(np.linspace(1, 1000, num=200))
    subset = np.concatenate(([0], subset), dtype=int, casting="unsafe")

    seg = read_file("outputs/peeling.seg")
    write_file("outputs/seg.subset.txt", seg[:, subset])

    os.system(commands[4])
    os.system(commands[5])
