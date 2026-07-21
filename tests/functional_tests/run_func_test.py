import os
import shutil
import src.tinypeel.tinypeel as tinypeel
import pytest


def read_file(file_path, test_alt_allele_prob=False, **kwargs):
    """Read a whitespace-separated file and return its records.

    :param file_path: Path to the input file.
    :type file_path: str
    :param test_alt_allele_prob: If True, read the metafounders and
    the corresponding alt allele probabilities. Defaults to False.
    :type test_alt_allele_prob: bool
    :param decimal_place: When provided, numeric fields are rounded to this
        many decimal places. Passed as a keyword argument via ``kwargs``.
    :type decimal_place: int or None
    :return: 2-D list of records where each row is a list starting with an id
        followed by numeric values converted to ``float`` (or rounded if
        ``decimal_place`` is given). When ``test_alt_allele_prob`` is True,
        returns a tuple ``(values, metafounder_header)``.
    :rtype: list or tuple
    """
    with open(file_path, "r") as file:
        values = [line.strip().split() for line in file]
    if test_alt_allele_prob:
        MF = values[0]
        values.pop(0)

    decimal_place = kwargs.get("decimal_place")
    if test_alt_allele_prob:
        # Alt allele probability files have only numeric probability columns.
        if decimal_place is None:
            values = [
                [float(data) for data in line] if line else line for line in values
            ]
        else:
            values = [
                [round(float(data), decimal_place) for data in line] if line else line
                for line in values
            ]
    elif decimal_place is not None:
        # Round the data columns while preserving the leading id column.
        values = [
            [line[0]] + [round(float(data), decimal_place) for data in line[1:]]
            if line
            else line
            for line in values
        ]
    else:
        # Convert data columns to float for comparison while preserving ids.
        values = [
            [line[0]] + [float(data) for data in line[1:]] if line else line
            for line in values
        ]
    if test_alt_allele_prob:
        return values, MF
    else:
        return values


def read_and_sort_file(file_path, test_alt_allele_prob=False, id_list=None, **kwargs):
    """Read a file, optionally filter by id list, and return sorted records.

    This wraps :func:`read_file` then filters out empty rows, optionally
    restricts to records with ids in ``id_list``, and sorts by id.

    :param file_path: Path to the input file.
    :type file_path: str
    :param test_alt_allele_prob: If True, read the metafounders and
    the corresponding alt allele probabilities. Defaults to False.
    :type test_alt_allele_prob: bool
    :param id_list: Optional iterable of ids to select. If provided, only rows
        whose first element (id) is in this list are kept.
    :type id_list: list or None
    :param decimal_place: See :func:`read_file` for rounding behavior.
    :type decimal_place: int or None
    :return: Sorted list of records, or (values, metafounder_header) when
        ``test_alt_allele_prob`` is True.
    :rtype: list or tuple
    """

    if test_alt_allele_prob:
        values, MF = read_file(file_path, test_alt_allele_prob, **kwargs)
    else:
        values = read_file(file_path, test_alt_allele_prob, **kwargs)

    if id_list is not None:
        # consider only the entries with id in id_list
        values = [row for row in values if row[0] in id_list]

    # remove the empty strings
    values = list(filter(None, values))

    # sort according to the id
    values.sort(key=lambda row: row[0])

    if test_alt_allele_prob:
        return values, MF
    else:
        return values


def delete_columns(two_d_list, col_del):
    """Delete columns from a 2-D list in-place.

    Columns are removed in the order specified by ``col_del``. Each entry in
    ``col_del`` is a 1-based column index; the function adjusts indices while
    deleting so that later deletions refer to the updated list shape.

    :param two_d_list: A list of rows (each a mutable sequence) to modify.
    :type two_d_list: list
    :param col_del: Sequence of 1-based column indices to remove.
    :type col_del: list[int]
    :return: None. The operation modifies ``two_d_list`` in-place.
    :rtype: None
    """
    for n in range(len(col_del)):
        for row in two_d_list:
            del row[col_del[n] - n - 1]


def read_geno_hap(file_path):
    """Read a genotype/haplotype file into a dictionary keyed by id.

    The input file is expected to be whitespace-separated with each non-empty
    line starting with an identifier followed by values. Multiple lines with
    the same identifier are grouped into a list of value-rows.

    :param file_path: Path to the genotype/haplotype file.
    :type file_path: str
    :return: Dictionary mapping identifier -> list of rows (each a list of
        string values following the identifier).
    :rtype: dict[str, list[list[str]]]
    """
    dic_file = {}
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip().split()
            if line == []:
                continue
            if line[0] not in dic_file:
                dic_file[line[0]] = [line[1:]]
            else:
                dic_file[line[0]].append(line[1:])
    return dic_file


def compare_geno_hap(output, true, total_error=2):
    """Compare two genotype/haplotype files and raise on too many mismatches.

    The function reads both files via :func:`read_geno_hap`, asserts that the
    set of ids match, and counts mismatches between corresponding entries. If
    the number of mismatches exceeds ``total_error`` a :class:`ValueError` is
    raised.

    :param output: Path to the output file produced by the code under test.
    :type output: str
    :param true: Path to the expected (true) file to compare against.
    :type true: str
    :param total_error: Maximum allowed number of mismatches before failing.
    :type total_error: int
    :return: None. Raises on failure.
    :rtype: None
    """
    outputs = read_geno_hap(output)
    trues = read_geno_hap(true)

    trues_id = sorted(trues.keys())
    outputs_id = sorted(outputs.keys())
    # check the ids are the same
    assert trues_id == outputs_id
    number_error = 0
    # check the content
    for i in trues_id:
        value_output = outputs[i]
        value_true = trues[i]
        for k in range(len(value_output)):
            rows_value_output = value_output[k]
            rows_value_true = value_true[k]
            number = len(rows_value_output)
            for j in range(number):
                # check the number of mismatches
                if rows_value_output[j] != rows_value_true[j]:
                    number_error += 1
                    print(f"the {i} {j}th genotype/haplotype{k} is different")
                if number_error > total_error:
                    raise ValueError("the number of error is larger than 2")


class TestClass:
    """Functional test runner helpers for AlphaPeel unit tests.

    This class provides helper methods and test cases that invoke the
    `tinypeel` module with different argument combinations and verify
    the produced output files against expected results located in the
    ``tests/functional_tests`` directory.
    """

    path = os.path.join("tests", "functional_tests")
    test_cases = None
    input_file_depend_on_test_cases = None

    # all the input file options for non-hybrid peeling except the binary file
    files_to_input = [
        "geno_file",
        "ped_file",
        "phased_geno_prob_file",
        "hap_file",
        "seq_file",
    ]
    # all the output files except the binary file and the parameter files
    files_to_check = [
        "hap_0.5",
        "geno_0.333",
        "dosage",
        "phased_geno_prob",
        "seg_prob",
    ]

    def mk_output_dir(self):
        """Create an empty output directory for the current test.

        If the directory already exists it is removed first.

        :return: None
        :rtype: None
        """
        if os.path.exists(self.output_path):
            shutil.rmtree(self.output_path)

        os.mkdir(self.output_path)

    def generate_arguments(self):
        """Build the ``argv`` list to pass to :mod:`tinypeel` for the test.

        The method inspects ``self.input_files`` and ``self.arguments`` and
        constructs the corresponding command-line arguments stored in
        ``self.argv``.

        :return: None. Result is stored in ``self.argv``.
        :rtype: None
        """
        self.argv = []

        for file in self.input_files:
            if (
                (self.test_cases is not None)
                and (self.input_file_depend_on_test_cases is not None)
                and (file in self.input_file_depend_on_test_cases)
            ):
                self.argv.append(f"-{file}")
                self.argv.append(
                    os.path.join(self.path, f"{file}-{self.test_cases}.txt")
                )
            else:
                self.argv.append(f"-{file}")
                self.argv.append(os.path.join(self.path, f"{file}.txt"))

        for key, value in self.arguments.items():
            self.argv.append(f"-{key}")
            if value is not None:
                self.argv.append(value)

        self.argv.append("-out_file")
        self.argv.append(os.path.join(self.output_path, self.output_file_prefix))

    def prepare_path(self):
        """Prepare test file paths and create the output directory.

        Sets ``self.path`` to the directory for the current test case and
        ``self.output_path`` to its ``outputs`` subdirectory, then ensures the
        output directory exists and is empty.

        :return: None
        :rtype: None
        """
        self.path = os.path.join(self.path, self.test_name)
        self.output_path = os.path.join(self.path, "outputs")
        self.mk_output_dir()

    def check_files(self):
        """Check whether the expected output files were created.

        :return: List of booleans indicating presence for each of the files in
            the fixed check list.
        :rtype: list[bool]
        """

        def check(file_type):
            return os.path.exists(
                os.path.join(
                    self.output_path, f"{self.output_file_prefix}.{file_type}.txt"
                )
            )

        files = [
            "dosage",
            "seg_prob",
            "alt_allele_prob",
            "geno_error_prob",
            "seq_error_prob",
            "phased_geno_prob",
            "rec_prob",
        ]
        return [check(file) for file in files]

    def test_files(self):
        """Functional test: various input formats produce the expected genotype.

        This test provides multiple unrelated input files (geno, hap, pedigree,
        phased probabilities and sequence) and checks that the generated
        genotype output matches the expected file.

        :return: None. Uses :mod:`pytest` assertions.
        """
        self.test_name = "test_files"
        self.prepare_path()

        self.input_files = self.files_to_input
        self.arguments = {
            "method": "multi",
            "geno_threshold": ".1",
            "geno": None,
            "est_geno_error_prob": None,
            "est_seq_error_prob": None,
            "seg_prob": None,
        }
        self.output_file_prefix = "files"
        self.output_file_to_check = "geno_0.333"

        self.generate_arguments()
        with pytest.warns(
            UserWarning,
            match="External phased genotype probability file included, but est_geno_error_prob flag used. The two options are incompatible. est_geno_error_prob set to false.",
        ):
            with pytest.warns(
                UserWarning,
                match="External phased genotype probability file included, but est_seq_error_prob flag used. The two options are incompatible. est_seq_error_prob set to false.",
            ):
                tinypeel.main(argv=self.argv)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )
        self.expected_file_path = os.path.join(
            self.path, f"true-{self.output_file_to_check}.txt"
        )

        self.output = read_and_sort_file(self.output_file_path)
        self.expected = read_and_sort_file(self.expected_file_path)

        # Produced genotype file correctly where multiple different files are inputted of unrelated individuals:
        # Geno_file, hap_file, ped_file, phased_geno_prob_file, and seq_file
        assert self.output == self.expected

    def test_subset(self):
        """Functional test: check if subsetting of SNPs using ``start_snp``/``stop_snp``
        flags works correctly.

        Ensures that selecting a subset of loci produces output equal to the
        corresponding chunk of the full expected output.
        """
        self.test_name = "test_subset"
        self.prepare_path()

        self.input_files = self.files_to_input
        self.arguments = {
            "method": "multi",
            "geno_threshold": ".1",
            "start_snp": "2",
            "stop_snp": "4",
            "geno": None,
            "seg_prob": None,
        }
        self.output_file_prefix = "subset"
        self.output_file_to_check = "geno_0.333"

        self.generate_arguments()
        tinypeel.main(argv=self.argv)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )
        self.expected_file_path = os.path.join(
            self.path, f"true-{self.output_file_to_check}.txt"
        )

        self.output = read_and_sort_file(self.output_file_path)
        self.expected = read_and_sort_file(self.expected_file_path)
        # Remove the first and last column inline with start_snp (2) and stop_snp (4) command
        delete_columns(self.expected, [2, 6])
        # Test start_snp and stop_snp commands across different input files
        # Compares the outputted genotype dosage file with the expected genotype dosage.
        assert self.output == self.expected

    def test_out_id_order(self):
        """Functional test: verify output ordering behavior controlled by flags.

        Tests the ``out_id_order`` option and verifies
        that the output ordering matches the expected order for several modes.
        """
        self.test_name = "test_out_id_order"
        self.prepare_path()

        self.input_files = self.files_to_input
        self.arguments = {
            "method": "multi",
            "geno_threshold": ".1",
            "geno": None,
            "out_id_order": None,
            "seg_prob": None,
        }

        methods = ["id", "pedigree", "genotypes", "sequence"]
        answer = {
            "id": "genotypes",
            "pedigree": "penetrance",
            "genotypes": "genotypes",
            "sequence": "seq",
        }

        self.output_file_to_check = "geno_0.333"

        for self.test_cases in methods:
            self.arguments["out_id_order"] = self.test_cases
            self.output_file_prefix = f"out_id_order.{self.test_cases}"

            self.generate_arguments()
            tinypeel.main(argv=self.argv)

            self.output_file_path = os.path.join(
                self.output_path,
                f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
            )

            self.output = read_file(self.output_file_path)
            # Total four individuals across the five inputted files
            assert len(self.output) == 4
            # Check the outputted order under different commands: id, pedigree, genotypes, sequence.
            assert self.output[0][0] == answer[self.test_cases]

        self.test_cases = "out_id_only"
        self.arguments["out_id_only"] = None
        self.output_file_prefix = f"out_id_order.{self.test_cases}"

        self.generate_arguments()
        tinypeel.main(argv=self.argv)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.output = read_file(self.output_file_path)
        # one observation as one individual in seq file
        assert len(self.output) == 1
        # First ID equal to "seq" (id used in sequence file)
        assert self.output[0][0] == "seq"

    def test_out_id_only(self):
        """Functional test: behavior of ``out_id_only`` option for half founders.

        Confirms that only non-dummy individuals are present in the output and
        that identifiers for parents (``MotherOf``/``FatherOf``) are excluded.
        """
        self.test_name = "test_out_id_only"
        self.prepare_path()

        self.input_files = ["geno_file", "ped_file"]
        self.arguments = {"method": "multi", "out_id_only": None, "seg_prob": None}
        self.output_file_prefix = "out_id_only"
        self.output_file_to_check = "dosage"

        self.generate_arguments()
        tinypeel.main(argv=self.argv)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.output = read_file(self.output_file_path)

        # Number of observations in dosage file is 6 as no dummy individuals in output.
        assert len(self.output) == 6

        for ind in self.output:
            assert "MotherOf" not in ind[0] and "FatherOf" not in ind[0]

    def test_est(self):
        """Functional test: exercise various estimation-related flags.

        Verifies that running with ``est_geno_error_prob``, ``est_seq_error_prob``,
        ``est_start_alt_allele_prob`` and ``rec_length`` executes and produces
        expected outputs where applicable.
        """
        self.test_name = "test_est"
        self.prepare_path()

        self.input_files = self.files_to_input.copy()
        self.input_files.remove("phased_geno_prob_file")
        self.input_file_depend_on_test_cases = self.input_files
        self.arguments = {"method": "multi", "geno_threshold": ".1", "geno": None}
        self.output_file_to_check = "geno_0.333"

        for self.test_cases in [
            "est_geno_error_prob",
            "est_seq_error_prob",
            "rec_length",
        ]:
            # TODO estrecombrate instead of just adding length
            if self.test_cases != "rec_length":
                self.arguments[self.test_cases] = None
            else:
                # Do we need to continue use this value for lengh
                # as it is the same as the default value
                self.arguments["rec_length"] = "1.0"
                self.input_files.append("phased_geno_prob_file")
            self.output_file_prefix = f"est.{self.test_cases}"

            self.generate_arguments()
            tinypeel.main(argv=self.argv)

            self.output_file_path = os.path.join(
                self.output_path,
                f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
            )
            self.expected_file_path = os.path.join(
                self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
            )

            self.output = read_and_sort_file(self.output_file_path)
            self.expected = read_and_sort_file(self.expected_file_path)
            # Checking AlphaPeel runs so compares outputted genotypes with expected.
            assert self.output == self.expected

            self.arguments.pop(self.test_cases)

        # test est_start_alt_allele_prob functionality
        self.test_cases = "est_start_alt_allele_prob"
        self.input_files = ["geno_file", "ped_file"]
        self.input_file_depend_on_test_cases = self.input_files
        self.arguments["est_start_alt_allele_prob"] = None
        self.output_file_prefix = f"est.{self.test_cases}"

        self.generate_arguments()
        tinypeel.main(argv=self.argv)

        self.output_file_to_check = "alt_allele_prob"
        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )
        self.expected_file_path = os.path.join(
            self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
        )

        self.output, MF = read_and_sort_file(
            self.output_file_path, test_alt_allele_prob=True, decimal_place=1
        )
        self.expected, MF = read_and_sort_file(
            self.expected_file_path, test_alt_allele_prob=True, decimal_place=1
        )
        assert self.output == self.expected

        self.output_file_to_check = "dosage"

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )
        self.expected_file_path = os.path.join(
            self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
        )

        self.output = read_and_sort_file(self.output_file_path, decimal_place=1)
        self.expected = read_and_sort_file(self.expected_file_path, decimal_place=1)

        # Check the outputted dosage file with expected when est_start_alt_allele_prob is used as starting point for estimation.
        assert self.output == self.expected

        self.arguments.pop(self.test_cases)

    def test_no(self):
        """Functional test: verify which output files are produced for flags.

        Iterates over a set of flags and asserts that the presence/absence of
        each expected output file matches a predefined expectation matrix.

        The set of flags include:
        - ``no_dosage``
        - ``seg_prob``
        - ``est_geno_error_prob``
        - ``est_seq_error_prob``
        - ``est_start_alt_allele_prob``
        - ``est_alt_allele_prob``
        - ``alt_allele_prob``
        - ``phased_geno_prob``
        - ``rec_prob``
        """
        self.test_name = "test_no"
        self.prepare_path()

        self.input_files = self.files_to_input.copy()
        self.input_files.remove("phased_geno_prob_file")
        self.arguments = {"method": "multi"}
        # whether the output files exist
        # 0: not exist
        # 1: exist
        expect = {
            "no_dosage": [0, 0, 0, 0, 0, 0, 0],
            "seg_prob": [1, 1, 0, 0, 0, 0, 0],
            "est_geno_error_prob": [1, 0, 0, 1, 0, 0, 0],
            "est_seq_error_prob": [1, 0, 0, 0, 1, 0, 0],
            "est_start_alt_allele_prob": [1, 0, 1, 0, 0, 0, 0],
            "est_alt_allele_prob": [1, 0, 1, 0, 0, 0, 0],
            "alt_allele_prob": [1, 0, 1, 0, 0, 0, 0],
            "phased_geno_prob": [1, 0, 0, 0, 0, 1, 0],
            "rec_prob": [1, 0, 0, 0, 0, 0, 1],
        }

        for self.test_cases in [
            "no_dosage",
            "seg_prob",
            "est_geno_error_prob",
            "est_seq_error_prob",
            "est_start_alt_allele_prob",
            "est_alt_allele_prob",
            "alt_allele_prob",
            "phased_geno_prob",
            "rec_prob",
        ]:
            self.arguments[self.test_cases] = None
            self.output_file_prefix = f"no.{self.test_cases}"

            self.generate_arguments()
            tinypeel.main(argv=self.argv)
            # When requested through commands, test the presents of file outputs:
            # no_dosage, output files: no output
            # seg_prob, output files: dosage, seg_prob
            # est_geno_error_prob, output files: dosage, geno_error_prob
            # est_seq_error_prob, output files: dosage, seq_error_prob
            # est_start_alt_allele_prob, output files: dosage, alt_allele_prob
            # est_alt_allele_prob, output files: dosage, alt_allele_prob
            # alt_allele_prob, output files: dosage, alt_allele_prob
            # phased_geno_prob, output files: dosage, phased_geno_prob
            # rec_prob, output files: dosage, rec_prob
            assert self.check_files() == expect[self.test_cases]

            self.arguments.pop(self.test_cases)

    def test_rec(self):
        """Functional test: recombination-related processing.

        This test suppose to check the recombination probability, but it
        currently is not implemented correctly.
        """
        self.test_name = "test_rec"
        self.prepare_path()

        self.input_files = ["ped_file"]
        self.arguments = {
            "method": "multi",
            "phased_geno_prob": None,
            "geno": None,
            "geno_threshold": ".1",
            "hap": None,
            "hap_threshold": ".1",
            "seg_prob": None,
        }

        # test for genotype input and sequence input separately
        for self.test_cases in ["geno_file", "seq_file"]:
            self.input_files.append(self.test_cases)
            self.output_file_prefix = f"rec.{self.test_cases}"
            if self.test_cases == "geno_file":
                self.arguments["geno_error_prob"] = "0"
            else:
                # set a small value for sequence error
                # as sequence error cannot be 0
                self.arguments["seq_error_prob"] = "0.00000001"

            self.generate_arguments()
            tinypeel.main(argv=self.argv)

            # bug in writeCalledPhase()
            # skip called phase file for now
            self.files_to_check.pop(0)

            for self.output_file_to_check in self.files_to_check:
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}.txt"
                )

                if self.test_cases == "seq_file":
                    # since sequence error is not 0, we need to round the output up
                    self.output = read_and_sort_file(
                        self.output_file_path, decimal_place=2
                    )
                else:
                    self.output = read_and_sort_file(self.output_file_path)
                self.expected = read_and_sort_file(self.expected_file_path)
                # Check outputted genotype, dosage, phased, and seg to expected files
                assert self.output == self.expected

            self.input_files.pop(-1)
            if self.test_cases == "geno_file":
                self.arguments.pop("geno_error_prob")

    def test_sex(self):
        """Functional test: sex chromosome handling.

        Verifies correct behaviour on X-chromosome data across several
        recombination/missing-data scenarios. Uses :func:`compare_geno_hap` to
        compare outputs with expected files and allows a larger mismatch
        tolerance for phased/haplotypes when appropriate.
        The true values to check against are wrong for test_sex,
        needs to rewrite.
        """
        self.test_name = "test_sex"
        self.prepare_path()

        self.arguments = {
            "method": "multi",
            "x_chr": None,
            "hap": None,
            "geno": None,
        }
        self.input_files = ["geno_file", "ped_file"]
        self.input_file_depend_on_test_cases = self.input_files

        for self.test_cases in [
            "no_recom",
            "no_recom_missing",
            "with_recom",
            "with_recom_missing",
        ]:
            self.output_file_prefix = f"sex.{self.test_cases}"
            self.output_file_to_check = ["geno_0.333", "hap_0.5"]

            self.generate_arguments()
            tinypeel.main(argv=self.argv)

            for check in self.output_file_to_check:
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_prefix}.{check}.txt"
                )
                # Compares outputted genotype files to expected.
                if check == "geno_0.333":
                    compare_geno_hap(
                        self.output_file_path, self.expected_file_path, total_error=2
                    )
                else:
                    compare_geno_hap(
                        self.output_file_path, self.expected_file_path, total_error=5
                    )

    def test_error(self):
        """Functional test: error-correction scenarios.

        The true values to check against for test_error is not written yet.
        """
        self.test_name = "test_error"
        self.prepare_path()

        # using default error rates: genotype error rate: 0.001
        #                            sequence error rate: 0.0001
        self.arguments = {"method": "multi", "seg_prob": None}
        self.input_files = ["geno_file", "seq_file", "ped_file"]
        self.input_file_depend_on_test_cases = ["geno_file", "seq_file"]

        for self.test_cases in ["a", "b", "c", "d"]:
            # test case a: somatic mutation at locus 5 of M1 in geno_file and seq_file
            #           b: germline mutation at locus 5 of M1 in geno_file and seq_file
            #           c: somatic mutation at locus 5 of M1 in seq_file only,
            #              with genotype value missing
            #           d: germline mutation at locus 5 of M1 in seq_file only,
            #              with genotype value missing

            self.output_file_prefix = f"error.{self.test_cases}"

            self.generate_arguments()
            tinypeel.main(argv=self.argv)

    def test_alt_allele_prob(self):
        """Functional test: alternative allele probability handling and metafounders.

        Tests many sub-cases related to input files, estimation, multiple
        metafounders, and error conditions. Each sub-case compares outputs to
        expected files or asserts that the correct exception is raised.
        """
        self.test_name = "test_alt_allele_prob"
        self.prepare_path()

        self.input_files = ["geno_file", "ped_file"]
        self.input_file_depend_on_test_cases = self.input_files
        self.arguments = {"method": "multi", "out_id_only": None, "out_digits": "2"}

        for self.test_cases in [
            "default",
            "alt_allele_prob_file_single",
            "alt_allele_prob_file_multiple",
            "multiple_metafounder_individual",
            "est_alt_allele_prob_single",
            "est_alt_allele_prob_multiple",
            "est_alt_allele_prob_multiple_individual",
            "one_metafounder_individual",
            "both",
            "incorrect_pedigree",
            "default_metafounder",
            "main_metafounder",
            "incorrect_main_metafounder",
            "incorrect_metafounder_in_file",
            "missing_metafounder_in_file",
            "extra_metafounder_in_file",
            "alt_allele_prob_missing",
            "alt_allele_prob_missing_9",
            "alt_allele_prob_missing_metafounder",
            "metafounder_order_in_output",
        ]:
            # test case default: Test the default values of the alternative allele frequency
            #                    without any input or estimation with multiple metafounders
            #           alt_allele_prob_file_single: Test the input option alt_allele_prob_file
            #                                   for a single metafounder
            #           alt_allele_prob_file_multiple: Test the input option alt_allele_prob_file
            #                                     for multiple metafounders
            #           multiple_metafounder_individual: Test case when an individual has multiple metafounders
            #                                       assigned and whether the average alternative allele
            #                                       probabilities are being used correctly
            #           est_alt_allele_prob_single: Test the option est_alt_allele_prob
            #                                       for a single metafounder
            #           est_alt_allele_prob_multiple: Test the option est_alt_allele_prob
            #                                         for multiple metafounders
            #           est_alt_allele_prob_multiple_individual: Test case when an individual has multiple metafounders
            #                                       assigned and whether the average alternative allele
            #                                       probabilities are being updated correctly with est_alt_allele_prob
            #           one_metafounder_individual: Test case when an individual has one metafounder. This will trigger an error
            #           both: Test the case when both alt_allele_prob_file and est_alt_allele_prob options are used,
            #                 whether the inputted alternative allele probabilities are used as
            #                 a starting point for alternative allele probabilities estimation
            #           incorrect_pedigree: Test case when a metafounder is written incorrectly as
            #                               not a founder in the pedigree file
            #           default_metafounder: Test case when 0 is being used as parents and
            #                                no main metafounder is being provided as input,
            #                                test whether 0 would be replaced by the default MF_1
            #           main_metafounder: Test if the input option main_metafounder is working
            #                             i.e the user defines the default metafounder where 0 is used.
            #           incorrect_main_metafounder: Test case when the input main_metafounder does not start with MF_,
            #                                       whether an error would be raised
            #           incorrect_metafounder_in_file: Test case when the names of input metafounders
            #                                          in the input alternative allele probability file do not start with MF_,
            #                                          whether an error would be raised
            #           missing_metafounder_in_file: Test case when a metafounder is present in the pedigree but missing in the
            #                                       input alternative allele probability file.
            #           extra_metafounder_in_file: Test case when an additional metafounder is present in the input alternative
            #                                       allele probability file.
            #           alt_allele_prob_missing: Test case when a value is missing in
            #                                     the inputted alternative allele probability file.
            #           alt_allele_prob_missing_9: Test case when a value is missing (with a missing value of 9) in the inputted alternative allele probability file.
            #           alt_allele_prob_missing_metafounder: Test case when a metafounder is listed but empty in the
            #                                       inputted alternative allele probability file.
            #           metafounder_order_in_output: Check the order of metafounders in the output is numerical when over 10 metafounders.

            self.output_file_prefix = f"alt_allele_prob.{self.test_cases}"

            if self.test_cases == "default":
                self.output_file_to_check = "alt_allele_prob"
                self.arguments[
                    "alt_allele_prob"
                ] = None  # To output the alt_allele_prob
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True, decimal_place=1
                )
                self.expected, MF = read_and_sort_file(
                    self.expected_file_path, test_alt_allele_prob=True, decimal_place=1
                )
                # Each metafounder has alt_allele_prob of 0.5 per marker
                assert self.output == self.expected

            elif self.test_cases == "alt_allele_prob_file_single":
                self.input_file_depend_on_test_cases.append("alt_allele_prob_file")

                self.output_file_to_check = "dosage"

                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )

                self.output = read_and_sort_file(self.output_file_path, decimal_place=1)
                self.expected = read_and_sort_file(self.expected_file_path)

                # Compares the outputted dosage file to the expected based on inputted alt_allele_prob file.
                assert self.output == self.expected

            elif self.test_cases == "alt_allele_prob_file_multiple":
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )

                self.output = read_and_sort_file(self.output_file_path, decimal_place=1)
                self.expected = read_and_sort_file(self.expected_file_path)
                # Compares the outputted dosage file to the expected based on inputted alt_allele_prob file.
                assert self.output == self.expected
            elif self.test_cases == "multiple_metafounder_individual":
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )

                self.output = read_and_sort_file(self.output_file_path, decimal_place=1)
                self.expected = read_and_sort_file(self.expected_file_path)
                # Compares the outputted dosage file to the expected based on inputted alt_allele_prob file
                assert self.output == self.expected
                # self.input_files.pop(-1)
                self.input_file_depend_on_test_cases.pop(-1)

            elif self.test_cases == "est_alt_allele_prob_single":
                self.arguments["est_alt_allele_prob"] = None
                self.output_file_to_check = "alt_allele_prob"

                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )
                self.expected, MF = read_and_sort_file(
                    self.expected_file_path, test_alt_allele_prob=True
                )
                # Compares alt_allele_prob output with expected when estimated by AlphaPeel for one metafounder
                assert self.output == self.expected

            elif self.test_cases == "est_alt_allele_prob_multiple":
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )
                self.expected, MF = read_and_sort_file(
                    self.expected_file_path, test_alt_allele_prob=True
                )
                # Compares alt_allele_prob output with expected when estimated by AlphaPeel for multiple metafounders
                assert self.output == self.expected
            elif self.test_cases == "est_alt_allele_prob_multiple_individual":
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )
                self.expected, MF = read_and_sort_file(
                    self.expected_file_path, test_alt_allele_prob=True
                )
                # Compares alt_allele_prob output with expected when estimated by AlphaPeel for multiple metafounders per individual
                assert self.output == self.expected
            elif self.test_cases == "one_metafounder_individual":
                self.generate_arguments()

                with pytest.raises(
                    ValueError,
                    match="Both parents must be metafounders if one is a metafounder. For individual D0 the parents were F0 and MF_2.\nConsider using a dummy individual for the metafounder.",
                ):
                    tinypeel.main(argv=self.argv)

            elif self.test_cases == "both":
                self.input_file_depend_on_test_cases.append("alt_allele_prob_file")
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )

                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )
                self.expected, MF = read_and_sort_file(
                    self.expected_file_path, test_alt_allele_prob=True
                )

                # check if the estimated alt_allele_prob is 0.5
                assert self.output == self.expected

                self.arguments.pop("est_alt_allele_prob")

            elif self.test_cases == "incorrect_pedigree":
                self.generate_arguments()

                with pytest.raises(
                    ValueError,
                    match="Individual MF_1 uses the prefix 'MF_' which is reserved for metafounders and cannot be used for an individual's id.",
                ):
                    tinypeel.main(argv=self.argv)

                self.input_file_depend_on_test_cases.pop(-1)

            elif self.test_cases == "default_metafounder":
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_to_check = "alt_allele_prob"
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )

                # check if there is only one metafounder
                assert len(MF) == 1
                # check if the name of the metafounder is MF_1
                assert MF[0] == "MF_1"

            elif self.test_cases == "main_metafounder":
                self.arguments["main_metafounder"] = "MF_test"
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )

                # check if there is only one metafounder
                assert len(MF) == 1
                # check if the name of the metafounder is MF_test
                assert MF[0] == "MF_test"

            elif self.test_cases == "incorrect_main_metafounder":
                self.arguments["main_metafounder"] = "test"
                self.generate_arguments()

                with pytest.raises(
                    ValueError, match="The main_metafounder must start with MF_."
                ):
                    tinypeel.main(argv=self.argv)

                self.arguments.pop("main_metafounder")

            elif self.test_cases == "incorrect_metafounder_in_file":
                self.input_files.append("alt_allele_prob_file")
                self.input_file_depend_on_test_cases.append("alt_allele_prob_file")
                self.generate_arguments()

                with pytest.raises(
                    ValueError,
                    match="Incorrect number of locus rows in the `alt_allele_prob_file`. Expected 5 rows but found 6.",
                ):
                    tinypeel.main(argv=self.argv)

            elif self.test_cases == "missing_metafounder_in_file":
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path,
                    f"true-{self.output_file_to_check}-{self.test_cases}.txt",
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )

                self.expected, MF = read_and_sort_file(
                    self.expected_file_path, test_alt_allele_prob=True
                )

                # Check that the default is assigned to any metafounders present in the pedigree but not the alt_allele_prob_file
                assert self.output == self.expected

            elif self.test_cases == "extra_metafounder_in_file":
                self.generate_arguments()
                with pytest.warns(
                    UserWarning,
                    match="MF_3 is not in the pedigree. The alternative allele probability for MF_3 has been ignored.",
                ):
                    tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path,
                    f"true-{self.output_file_to_check}-{self.test_cases}.txt",
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )

                self.expected, MF = read_and_sort_file(
                    self.expected_file_path, test_alt_allele_prob=True
                )
                # Check that the extra is removed from the alt_allele_prob output
                assert self.output == self.expected

            elif self.test_cases == "alt_allele_prob_missing":
                self.generate_arguments()
                with pytest.raises(
                    ValueError,
                    match="In the `alt_allele_prob_file`, locus row 3 has 1 values but header has 2 metafounders. If the alternative allele probability is unknown, please use default of 0.5.",
                ):
                    tinypeel.main(argv=self.argv)

            elif self.test_cases == "alt_allele_prob_missing_9":
                self.generate_arguments()
                with pytest.raises(
                    ValueError,
                    match=r"Invalid value 9.0 for alternative allele probability for metafounder MF_2 at locus 2. \nValues must be between 0 and 1. Set to 0.5 \(default\) if unknown.",
                ):
                    tinypeel.main(argv=self.argv)

            elif self.test_cases == "alt_allele_prob_missing_metafounder":
                self.generate_arguments()
                with pytest.raises(
                    ValueError,
                    match="In the `alt_allele_prob_file`, locus row 1 has 1 values but header has 2 metafounders. If the alternative allele probability is unknown, please use default of 0.5.",
                ):
                    tinypeel.main(argv=self.argv)

            elif self.test_cases == "metafounder_order_in_output":
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )

                # Check that the metafounders are in the correct order in the output
                assert MF == [
                    "MF_1",
                    "MF_2",
                    "MF_3",
                    "MF_4",
                    "MF_5",
                    "MF_6",
                    "MF_7",
                    "MF_8",
                    "MF_9",
                    "MF_10",
                    "MF_11",
                    "MF_12",
                ]

    def test_pheno(self):
        """Functional test: phenotype and penetrance handling.

        Verifies that phenotype probabilities are produced only when the
        appropriate penetrance file is provided and that dosage updates from
        phenotype information behave as expected.
        """
        self.test_name = "test_pheno"
        self.prepare_path()

        self.input_files = ["geno_file", "ped_file"]
        self.input_file_depend_on_test_cases = self.input_files
        self.arguments = {
            "method": "single",
            "out_id_only": None,
        }

        for self.test_cases in [
            "pheno_probs_no_penetrance",
            "pheno_probs_with_penetrance",
            "pheno_file_with_penetrance",
            "repeat_pheno_record",
            "multi_pheno_state",
            "pheno_file_with_multi_loci_geno_file",
            "pheno_file_only",
        ]:
            self.output_file_prefix = f"pheno.{self.test_cases}"

            if self.test_cases == "pheno_probs_no_penetrance":
                # This will give a warning and not print phenotype probabilities
                self.arguments["pheno_prob"] = None
                self.generate_arguments()

                with pytest.warns(
                    UserWarning,
                    match="Phenotype probabilities are not available. Please provide a penetrance file with -pheno_penetrance_prob_file. -pheno_prob will be ignored.",
                ):
                    tinypeel.main(argv=self.argv)

                self.output_file_to_check = "pheno_prob"
                # Check the pheno_prob file does not exist

                test = os.path.exists(
                    os.path.join(
                        self.output_path,
                        f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                    )
                )

                expect = False
                assert test == expect

            elif self.test_cases == "pheno_probs_with_penetrance":
                # This will print phenotype probabilities
                self.input_file_depend_on_test_cases.append(
                    "pheno_penetrance_prob_file"
                )

                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_to_check = "pheno_prob"
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )
                self.output = read_and_sort_file(self.output_file_path)
                self.expected = read_and_sort_file(self.expected_file_path)
                # Compares the outputted pheno_probs file to the expected based on inputted pheno_penetrance_prob_file.
                assert self.output == self.expected

            elif self.test_cases == "pheno_file_with_penetrance":
                # This will update the dosage file from pheno data and print phenotype probabilities
                self.input_file_depend_on_test_cases.append("pheno_file")

                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_to_check = "pheno_prob"
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )
                self.output = read_and_sort_file(self.output_file_path)
                self.expected = read_and_sort_file(self.expected_file_path)
                # Compares the outputted pheno_probs file to the expected based on inputted pheno_penetrance_prob_file.
                assert self.output == self.expected

                self.output_file_to_check = "dosage"
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )
                self.output = read_and_sort_file(self.output_file_path)
                self.expected = read_and_sort_file(self.expected_file_path)
                # Compares the outputted dosage file to the expected based on inputted pheno_penetrance_prob_file.
                assert self.output == self.expected

            elif self.test_cases == "repeat_pheno_record":
                # This will update the dosage and pheno_prob file
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_to_check = "pheno_prob"
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )
                self.output = read_and_sort_file(self.output_file_path)
                self.expected = read_and_sort_file(self.expected_file_path)
                # Compares the outputted pheno_probs file to the expected based on inputted pheno_penetrance_prob_file.
                assert self.output == self.expected

                self.output_file_to_check = "dosage"
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )
                self.output = read_and_sort_file(self.output_file_path)
                self.expected = read_and_sort_file(self.expected_file_path)
                # Compares the outputted dosage file to the expected based on inputted pheno_penetrance_prob_file.
                assert self.output == self.expected

            elif self.test_cases == "multi_pheno_state":
                # This will update the dosage and pheno_prob file
                self.generate_arguments()
                tinypeel.main(argv=self.argv)

                self.output_file_to_check = "pheno_prob"
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )
                self.output = read_and_sort_file(self.output_file_path)
                self.expected = read_and_sort_file(self.expected_file_path)
                # Compares the outputted pheno_probs file to the expected based on inputted pheno_penetrance_prob_file.
                assert self.output == self.expected

            elif self.test_cases == "pheno_file_with_multi_loci_geno_file":
                # This will flag an error and exit the program (at the moment)
                self.generate_arguments()
                with pytest.raises(
                    ValueError,
                    match="Currently phenotype information can only be used with a single locus genotype input. Please either remove the pheno_file or use a single locus genotype input.",
                ):
                    tinypeel.main(argv=self.argv)

                self.input_file_depend_on_test_cases.pop(-2)

            elif self.test_cases == "pheno_file_only":
                # This will flag an error and exit the program
                self.generate_arguments()
                with pytest.raises(
                    ValueError,
                    match="To use phenotype information, please provide a phenotype penetrance via '-pheno_penetrance_file'.",
                ):
                    tinypeel.main(argv=self.argv)

    def test_map_input(self):
        """Functional test: behaviour when providing a genetic map file.

        Compares outputs produced when a specific map file is supplied to outputs from
        the same run without the map file to ensure equivalence for a given
        start/stop SNP range.
        """
        self.test_name = "test_map_input"
        self.prepare_path()

        self.arguments = {"method": "multi", "start_snp": "2", "stop_snp": "5"}
        self.output_file_to_check = "dosage"

        # without map file input
        self.input_files = ["geno_file", "ped_file"]
        self.output_file_prefix = "map_input.no_map_file"

        self.generate_arguments()
        tinypeel.main(argv=self.argv)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.first_output = read_and_sort_file(self.output_file_path)

        # with map file input
        self.input_files.append("map_file")
        self.output_file_prefix = "map_input.with_map_file"

        self.generate_arguments()
        tinypeel.main(argv=self.argv)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.second_output = read_and_sort_file(self.output_file_path)

        # the two outputs should match
        assert self.first_output == self.second_output

    def test_founder_status(self):
        """Functional test: regression test for founder status bug (tinyhouse#165).

        Ensures the previous bug is fixed by asserting a specific numeric value
        in the produced dosage output.
        """
        self.test_name = "test_founder_status"
        self.prepare_path()

        self.input_files = ["ped_file", "geno_file"]
        self.arguments = {"method": "multi"}

        self.output_file_prefix = "founder_status"
        self.output_file_to_check = "dosage"

        self.generate_arguments()
        tinypeel.main(argv=self.argv)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.output = read_and_sort_file(self.output_file_path)

        assert round(float(self.output[1][1])) == 2
