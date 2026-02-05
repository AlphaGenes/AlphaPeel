import os
import shutil


def read_file(file_path, test_alt_allele_prob=False, **kwargs):
    """
    INPUT:
    file_path: str, the path of the file to be read
    decimal_place(optional): if provided, round the data with the given
    decimal places
    OUTPUT:
    values: 2d list of str, store the values of the records
    """
    with open(file_path, "r") as file:
        values = [line.strip().split() for line in file]
    if test_alt_allele_prob:
        MF = values[0]
        values.pop(0)

    if "decimal_place" in kwargs.keys():
        # round the data if data and rounding decimal place exist
        values = [
            [line[0]]
            + [round(float(data), kwargs["decimal_place"]) for data in line[1:]]
            if line
            else line
            for line in values
        ]
    else:
        # convert data to float for comparison if data exists
        values = [
            [line[0]] + [float(data) for data in line[1:]] if line else line
            for line in values
        ]
    if test_alt_allele_prob:
        return values, MF
    else:
        return values


def read_and_sort_file(file_path, test_alt_allele_prob=False, id_list=None, **kwargs):
    """
    INPUT:
    file_path: str, the path of the file to be read
    id_list: None or list of int, the ids of the records to be read
    OUTPUT:
    values: 2d list of str, store the sorted values of the (selected) records
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
    """
    Delete the columns in col_del of two_d_list
    """
    for n in range(len(col_del)):
        for row in two_d_list:
            del row[col_del[n] - n - 1]


def read_geno_hap(file_path):
    """
    Read the geno/hap file and return a dictionary
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
    """
    Compare the output file with the true file
    the error tolerance is mismatch genotype <= total_error
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
    path = os.path.join("tests", "functional_tests")
    test_cases = None
    input_file_depend_on_test_cases = None

    # all the input file options for non-hybrid peeling except the binary file
    files_to_input = ["geno_file", "ped_file", "penetrance", "hap_file", "seq_file"]
    # all the output files except the binary file and the parameter files
    files_to_check = [
        "hap_0.5",
        "geno_0.3333333333333333",
        "dosage",
        "phased_geno_prob",
        "seg_prob",
    ]

    def mk_output_dir(self):
        """
        Prepare a empty folder at the input path
        """
        if os.path.exists(self.output_path):
            shutil.rmtree(self.output_path)

        os.mkdir(self.output_path)

    def generate_command(self):
        """
        generate the command for the test
        """
        self.command = "AlphaPeel "

        for file in self.input_files:
            if (
                (self.test_cases is not None)
                and (self.input_file_depend_on_test_cases is not None)
                and (file in self.input_file_depend_on_test_cases)
            ):
                self.command += f"-{file} {os.path.join(self.path, f'{file}-{self.test_cases}.txt')} "
            else:
                self.command += f"-{file} {os.path.join(self.path, f'{file}.txt')} "

        for key, value in self.arguments.items():
            if value is not None:
                self.command += f"-{key} {value} "
            else:
                self.command += f"-{key} "

        self.command += (
            f"-out_file {os.path.join(self.output_path, self.output_file_prefix)}"
        )

    def prepare_path(self):
        """
        Initialize the paths for the test
        """
        self.path = os.path.join(self.path, self.test_name)
        self.output_path = os.path.join(self.path, "outputs")
        self.mk_output_dir()

    def check_files(self):
        """
        Check the existence of the output files
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
        ]
        return [check(file) for file in files]

    def test_files(self):
        """
        Can we read in unrelated individuals from multiple file formats and
        output the values to a normal dosage file
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
        self.output_file_to_check = "geno_0.3333333333333333"

        self.generate_command()
        os.system(self.command)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )
        self.expected_file_path = os.path.join(
            self.path, f"true-{self.output_file_to_check}.txt"
        )

        self.output = read_and_sort_file(self.output_file_path)
        self.expected = read_and_sort_file(self.expected_file_path)

        # Produced dosage file correctly where multiple different files are inputted of unrelated individuals:
        # Geno_file, hap_file, ped_file, penetrance, and seq_file
        assert self.output == self.expected

    def test_subset(self):
        """
        Can we read in a subset of values as in Test 1 output them and
        make sure it's the same chunk? (=testing start_snp and stop_snp)
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
        self.output_file_to_check = "geno_0.3333333333333333"

        self.generate_command()
        os.system(self.command)

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
        """
        Can we read in values and return them in the correct order.
        Check id, pedigree, genotypes, sequence, segregation. Also check out_id_only.
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

        self.output_file_to_check = "geno_0.3333333333333333"

        for self.test_cases in methods:
            self.arguments["out_id_order"] = self.test_cases
            self.output_file_prefix = f"out_id_order.{self.test_cases}"

            self.generate_command()
            os.system(self.command)

            self.output_file_path = os.path.join(
                self.output_path,
                f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
            )

            self.output = read_file(self.output_file_path)
            # Total four individuals across the five inputted files
            assert len(self.output) == 4
            # Check the outputted order under different commands: id, pedigree, genotypes, sequence.
            assert self.output[0][0] == answer[self.test_cases]

            self.command = "AlphaPeel "

        self.test_cases = "out_id_only"
        self.arguments["out_id_only"] = None
        self.output_file_prefix = f"out_id_order.{self.test_cases}"

        self.generate_command()
        os.system(self.command)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.output = read_file(self.output_file_path)
        # one observation as one individual in seq file
        assert len(self.output) == 1
        # First ID equal to "seq" (id used in sequence file)
        assert self.output[0][0] == "seq"

    def test_est(self):
        """
        Check -est_geno_error_prob, -est_seq_error_prob, -est_alt_allele_prob, -rec_length just to make sure it runs.
        """
        self.test_name = "test_est"
        self.prepare_path()

        self.input_files = self.files_to_input
        self.input_file_depend_on_test_cases = self.files_to_input
        self.arguments = {"method": "multi", "geno_threshold": ".1", "geno": None}
        self.output_file_to_check = "geno_0.3333333333333333"

        for self.test_cases in [
            "est_geno_error_prob",
            "est_seq_error_prob",
            "est_alt_allele_prob",
            "rec_length",
        ]:
            # TODO estrecombrate instead of just adding length
            if self.test_cases != "rec_length":
                self.arguments[self.test_cases] = None
            else:
                # Do we need to continue use this value for lengh
                # as it is the same as the default value
                self.arguments["rec_length"] = "1.0"
            self.output_file_prefix = f"est.{self.test_cases}"

            self.generate_command()
            print(self.command)
            os.system(self.command)

            self.output_file_path = os.path.join(
                self.output_path,
                f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
            )
            self.expected_file_path = os.path.join(
                self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
            )

            self.output = read_and_sort_file(self.output_file_path)
            self.expected = read_and_sort_file(self.expected_file_path)
            # Checking AlphaPeel runs so compares outputted genotypes, haplotypes, seq, and penetrance with expected.
            assert self.output == self.expected

            self.arguments.pop(self.test_cases)
            self.command = "AlphaPeel "

    def test_no(self):
        """
        Check to make sure the no_dosage, seg_prob, no_param, phased_geno_prob
        flags work.
        """
        self.test_name = "test_no"
        self.prepare_path()

        self.input_files = self.files_to_input
        self.arguments = {"method": "multi"}
        # whether the output files exist
        # 0: not exist
        # 1: exist
        expect = {
            "no_dosage": [0, 0, 0, 1, 1, 0],
            "seg_prob": [1, 1, 0, 1, 1, 0],
            "alt_allele_prob": [1, 0, 1, 1, 1, 0],
            "no_param": [1, 0, 0, 0, 0, 0],
            "phased_geno_prob": [1, 0, 0, 1, 1, 1],
        }

        for self.test_cases in [
            "no_dosage",
            "seg_prob",
            "alt_allele_prob",
            "no_param",
            "phased_geno_prob",
        ]:
            self.arguments[self.test_cases] = None
            self.output_file_prefix = f"no.{self.test_cases}"

            self.generate_command()
            os.system(self.command)
            # When requested through commands, test the presents of file outputs:
            # no_dosage, output files: alt_allele_prob, geno_error_prob, seg_error_prob
            # seg_prob, output files: dosage, seg_prob, alt_allele_prob, geno_error_prob, seg_error_prob
            # no_param, output file: dosage
            # phased_geno_prob, output files: dosage, alt_allele_prob, geno_error_prob, seg_error_prob, phased_geno_prob
            assert self.check_files() == expect[self.test_cases]

            self.arguments.pop(self.test_cases)
            self.command = "AlphaPeel "

    def test_rec(self):
        """
        Run the test of the recombination functionality of AlphaPeel
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

            self.generate_command()
            os.system(self.command)

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
                # Check outputted dosage, phased, and seg to expected files
                assert self.output == self.expected

            self.input_files.pop(-1)
            self.command = "AlphaPeel "
            if self.test_cases == "geno_file":
                self.arguments.pop("geno_error_prob")

    # the true values to check against are wrong for test_sex
    # needs to rewrite
    def test_sex(self):
        """
        Run the test of the sex chromosome functionality of AlphaPeel
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
            self.output_file_to_check = ["geno_0.3333333333333333", "hap_0.5"]

            self.generate_command()
            os.system(self.command)

            for check in self.output_file_to_check:
                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{check}.txt",
                )
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_prefix}.{check}.txt"
                )
                # Compares outputted genotype files to expected.
                if check == "geno_0.3333333333333333":
                    compare_geno_hap(
                        self.output_file_path, self.expected_file_path, total_error=2
                    )
                else:
                    compare_geno_hap(
                        self.output_file_path, self.expected_file_path, total_error=5
                    )

            self.command = "AlphaPeel "

    # the true values to check against for test_error is not written yet
    def test_error(self):
        """
        Run the test of the correcting errors functionality of AlphaPeel
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
            # self.output_file_to_check = "genotypes"

            self.generate_command()
            os.system(self.command)

            # self.output_file_path = os.path.join(
            #     self.output_path,
            #     f"{self.output_file_prefix}.{self.output_file_to_check}.txt"
            #     )
            # self.expected_file_path = os.path.join(
            #     self.path,
            #     f"true-{self.output_file_to_check}-{self.test_cases}.txt"
            #     )

            # self.output = read_and_sort_file(self.output_file_path)
            # self.expected = read_and_sort_file(self.expected_file_path)

            # assert self.output == self.expected
            self.command = "AlphaPeel "

    def test_out_id_only(self):
        """
        Run the test for the half-founders with out_id_only option
        """
        self.test_name = "test_out_id_only"
        self.prepare_path()

        self.input_files = ["geno_file", "ped_file"]
        self.arguments = {"method": "multi", "out_id_only": None, "seg_prob": None}
        self.output_file_prefix = "out_id_only"
        self.output_file_to_check = "dosage"

        self.generate_command()
        os.system(self.command)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.output = read_file(self.output_file_path)

        # Number of observations in dosage file is 6 as no dummy individuals in output.
        assert len(self.output) == 6

        for ind in self.output:
            assert "MotherOf" not in ind[0] and "FatherOf" not in ind[0]

    def test_alt_allele_prob(self):
        """
        Run the test for alt_allele_prob and est_alt_allele_prob with multiple metafounders
        """
        self.test_name = "test_alt_allele_prob"
        self.prepare_path()

        self.input_files = ["geno_file", "ped_file"]
        self.input_file_depend_on_test_cases = self.input_files
        self.arguments = {"method": "multi", "out_id_only": None}

        for self.test_cases in [
            "default",
            "alt_allele_prob_file_single",
            "alt_allele_prob_file_multiple",
            "update_alt_allele_prob_single",
            "update_alt_allele_prob_multiple",
            "both",
            "incorrect_pedigree",
            "default_metafounder",
            "main_metafounder",
            "incorrect_main_metafounder",
            "incorrect_metafounder_in_file",
            "missing_metafounder_in_file",
            "extra_metafounder_in_file",
            "metafounder_order_in_output",
        ]:
            # test case default: Test the default values of the alternative allele frequency
            #                    without any input or estimation with multiple metafounders
            #           alt_allele_prob_file_single: Test the input option alt_allele_prob_file
            #                                   for a single metafounder
            #           alt_allele_prob_file_multiple: Test the input option alt_allele_prob_file
            #                                     for multiple metafounders
            #           update_alt_allele_prob_single: Test the option alt_allele_prob_method
            #                                       for a single metafounder
            #           update_alt_allele_prob_multiple: Test the option alt_allele_prob_method
            #                                         for multiple metafounders
            #           both: Test the case when both options are used,
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
            #           missing_metafounder_in_file: Test case when a metafounder is present in the pedigree but missing in the input alternative allele probability file.
            #           extra_metafounder_in_file: Test case when an additional metafounder is present in the input alternative allele probability file.
            #           metafounder_order_in_output: Check the order of metafounders in the output is numerical when over 10 metafounders.

            self.output_file_prefix = f"alt_allele_prob.{self.test_cases}"

            if self.test_cases == "default":
                self.output_file_to_check = "alt_allele_prob"
                self.arguments[
                    "alt_allele_prob"
                ] = None  # To output the alt_allele_prob
                self.generate_command()
                os.system(self.command)

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
                # self.input_files.append("alt_allele_prob_file")
                self.input_file_depend_on_test_cases.append("alt_allele_prob_file")

                self.output_file_to_check = "dosage"

                self.generate_command()
                os.system(self.command)

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
                self.generate_command()
                os.system(self.command)

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

                # self.input_files.pop(-1)
                self.input_file_depend_on_test_cases.pop(-1)

            elif self.test_cases == "update_alt_allele_prob_single":
                self.arguments["update_alt_allele_prob"] = None
                self.output_file_to_check = "alt_allele_prob"

                self.generate_command()
                os.system(self.command)

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

            elif self.test_cases == "update_alt_allele_prob_multiple":
                self.generate_command()
                os.system(self.command)

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

            elif self.test_cases == "both":
                self.input_file_depend_on_test_cases.append("alt_allele_prob_file")
                self.expected_file_path = os.path.join(
                    self.path, f"true-{self.output_file_to_check}-{self.test_cases}.txt"
                )

                self.generate_command()
                os.system(self.command)

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

                self.arguments.pop("update_alt_allele_prob")

            elif self.test_cases == "incorrect_pedigree":
                self.generate_command()
                exit_code = os.system(self.command)

                # check if error message is in the output
                assert exit_code in [512, 2]

                self.input_file_depend_on_test_cases.pop(-1)

            elif self.test_cases == "default_metafounder":
                self.generate_command()
                os.system(self.command)

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
                self.generate_command()
                os.system(self.command)

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
                self.generate_command()
                exit_code = os.system(self.command)

                # check if error message is in the output
                assert exit_code in [512, 2]

                self.arguments.pop("main_metafounder")

            elif self.test_cases == "incorrect_metafounder_in_file":
                self.input_files.append("alt_allele_prob_file")
                self.input_file_depend_on_test_cases.append("alt_allele_prob_file")

                self.generate_command()
                exit_code = os.system(self.command)

                # check if error message is in the output
                assert exit_code in [512, 2]

            elif self.test_cases == "missing_metafounder_in_file":
                self.generate_command()
                os.system(self.command)

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
                self.generate_command()
                os.system(self.command)

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

                self.input_files.pop(-1)
                self.input_file_depend_on_test_cases.pop(-1)

            elif self.test_cases == "metafounder_order_in_output":
                self.generate_command()
                os.system(self.command)

                self.output_file_path = os.path.join(
                    self.output_path,
                    f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
                )

                self.output, MF = read_and_sort_file(
                    self.output_file_path, test_alt_allele_prob=True
                )

                # Check that the metafounders are in the correct order in the output
                print(MF)
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

            self.command = "AlphaPeel "

    def test_pheno(self):
        """
        Testing of the phenotype functionality in AlphaPeel
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
                self.generate_command()
                os.system(self.command)

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

                self.command = "AlphaPeel "

            elif self.test_cases == "pheno_probs_with_penetrance":
                # This will print phenotype probabilities
                self.input_file_depend_on_test_cases.append("pheno_penetrance_file")

                self.generate_command()
                os.system(self.command)

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
                # Compares the outputted pheno_probs file to the expected based on inputted pheno_penetrance_file.
                assert self.output == self.expected

                self.command = "AlphaPeel "

            elif self.test_cases == "pheno_file_with_penetrance":
                # This will update the dosage file from pheno data and print phenotype probabilities
                self.input_file_depend_on_test_cases.append("pheno_file")

                self.generate_command()
                os.system(self.command)

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
                # Compares the outputted pheno_probs file to the expected based on inputted pheno_penetrance_file.
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
                # Compares the outputted dosage file to the expected based on inputted pheno_penetrance_file.
                assert self.output == self.expected

                self.command = "AlphaPeel "

            elif self.test_cases == "repeat_pheno_record":
                # This will update the dosage and pheno_prob file
                self.generate_command()
                os.system(self.command)

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
                # Compares the outputted pheno_probs file to the expected based on inputted pheno_penetrance_file.
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
                # Compares the outputted dosage file to the expected based on inputted pheno_penetrance_file.
                assert self.output == self.expected

                self.command = "AlphaPeel "

            elif self.test_cases == "multi_pheno_state":
                # This will update the dosage and pheno_prob file
                self.generate_command()
                print(self.command)
                os.system(self.command)

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
                # Compares the outputted pheno_probs file to the expected based on inputted pheno_penetrance_file.
                assert self.output == self.expected

                self.command = "AlphaPeel "

            elif self.test_cases == "pheno_file_with_multi_loci_geno_file":
                # This will flag an error and exit the program (at the moment)
                self.generate_command()
                exit_code = os.system(self.command)
                # check if error message is in the output
                assert exit_code in [256, 512, 2]

                self.input_file_depend_on_test_cases.pop(-2)

                self.command = "AlphaPeel "

            elif self.test_cases == "pheno_file_only":
                # This will flag an error and exit the program
                self.generate_command()
                exit_code = os.system(self.command)
                # check if error message is in the output
                assert exit_code in [256, 512, 2]

                self.command = "AlphaPeel "

    def test_map_input(self):
        """
        Run the test for the input map file
        """
        self.test_name = "test_map_input"
        self.prepare_path()

        self.arguments = {"method": "multi"}
        self.output_file_to_check = "dosage"

        # without map file input
        self.input_files = ["geno_file", "ped_file"]
        self.output_file_prefix = "map_input.no_map_file"

        self.generate_command()
        os.system(self.command)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.first_output = read_and_sort_file(self.output_file_path)

        # with map file input
        self.input_files.append("map_file")
        self.output_file_prefix = "map_input.with_map_file"

        self.generate_command()
        os.system(self.command)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.second_output = read_and_sort_file(self.output_file_path)

        # the two outputs should match
        assert self.first_output == self.second_output

    def test_prev_bug(self):
        """
        Run the test for the previous bug described in tinyhouse#165
        """
        self.test_name = "test_prev_bug"
        self.prepare_path()

        self.input_files = ["ped_file", "geno_file"]
        self.arguments = {"method": "multi"}

        self.output_file_prefix = "prev_bug"
        self.output_file_to_check = "dosage"

        self.generate_command()
        os.system(self.command)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.output = read_and_sort_file(self.output_file_path)

        assert round(float(self.output[1][1])) == 2

    # TODO test_plink for PLINK
    #      a. binary PLINK output
    #      b. binary output + input
    #      c. ped_file
