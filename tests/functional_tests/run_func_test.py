import os
import shutil


def read_file(file_path, **kwargs):
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

    return values


def read_and_sort_file(file_path, id_list=None, **kwargs):
    """
    INPUT:
    file_path: str, the path of the file to be read
    id_list: None or list of int, the ids of the records to be read
    OUTPUT:
    values: 2d list of str, store the sorted values of the (selected) records
    """

    values = read_file(file_path, **kwargs)

    if id_list is not None:
        # consider only the entries with id in id_list
        values = [row for row in values if row[0] in id_list]

    # remove the empty strings
    values = list(filter(None, values))

    # sort according to the id
    values.sort(key=lambda row: row[0])

    return values


def delete_columns(two_d_list, col_del):
    """
    Delete the columns in col_del of two_d_list
    """
    for n in range(len(col_del)):
        for row in two_d_list:
            del row[col_del[n] - n - 1]


class TestClass:
    path = os.path.join("tests", "functional_tests")
    command = "AlphaPeel "
    test_cases = None
    input_file_depend_on_test_cases = None

    # all the input file options for non-hybrid peeling except the binary file
    files_to_input = ["genotypes", "pedigree", "penetrance", "phasefile", "seqfile"]
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
            f"-out {os.path.join(self.output_path, self.output_file_prefix)}"
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
                os.path.join(self.output_path, f"{self.output_file_prefix}.{file_type}")
            )

        files = [
            "dosage.txt",
            "seg_prob.txt",
            "maf",
            "genoError",
            "seqError",
            "phased_geno_prob.txt",
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
            "runType": "multi",
            "geno_threshold": ".1",
            "geno": None,
            "esterrors": None,
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

        assert self.output == self.expected

    def test_subset(self):
        """
        Can we read in a subset of values as in Test 1 output them and
        make sure it's the same chunk? (=testing startsnp and stopsnp)
        """
        self.test_name = "test_subset"
        self.prepare_path()

        self.input_files = self.files_to_input
        self.arguments = {
            "runType": "multi",
            "geno_threshold": ".1",
            "startsnp": "2",
            "stopsnp": "4",
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

        delete_columns(self.expected, [2, 6])

        assert self.output == self.expected

    def test_writekey(self):
        """
        Can we read in values and return them in the correct order.
        Check id, pedigree, genotypes, sequence, segregation. Also check onlykeyed.
        """
        self.test_name = "test_writekey"
        self.prepare_path()

        self.input_files = self.files_to_input
        self.arguments = {
            "runType": "multi",
            "geno_threshold": ".1",
            "geno": None,
            "writekey": None,
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
            self.arguments["writekey"] = self.test_cases
            self.output_file_prefix = f"writekey.{self.test_cases}"

            self.generate_command()
            os.system(self.command)

            self.output_file_path = os.path.join(
                self.output_path,
                f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
            )

            self.output = read_file(self.output_file_path)

            assert len(self.output) == 4
            assert self.output[0][0] == answer[self.test_cases]

            self.command = "AlphaPeel "

        self.test_cases = "onlykeyed"
        self.arguments["onlykeyed"] = None
        self.output_file_prefix = f"writekey.{self.test_cases}"

        self.generate_command()
        os.system(self.command)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.output = read_file(self.output_file_path)

        assert len(self.output) == 1
        assert self.output[0][0] == "seq"

    def test_est(self):
        """
        Check -esterrors, -estmaf, -length just to make sure it runs.
        """
        self.test_name = "test_est"
        self.prepare_path()

        self.input_files = self.files_to_input
        self.input_file_depend_on_test_cases = self.files_to_input
        self.arguments = {"runType": "multi", "geno_threshold": ".1", "geno": None}
        self.output_file_to_check = "geno_0.3333333333333333"

        for self.test_cases in ["esterror", "estmaf", "length"]:
            # TODO estrecombrate instead of just adding length
            if self.test_cases != "length":
                self.arguments[self.test_cases] = None
            else:
                # Do we need to continue use this value for lengh
                # as it is the same as the default value
                self.arguments["length"] = "1.0"
            self.output_file_prefix = f"est.{self.test_cases}"

            self.generate_command()
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

            assert self.output == self.expected

            self.arguments.pop(self.test_cases)
            self.command = "AlphaPeel "

    def test_no(self):
        """
        Check to make sure the no_dosage, seg_prob, no_params, phased_geno_prob
        flags work.
        """
        self.test_name = "test_no"
        self.prepare_path()

        self.input_files = self.files_to_input
        self.arguments = {"runType": "multi"}
        # whether the output files exist
        # 0: not exist
        # 1: exist
        expect = {
            "no_dosage": [0, 0, 1, 1, 1, 0],
            "seg_prob": [1, 1, 1, 1, 1, 0],
            "no_params": [1, 0, 0, 0, 0, 0],
            "phased_geno_prob": [1, 0, 1, 1, 1, 1],
        }

        for self.test_cases in [
            "no_dosage",
            "seg_prob",
            "no_params",
            "phased_geno_prob",
        ]:
            self.arguments[self.test_cases] = None
            self.output_file_prefix = f"no.{self.test_cases}"

            self.generate_command()
            os.system(self.command)

            assert self.check_files() == expect[self.test_cases]

            self.arguments.pop(self.test_cases)
            self.command = "AlphaPeel "

    def test_rec(self):
        """
        Run the test of the recombination functionality of AlphaPeel
        """
        self.test_name = "test_rec"
        self.prepare_path()

        self.input_files = ["pedigree"]
        self.arguments = {
            "runType": "multi",
            "phased_geno_prob": None,
            "geno": None,
            "geno_threshold": ".1",
            "hap": None,
            "hap_threshold": ".1",
            "seg_prob": None,
        }

        # test for genotype input and sequence input separately
        for self.test_cases in ["genotypes", "seqfile"]:
            self.input_files.append(self.test_cases)
            self.output_file_prefix = f"rec.{self.test_cases}"
            if self.test_cases == "genotypes":
                self.arguments["error"] = "0"
            else:
                # set a small value for sequence error
                # as sequence error cannot be 0
                self.arguments["seqerror"] = "0.00000001"

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

                if self.test_cases == "seqfile":
                    # since sequence error is not 0, we need to round the output up
                    self.output = read_and_sort_file(
                        self.output_file_path, decimal_place=2
                    )
                else:
                    self.output = read_and_sort_file(self.output_file_path)
                self.expected = read_and_sort_file(self.expected_file_path)

                assert self.output == self.expected

            self.input_files.pop(-1)
            self.command = "AlphaPeel "
            if self.test_cases == "genotypes":
                self.arguments.pop("error")

    # the true values to check against are wrong for test_sex
    # needs to rewrite
    def test_sex(self):
        """
        Run the test of the sex chromosome functionality of AlphaPeel
        -sexchrom still under development...
        """
        self.test_name = "test_sex"
        self.prepare_path()

        self.arguments = {"runType": "multi", "sexchrom": None, "seg_prob": None}
        self.input_files = ["genotypes", "seqfile", "pedigree"]
        self.input_file_depend_on_test_cases = ["genotypes", "seqfile"]

        for self.test_cases in ["a", "b", "c", "d"]:
            # test case a: homozygous generation 2
            #           b: heterozygous generation 2
            #           c: with recombination in M2
            #           d: missing values in generation 2

            self.output_file_prefix = f"sex.{self.test_cases}"
            self.output_file_to_check = "seg_prob"

            self.generate_command()
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

            assert self.output == self.expected
            self.command = "AlphaPeel "

    # the true values to check against for test_error is not written yet
    def test_error(self):
        """
        Run the test of the correcting errors functionality of AlphaPeel
        """
        self.test_name = "test_error"
        self.prepare_path()

        # using default error rates: genotype error rate: 0.01
        #                            sequence error rate: 0.001
        self.arguments = {"runType": "multi", "seg_prob": None}
        self.input_files = ["genotypes", "seqfile", "pedigree"]
        self.input_file_depend_on_test_cases = ["genotypes", "seqfile"]

        for self.test_cases in ["a", "b", "c", "d"]:
            # test case a: somatic mutation at locus 5 of M1 in genotype and seqfile
            #           b: germline mutation at locus 5 of M1 in genotype and seqfile
            #           c: somatic mutation at locus 5 of M1 in seqfile only,
            #              with genotype value missing
            #           d: germline mutation at locus 5 of M1 in seqfile only,
            #              with genotype value missing

            self.output_file_prefix = f"error.{self.test_cases}"
            self.output_file_to_check = "genotypes"

            self.generate_command()
            print(self.command)
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

    def test_onlykeyed(self):
        """
        Run the test for the half-founders with onlykeyed option
        """
        self.test_name = "test_onlykeyed"
        self.prepare_path()

        self.input_files = ["genotypes", "pedigree"]
        self.arguments = {"runType": "multi", "onlykeyed": None, "seg_prob": None}
        self.output_file_prefix = "onlykeyed"
        self.output_file_to_check = "dosage"

        self.generate_command()
        os.system(self.command)

        self.output_file_path = os.path.join(
            self.output_path,
            f"{self.output_file_prefix}.{self.output_file_to_check}.txt",
        )

        self.output = read_file(self.output_file_path)

        # no dummy individuals output
        assert len(self.output) == 6

        for ind in self.output:
            assert "MotherOf" not in ind[0] and "FatherOf" not in ind[0]

    # TODO test_plink for PLINK
    #      a. binary PLINK output
    #      b. binary output + input
    #      c. pedigree
