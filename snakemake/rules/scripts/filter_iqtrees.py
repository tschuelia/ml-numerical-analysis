import regex

# define some regex stuff
blanks = r"\s+"  # matches >=1  subsequent whitespace characters
sign = r"[-+]?"  # contains either a '-' or a '+' symbol or none of both
float_re = r"\d+(?:\.\d+)?(?:[e][-+]?\d+)?"  # matches ints or floats of forms '1.105' or '1.105e-5' or '1.105e5' or '1.105e+5'

tree_id_re = r"\d+"  # tree ID is an int
llh_re = rf"{sign}{float_re}"  # likelihood is a signed floating point
deltaL_re = rf"{sign}{float_re}"  # deltaL is a signed floating point
test_result_re = rf"{float_re}{blanks}{sign}"  # test result entry is of form '0.123 +'

# a table entry in the .iqtree file looks for example like this:
# 5 -5708.931281 1.7785e-06  0.0051 -  0.498 +  0.987 +  0.498 +  0.987 +      0.05 +    0.453 +
table_entry = rf"({tree_id_re}){blanks}({llh_re}){blanks}({deltaL_re}){blanks}(?:({test_result_re}){blanks})*"
table_entry_re = regex.compile(table_entry)

START_STRING = "USER TREES"
END_STRING = "TIME STAMP"


def get_relevant_section(input_file):
    """
    Returns the content of input_file between START_STRING and END_STRING.
    Throws ValueError if the section between START_STRING and END_STRING is empty.

    input_file: the iqtree test summary
    """
    with open(input_file) as f:
        content = f.readlines()

    # remove newline character and trailing white space
    content = [l.strip("\n").strip() for l in content]

    # now let's find the relevant lines
    # the relevant lines are only between the start and end string
    start = 0
    end = 0

    for i, line in enumerate(content):
        if START_STRING in line:
            start = i
        if END_STRING in line:
            end = i

    if start == end:
        raise ValueError(
            f"The section between START_STRING {START_STRING} and END_STRING {END_STRING} is empty. Please check the input file {input_file}."
        )

    return content[start:end]


def get_cleaned_table_entries(table_section):
    """
    Returns the cleaned table entries in the given section.
    If the format of the table entries in future IQTREE versions change, make sure to change the defined regex above.
    """
    entries = []
    for line in table_section:
        # match the line against the regex defined above for a table entry
        m = regex.match(table_entry_re, line)
        if m:
            # if a match was found: capture the results in variables
            tree_id, llh, deltaL, result_group = m.groups()
            # to capture all test results individually we have to explicitly unpack it
            test_results = m.captures(4)
            entry = (tree_id, llh, deltaL, test_results)
            entries.append(entry)

    if not entries:
        raise ValueError(
            "No line in the given section matches the regex. Compare the regex and the given section. Maybe the format has changed."
        )

    return entries


def get_indices_of_trees_passing_all_tests(entries):
    def _has_significant_exclusion(test_results):
        return any("-" in t for t in test_results)

    return [
        int(tree_id)
        for (tree_id, _, _, test_results) in entries
        if not _has_significant_exclusion(test_results)
    ]


def save_filtered_trees(ml_trees, indices, save_path):
    with open(ml_trees) as f:
        trees = f.readlines()

    if len(trees) < max(indices):
        raise ValueError(
            "The tree indices and the number of raxml-ng ml_trees do not match. Make sure you input the correct files."
        )

    filtered_trees = [trees[i - 1] for i in indices]
    filtered_trees = "".join(filtered_trees)
    with open(save_path, "w") as w:
        w.write(filtered_trees)


# first get the section containing the result table
section = get_relevant_section(snakemake.input.summary)
entries = get_cleaned_table_entries(section)

# for now: check which trees passed ALL tests (only + in test results)
indices = get_indices_of_trees_passing_all_tests(entries)
save_filtered_trees(snakemake.input.ml_trees, indices, snakemake.output.outfile)