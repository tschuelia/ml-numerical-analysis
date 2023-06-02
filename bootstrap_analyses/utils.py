import re


def get_support_values(treefile):
    newick = open(treefile).readline().strip()
    r = re.compile(r"(\d+):")
    m = re.findall(r, newick)
    support_values = [float(el) for el in m]

    if len(support_values) == 1:
        support_values = support_values + support_values

    return support_values


def get_value_from_line(line: str, search_string: str) -> float:
    line = line.strip()
    if search_string in line:
        _, value = line.rsplit(" ", 1)
        return float(value)

    raise ValueError(
        f"The given line '{line}' does not contain the search string '{search_string}'."
    )


def get_rfdist_info(log_file):
    rel_rfdist = None
    num_topos = None

    for line in open(log_file).readlines():
        line = line.strip()

        if "Average relative RF distance in this tree set:" in line:
            rel_rfdist = get_value_from_line(
                line, "Average relative RF distance in this tree set:"
            )
        elif "Number of unique topologies in this tree set:" in line:
            num_topos = get_value_from_line(
                line, "Number of unique topologies in this tree set:"
            )

    if rel_rfdist is None or num_topos is None:
        raise ValueError("Error parsing raxml-ng log.")

    return num_topos, rel_rfdist


def get_raxmlng_elapsed_time(log_file):
    for line in open(log_file).readlines():
        line = line.strip()
        if not "Elapsed time:" in line:
            continue
        # two cases now:
        # either the run was cancelled an rescheduled
        if "restarts" in line:
            # line looks like this: "Elapsed time: 5562.869 seconds (this run) / 91413.668 seconds (total with restarts)"
            _, right = line.split("/")
            value = right.split(" ")[1]
            return float(value)

        # ...or the run ran in one sitting...
        else:
            # line looks like this: "Elapsed time: 63514.086 seconds"
            value = line.split(" ")[2]
            return float(value)

    raise ValueError(
        f"The given input file {log_file} does not contain the elapsed time."
    )


def get_number_of_bs_replicates(log_file):
    number_replicates = 1000

    for line in open(log_file):
        if "Bootstrapping converged" in line:
            r = re.compile(r"(\d+)\s+replicates")
            m = re.findall(r, line)
            number_replicates = int(m[-1])

    return number_replicates