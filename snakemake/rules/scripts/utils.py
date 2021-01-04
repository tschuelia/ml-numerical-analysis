import json
import regex


def get_parameter_value(filename: str, param_identifier: str) -> float:
    with open(filename) as f:
        data = json.load(f)

    if not param_identifier in data:
        raise ValueError(
            f"The given parameter identifier {param_identifier} is not stored in the given file {filename}."
        )

    return data[param_identifier]


def _get_value_from_file(input_file, search_string):
    with open(input_file) as f:
        lines = f.readlines()

    for l in lines:
        l = l.strip()
        if search_string in l:
            # a line looks for example like this
            # Number of unique topologies in this tree set: 100
            _, value = l.rsplit(" ", 1)
            return float(value)

    raise ValueError(
        f'The given input file {input_file} does not contain the search string "{search_string}".'
    )


def get_raxml_llh(raxml_file: str) -> float:
    STR = "Final LogLikelihood:"
    return _get_value_from_file(raxml_file, STR)


def get_iqtree_llh(iqtree_file: str) -> float:
    STR = "BEST SCORE FOUND :"
    return _get_value_from_file(iqtree_file, STR)


def get_raxml_abs_rf_distance(log_file: str) -> float:
    STR = "Average absolute RF distance in this tree set:"
    return _get_value_from_file(log_file, STR)


def get_raxml_rel_rf_distance(log_file: str) -> float:
    STR = "Average relative RF distance in this tree set:"
    return _get_value_from_file(log_file, STR)


def get_raxml_num_unique_topos(log_file: str) -> int:
    STR = "Number of unique topologies in this tree set:"
    return _get_value_from_file(log_file, STR)