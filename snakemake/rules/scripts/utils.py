import json

RAXML_INDICATOR = "Final LogLikelihood"
IQTREE_INDICATOR = "BEST SCORE FOUND"


def get_parameter_value(filename: str, param_identifier: str) -> float:
    with open(filename) as f:
        data = json.load(f)

    if not param_identifier in data:
        raise ValueError(
            f"The given parameter identifier {param_identifier} is not stored in the given file {filename}."
        )

    return data[param_identifier]


def get_raxml_llh(raxml_file: str) -> float:
    # read final llh as computed in raxml-ng run
    with open(raxml_file) as f:
        content = f.readlines()
    if not any(RAXML_INDICATOR in l for l in content):
        raise ValueError(
            f"The given input file {raxml_file} does not contain a final llh. Make sure the files are correct and the indicator string has not changed."
        )
    llh_raxml = 0

    for line in content:
        if RAXML_INDICATOR in line:
            _, llh = line.split(":")
            llh_raxml = llh.strip()
            break

    return llh_raxml


def get_iqtree_llh(iqtree_file: str) -> float:
    # read final llh as reevaluated by iqtree
    with open(iqtree_file) as f:
        content = f.readlines()
    if not any(IQTREE_INDICATOR in l for l in content):
        raise ValueError(
            f"The given input file {iqtree_file} does not contain a best score. Make sure the files are correct and the indicator string has not changed."
        )

    llh_iqtree = 0

    for line in content:
        if IQTREE_INDICATOR in line:
            _, llh = line.split(":")
            llh_iqtree = llh.strip()
            break

    return llh_iqtree
