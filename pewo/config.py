"""
This is a config helper module. It contains helper functions
to retrieve values from the config files.
"""

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


from enum import Enum
from typing import Any, Dict
from pathlib import Path
from pewo.software import PlacementSoftware, AlignmentSoftware, CustomScripts


class ReadGen(Enum):
    ART = "ART",
    PARTITION = "partition"

class Mode(Enum):
    ACCURACY = 0,
    LIKELIHOOD = 1,
    RESOURCES = 2


def get_work_dir(config: Dict) -> str:
    """
    Returns working directory path, turned into an absolute path.
    This is the root directory of PEWO output.
    """
    return str(Path(config["workdir"]).resolve())


def get_read_generator(config: Dict) -> str:
    """
    Returns working directory path. This is the root directory of PEWO output.
    """
    if "read_generator" in config:
        if not any(config["read_generator"] for val in ReadGen):
            raise(RuntimeError(f"Uknown read generator: "
                               f"{config['read_generator']}"))
        return config["read_generator"]

    return ReadGen.PARTITION.value


def is_supported(software: Any) -> bool:
    """
    Checks if software is supported. Takes anything as input, returns True
    if the input parameter is PlacementSoftware, AlignmentSoftware or
    a custom script name.
    """
    return type(software) == PlacementSoftware or \
           type(software) == AlignmentSoftware or \
           type(software) == CustomScripts


def software_tested(config: Dict, software: PlacementSoftware) -> bool:
    """
    Checks if given software is being tested.
    """
    return software.value in config["test_soft"]


def get_mode(config: Dict) -> Mode:
    if "mode" in config:
        mode_dict = dict((m.name.lower(), m) for m in Mode)
        mode_name = config["mode"].lower()
        assert mode_name in mode_dict, f"Wrong mode value: {mode_name}"
        return mode_dict[mode_name]

    raise RuntimeError(f"PEWO mode not specified in the config file. See config.yaml for details")


def get_query_user(config: Dict) -> str:
    """
    Returns if PEWO should generate reads from the input tree.
    """
    if "query_user" not in config:
        return ''

    return str(Path(config["query_user"]).resolve())


def query_user(config: Dict) -> bool:
    """
    Returns if PEWO should generate reads from the input tree.
    """
    return "query_user" in config and config["query_user"]
