"""Script to generate a cupalaic prismatoid based on a {n/m}-gram"""
import argparse
import logging
import os
from pathlib import Path
import sys

from orbitit.base import Orbitit
from orbitit import orbit

DESCRIPTION = """Load an orbit shape in JSON format and save in OFF file format"""

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)-8s - %(message)s',
    datefmt='%m-%d %H:%M',
)
LOGGER = logging.getLogger("convert JSON to OFF")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "json_files",
        nargs="+",
        help="A list of JSON files to convert.",
    )
    parser.add_argument(
        "-o", "--out_dir",
        default=".",
        help="path to directory to save the resulting OFF file(s).",
    )
    ARGS = parser.parse_args()

    if os.path.exists(ARGS.out_dir):
        if not os.path.isdir(ARGS.out_dir):
            raise ValueError(
                f"The output directory {ARGS.out_dir} exists, but isn't a valid directory"
            )
    else:
        os.mkdir(ARGS.out_dir)

    out_dir = Path(ARGS.out_dir)
    for json_file in ARGS.json_files:
        shape = Orbitit.from_json_file(json_file)
        path = Path(json_file).stem
        new_file = out_dir / (path + ".off")
        with open(new_file, "w") as fd:
            fd.write(shape.to_off())
            logging.info("written %s", new_file)
