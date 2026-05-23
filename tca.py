"""Generate Tri-Composite Antiprism."""
import argparse
import logging
import os
import sys

from orbitit import geom_3d, geomtypes, n_gons

# TODO: update text below
DESCRIPTION = """Generate an off file for a tri-composite antiprism.
"""
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s - %(message)s",
    datefmt="%m-%d %H:%M",
)
LOGGER = logging.getLogger("TCA")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("n", type=int, help="Number of vertices of the {n/m}-polygon.")
    parser.add_argument(
        "-m",
        type=int,
        default=1,
        help="Vertex offset of the primary base polygon {n/m}. Use m > n/2 for "
        "connecting faces that fold inwards. ",
    )
    parser.add_argument(
        "-b", "--file_base_name",
        default="tca_",
        help="A header to name of the file. This will be used as a base for the OFF file. "
        "It will be appended by n_m.off.",
    )
    parser.add_argument(
        "-t", "--file_tail_name",
        default="",
        help="A string to append to name of the file. This will be appended to 'n_m'.",
    )
    parser.add_argument(
        "-o", "--out_dir",
        default=".",
        help="path to directory to save the resulting OFF file(s).",
    )
    parser.add_argument(
        "-H", "--allow_holes",
        action="store_true",
        help="If specified the {n/m} polygon will saved using n vertices, which results in holes "
        "in orbitit for the parts that have even coverage due to the stencil buffer. "
        "If not specified the n/q polygons and the crossed rectangles will be replaced by "
        "their outline, which results in OFF files where edges are broken and hence the will "
        "not have an even amount of faces joining in each edge, which might result in warnings "
        "or errors for some 3D programs.",
    )
    parser.add_argument(
        "-w", "--overwrite",
        action="store_true",
        help="If specified an existing file will be overwritten without asking. Otherwise the "
        "the script will ask interactively whether to overwrite an existing file.",
    )
    parser.add_argument(
        "-x", "--x-rotate",
        metavar="DEG",
        type=float,
        help="Rotate the model a certain amount of degrees around the x-axis.",
    )
    ARGS = parser.parse_args()

    if os.path.exists(ARGS.out_dir):
        if not os.path.isdir(ARGS.out_dir):
            raise ValueError(
                f"The output directory {ARGS.out_dir} exists, but isn't a valid directory"
            )
    else:
        os.mkdir(ARGS.out_dir)

    shape = n_gons.ThreeAntiPrisms(ARGS.n, ARGS.m, use_outline=not ARGS.allow_holes)

    if ARGS.x_rotate:
        shape.transform(
            geomtypes.Rot3(
                angle=geom_3d.DEG2RAD * ARGS.x_rotate,
                axis=geomtypes.Vec3([1, 0, 0]),
            )
        )

    file = f"{ARGS.out_dir}/{ARGS.file_base_name}{ARGS.n}_{ARGS.m}{ARGS.file_tail_name}.off"
    if os.path.exists(file) and not ARGS.over_write:
        LOGGER.warning("%s exists, use '--overwrite' to overwrite.", file)
        sys.exit(1)
    with open(file, "w") as fd:
        fd.write(shape.to_off())
        LOGGER.info("writtten %s", file)
