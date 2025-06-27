# Script to generate a cupalaic prismatoid based on a heptagram
import argparse
from pathlib import Path
import sys

from math import cos, pi, sin, sqrt
from orbitit import geom_3d as geom
from orbitit.colors import STD_COLORS as cols
from orbitit import geomtypes

parser = argparse.ArgumentParser(
    description="Generate an off file for {7/2} based pseudo-cupolaic prismatoids",
)
parser.add_argument("filename", help="Path to the off file to be save.")
parser.add_argument(
    "-k", "--allow_holes",
    action="store_true",
    help="If specified the {n/m} polygon will saved using n vertices, which results in holes in "
    "orbitit for the parts that have even coverage due to the stencil buffer.",
)
parser.add_argument(
    "-s", "--crossed_squares",
    action="store_true",
    help="If specified the {n/m} polygon will be kept as is, which results in holes in orbitit for "
)
args = parser.parse_args()

# TODO: generalise this script for any n/m

# Vertices
# This assumes the side of the heptagon has length 2
RADIUS = 1 / sin(pi / 7)
TWO_PI = 2 * pi
DIAGONAL = 2 * RADIUS * sin(TWO_PI / 7)

if not args.crossed_squares:
    HEIGHT = sqrt(DIAGONAL**2 - 4) / 2  # 4 = squared edge length
else:
    HEIGHT = 1  # half edge length
#  Y ^
#    |           v0
#    |    v6            v1
#    |
#    |
#    |  v5                v2
#    |
#    |
#    |       v4     v3
vs = [
    geom.vec(RADIUS * cos(i * TWO_PI / 7), RADIUS * sin(i * TWO_PI / 7), HEIGHT) for i in range(7)
]
vs.extend([geom.vec(v[0], v[1], -v[2]) for v in vs])
# Sides and colours
faces = [
    [0, 5, 3, 1, 6, 4, 2],
    [0, 2, 8],
    [1, 3, 9],
    [2, 4, 10],
    [3, 5, 11],
    [4, 6, 12],
    [5, 0, 13],
    [6, 1, 7],
    [0, 7, 1, 8],
    [1, 8, 2, 9],
    [2, 9, 3, 10],
    [3, 10, 4, 11],
    [4, 11, 5, 12],
    [5, 12, 6, 13],
    [6, 13, 0, 7],
]

if not args.keep_n_gram:
    heptagram = geom.Face([vs[faces[0][i]] for i in range(7)])
    with geomtypes.FloatHandler(8):
        extended_heptagram = heptagram.outline
        offset = len(vs)
        vs.extend(extended_heptagram)
        faces[0] = [i + offset for i in range(len(extended_heptagram))]

shape = geom.SimpleShape(
    vs,
    faces,
    colors=(cols, [0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]),
    name="7/2 pseudo-cupolaic prismatoid",
)
filepath = Path(args.filename)
if filepath.is_file():
    yes_or_no = input("File exists. Overwrite? y/N\n")
    if yes_or_no == "" or yes_or_no.lower()[0] != "y":
        sys.exit(1)

with open(args.filename, "w") as fd:
    fd.write(shape.to_off())
