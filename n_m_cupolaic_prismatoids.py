# Script to generate a cupalaic prismatoid based on a {n/m}-gram
import argparse
import logging
from pathlib import Path
import sys

from math import cos, gcd, pi, sin, sqrt
from orbitit import geom_3d as geom
from orbitit.colors import STD_COLORS as cols
from orbitit import geomtypes

DESCRIPTION = """Generate an off file for {n/m} based pseudo-cupolaic prismatoids

The result will be a polyhedron with a {n/m}-gram in the bottom and attached to the edges there will
be triangles. which will be equilateral by default. The polyhedraon will be closed by adding
crossed rectangles, for which the crossing edges are shared with the triangles and the parallel
edges are shared with a neighbouring crossed rectangle.

Note all combinations of n and m have solutions.

In two dimensions there is not difference between for instance a {7/3} and a {7/4}. It is the
convention for e.g. anti-prisms to use {7/4} if retrograde triangles are added. Here only one option
is valid, e.g. for {7/3} you can only get retrograde triangles and for {7/2} you can only get normal
triangles.

The script will however still expect to follow the convention and return an error if the wrong value
of 'm' is used.
"""
LOGGER = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description=DESCRIPTION)
parser.add_argument("n", type=int, help="Number of vertices of the {n/m}-gram.")
parser.add_argument("m", type=int, help="Number of full rounds the {n/m}-gram makes.")
parser.add_argument("filename", help="Path to the off file to be save.")
parser.add_argument(
    "-H", "--allow_holes",
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

if args.n < 5:
    LOGGER.error("n must be bigger than 4")
    sys.exit(1)

odd_n = args.n % 2 != 0
m_even = args.m % 2 == 0
m_first_half = args.m < args.n / 2

no_of_compounds = gcd(args.n, args.m)
no_of_vs_x_gram = args.n // no_of_compounds
if no_of_vs_x_gram == 2:
    LOGGER.error("The provided values for n and m lead to digons: these aren't supported")
    sys.exit(1)

if not m_even:
    LOGGER.error(
        "No solution found for these n and m, try m = %d, using retrograde triangles",
        args.n - args.m,
    )
    sys.exit(1)

if m_first_half:
    v_distance = args.m
else:
    # The smallest value that expresses the vertex jump to make in the n/m-gram without taken into
    # consideration the direction. This is used for the polygram at the bottom, while for even m the
    # vertices at the top follow the args.m
    v_distance = args.n - args.m

# Vertices
# This assumes the side of the n-gon has length 2
RADIUS = 1 / sin(pi / args.n)
TWO_PI = 2 * pi
DIAGONAL = 2 * RADIUS * sin(v_distance * pi / args.n)
# TOP_OFFSET expresses which diagonal/edges the crossed rectangles: if the crossed rectangles create
# a n/x gram, it is the value of x.
TOP_OFFSET = args.m // 2
# For the second half we turn in the opposite direction
if not m_first_half:
    TOP_OFFSET = -TOP_OFFSET

# This gives the following length for the crossed rectangle width (the larger value)
BOWTIE_DIAGONAL = abs(2 * RADIUS * sin(TOP_OFFSET * pi / args.n))

if not args.crossed_squares:
    assert geomtypes.FloatHandler.gt(DIAGONAL, BOWTIE_DIAGONAL), (
        "It is not possible to get equilateral triangles, e.g. the polyhedron might become "
        "completely flat. Try the option --crossed_squares instead"
    )
    HALF_HEIGHT = sqrt(DIAGONAL**2 - BOWTIE_DIAGONAL**2) / 2
else:
    HALF_HEIGHT = BOWTIE_DIAGONAL / 2  # half edge length

# E.g. {7}
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
    geom.vec(
        RADIUS * cos(i * TWO_PI / args.n),
        RADIUS * sin(i * TWO_PI / args.n),
        HALF_HEIGHT,
    ) for i in range(args.n)
]
vs.extend([geom.vec(v[0], v[1], -v[2]) for v in vs])

# Faces and colours
# The {n/2}-gram
n_gram_col = 3
if no_of_compounds > 1:
    # handle e.g. that {9/3} consists of three triangles
    faces = [
        # opposite distance to make sure the normal points outward
        [(m - i * v_distance) % args.n for i in range(no_of_vs_x_gram)]
        for m in range(no_of_compounds)
    ]
    col_i = [n_gram_col for m in range(no_of_compounds)]
else:
    faces = [
        [(-i * v_distance) % args.n for i in range(args.n)]
    ]
    col_i = [n_gram_col]

# The equilateral triangles
faces.extend(
    [
        [i, (i + v_distance) % args.n, (i + TOP_OFFSET) % args.n + args.n] for i in range(args.n)
    ]
)
col_i.extend([2 for _ in range(args.n)])

# The crossed rectangles:
faces.extend(
    [
        [
            i,
            i + args.n,
            (i + TOP_OFFSET) % args.n,
            (i + TOP_OFFSET) % args.n + args.n
        ] for i in range(args.n)
    ]
)
col_i.extend([1 for _ in range(args.n)])

if not args.allow_holes:
    for face_index in range(no_of_compounds):
        n_m_gram = geom.Face([vs[faces[face_index][i]] for i in range(no_of_vs_x_gram)])
        # FIXME: magical constant 8
        with geomtypes.FloatHandler(8):
            extended_n_m_gram = n_m_gram.outline
            offset = len(vs)
            vs.extend(extended_n_m_gram)
            faces[face_index] = [i + offset for i in range(len(extended_n_m_gram))]

shape = geom.SimpleShape(
    vs,
    faces,
    colors=(cols, col_i),
    name=f"{args.n}/2 pseudo-cupolaic prismatoid",
)
filepath = Path(args.filename)
if filepath.is_file():
    yes_or_no = input("File exists. Overwrite? y/N\n")
    if yes_or_no == "" or yes_or_no.lower()[0] != "y":
        sys.exit(1)

with open(args.filename, "w") as fd:
    fd.write(shape.to_off())
