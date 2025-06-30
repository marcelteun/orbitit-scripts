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
bowties, for which the crossing edges are shared with the triangles and the parallel edges are
shared with a neighbouring bowtie.

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

if odd_n:
    if not m_even:
        LOGGER.error(
            "No solution found for these n and m, try m = %d, using retrograde triangles",
            args.n - args.m,
        )
        sys.exit(1)
    if not m_first_half:
        # The code below is for historical reasons: the script only used to support m for m < n/2
        # Then for odd m it would use retrograde triangles by adjusting TOP_OFFSET. This special
        # handling and isn't needed, though in that case the polygrams (at least) would turn inside
        # out.
        args.m = args.n - args.m
        m_even = not m_even
        # not really needed, but better be correct:
        m_first_half = not m_first_half
else:
    if not m_even:
        # TODO: update log: same as abve
        LOGGER.error("Only even m supported for even n (at the moment)")
        sys.exit(1)
    else:
        # FIXME: Fix m for m > n / 2
        if no_of_vs_x_gram == 2:
            LOGGER.error("Values n and m lead to digons: these aren't supported")
            sys.exit(1)

# Vertices
# This assumes the side of the n-gon has length 2
RADIUS = 1 / sin(pi / args.n)
TWO_PI = 2 * pi
DIAGONAL = 2 * RADIUS * sin(args.m * pi / args.n)
# TOP_OFFSET expresses which diagonal/edges the bowties: if the bowties create a n/x gram, it is the
# value of x.
if m_even:
    TOP_OFFSET = args.m // 2
else:
    TOP_OFFSET = -(args.n - args.m) // 2
# This gives the following length for the bowtie width
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
# I think the if statement should be gcd(args.n, args.m) > 1
if no_of_compounds > 1:
    # handle e.g. that {9/3} consists of three triangles
    faces = [
        [(m - i * args.m) % args.n for i in range(no_of_vs_x_gram)]
        for m in range(args.m)
    ]
    col_i = [n_gram_col for m in range(args.m)]
else:
    faces = [
        [(-i * args.m) % args.n for i in range(args.n)]
    ]
    col_i = [n_gram_col]

# The equilateral triangles
faces.extend(
    [
        [i, (i + args.m) % args.n, (i + TOP_OFFSET) % args.n + args.n] for i in range(args.n)
    ]
)
col_i.extend([2 for _ in range(args.n)])

# The bow ties:
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

# TODO: handle 10 / 4 here
if no_of_compounds == 1 and not args.allow_holes:
    n_m_gram = geom.Face([vs[faces[0][i]] for i in range(args.n)])
    with geomtypes.FloatHandler(8):
        extended_n_m_gram = n_m_gram.outline
        offset = len(vs)
        vs.extend(extended_n_m_gram)
        faces[0] = [i + offset for i in range(len(extended_n_m_gram))]

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
