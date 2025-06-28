# Script to generate a cupalaic prismatoid based on a {n/m}-gram
import argparse
from pathlib import Path
import sys

from math import cos, pi, sin, sqrt
from orbitit import geom_3d as geom
from orbitit.colors import STD_COLORS as cols
from orbitit import geomtypes

parser = argparse.ArgumentParser(
    description="Generate an off file for {n/m} based pseudo-cupolaic prismatoids",
)
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
    print("n must be bigger than 4")
    sys.exit(1)

if 2 * args.m > args.n:
    print("Value of 'm' too big. Make sure that 2 * m < n")
    sys.exit(1)

# Vertices
# This assumes the side of the n-gon has length 2
RADIUS = 1 / sin(pi / args.n)
TWO_PI = 2 * pi
DIAGONAL = 2 * RADIUS * sin(args.m * pi / args.n)
odd_n = args.n % 2 != 0
m_even = args.m % 2 == 0
# TOP_OFFSET expresses which diagonal/edges the bowties: if the bowties create a n/x gram, it is the
# value of x.
if m_even:
    TOP_OFFSET = args.m // 2
else:
    TOP_OFFSET = -(args.n - args.m) // 2
# This gives the following length for the bowtie width
BOWTIE_DIAGONAL = abs(2 * RADIUS * sin(TOP_OFFSET * pi / args.n))

if not args.crossed_squares:
    HALF_HEIGHT = sqrt(DIAGONAL**2 - BOWTIE_DIAGONAL**2) / 2
    assert geomtypes.FloatHandler.gt(DIAGONAL, BOWTIE_DIAGONAL), (
        "It is not possible to get equilateral triangles, since the polyhedron will become "
        "completely flat. Try the option --crossed_squares instead"
    )
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
# FIXME: don't split in odd / even but in args.n % args.m
if odd_n:
    if args.n % args.m == 0:
        # handle e.g. that {9/3} consists of three triangles
        faces = [
            [(m - i * args.m) % args.n for i in range(args.n // args.m)]
            for m in range(args.m)
        ]
        col_i = [n_gram_col for m in range(args.m)]
    else:
        faces = [
            [(-i * args.m) % args.n for i in range(args.n)]
        ]
        col_i = [n_gram_col]
else:
    assert args.m == 2, "Only m=2 supported for even n (at the moment)"
    # handle e.g. that {8/2} consists of two squares
    faces = [
        [(m - i * args.m) % args.n for i in range(args.n // args.m)]
        for m in range(args.m)
    ]
    col_i = [n_gram_col for m in range(args.m)]

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

if args.n % args.m != 0 and not args.allow_holes:
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
