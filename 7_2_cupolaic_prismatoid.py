# Script to generate a cupalaic prismatoid based on a heptagram
from math import cos, pi, sin, sqrt
from orbitit import geom_3d as geom
from orbitit.colors import STD_COLORS as cols
from orbitit import geomtypes

# TODO: use optional parameters for these.

# Set to True and the concave points of the heptagram will be added to the vertices and polygon will
# be the outline with 14 vertices. This prevents the stencil buffer to remove the central heptagon
# when drawing the heptagram.
REPLACE_HEPTAGRAM = True

# If set to True then the triangles will be equilateral. Otherwise the crossed square, the bow-ties,
# will fit into a square.
USE_EQUILATERAL_TRIANGLES = True

# Vertices
# This assumes the side of the heptagon has length 2
radius = 1 / sin(pi / 7)
two_pi = 2 * pi
rho = 2 * radius * sin(two_pi / 7)

if USE_EQUILATERAL_TRIANGLES:
    z = sqrt(rho**2 - 4) / 2  # 4 = squared edge length
else:
    z = 1  # half edge length
print(f"rho = {rho}")
#  Y ^
#    |           v0
#    |    v6            v1
#    |
#    |
#    |  v5                v2
#    |
#    |
#    |       v4     v3
vs = [geom.vec(radius * cos(i * two_pi / 7), radius * sin(i * two_pi / 7), z) for i in range(7)]
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

if REPLACE_HEPTAGRAM:
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
)
print(shape.to_off())
