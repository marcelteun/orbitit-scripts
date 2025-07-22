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
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)-8s - %(message)s',
    datefmt='%m-%d %H:%M',
)
LOGGER = logging.getLogger("pseudo copulaic prismatoid generator")
TWO_PI = 2 * pi


class Object:
    """Just a class to be able to add attributes to an object."""


class PseudoCupolaicPrismatoid(geom.SimpleShape):
    """A pseudo cupolaic prismatoid consists of a {n/m} with vertical crossed rectangles and either
    (equilateral) triangles for even m or isosceles trapezoids for odd m. Note that even crossed
    rectangles and isosceles triangles can be seen as special kinds of trapezoids.
    """

    # When taking the outline of a concave polygon, use this exponent to decide whether two floats
    # are equal
    exp_tol_eq_float = 8

    # colour indices from orbitit.colors.STD_COLORS
    n_gram_col = 3
    crossed_col = 1
    slanted_col = 2

    def __init__(self, n, m, crossed_squares, allow_holes):
        """Initialise object

        n: amount of vertices in the {n/m}-gram
        m: amount of times going around following the edges of the {n/m}-gram
        crossed_squares: whether to use crossed squares. If False then either equilateral triangles
            (for even m) or trisosceles trapezoids (odd m) will be used
        allow_holes: if set to True then the {n/n}-gram(s) will be replaced by a polygon following
            the outline. If set to False the polygon will follow the n edges, which might not be
            shown well in a 3D player, e.g. holes might appear at parts that have even coverage.
        """
        if n < 3:
            raise ValueError("n must be bigger than 3")
        if not 0 < m < n:
            raise ValueError(f"Make sure that 0 < m < {n} got m={m}")

        self.n = n
        self.m = m
        self.crossed_squares = crossed_squares
        self.allow_holes = allow_holes

        self.n_is_even = n % 2 == 0
        self.m_is_even = m % 2 == 0
        self.m_first_half = m < n / 2

        # any data for this pseudo cupolaic prismatoid
        self._pcp_data = Object()
        self._pcp_data.radius = 1 / sin(pi / n)

        # properties of the {n/m}-gram (placed at the bottom)
        self.bottom = Object()
        # {n/m} may be a compound, e.g. for {8/2}
        self.bottom.no_of_compounds = gcd(n, m)
        # The number of vertices per sub-polygon:
        self.bottom.no_of_vs_x_gram = n // self.bottom.no_of_compounds
        self.bottom.vs_offset = 0

        if self.bottom.no_of_vs_x_gram == 2 and self.m_is_even:
            LOGGER.warning("Note, the provided values for n (%d) and m (%d) lead to digons", n, m)

        if self.m_first_half:
            self.bottom.v_distance = m
        else:
            # The smallest value that expresses the vertex jump to make in the n/m-gram without
            # taken into consideration the direction. This is used for the polygram at the bottom,
            # while for even m the vertices at the top follow the m
            self.bottom.v_distance = n - m

        # Vertices
        # This assumes the side of the n-gon has length 2
        self.bottom.diagonal = 2 * self._pcp_data.radius * sin(self.bottom.v_distance * pi / n)

        # Only used for even m:
        # The crossed rectangles also follow some diagonal in the {n}-gram which results in a
        # {n/to_top_offset}-gram
        self._pcp_data.to_top_offset = m // 2
        # For the second half we turn in the opposite direction
        if not self.m_first_half:
            self._pcp_data.to_top_offset = -self._pcp_data.to_top_offset

        # In case of odd m the top becomes a {n/x} gram, which might be a compound as well,
        # initialise to 0 for now.
        self.top = Object()
        self.top.no_of_compounds = 0
        self.top.no_of_vs_x_gram = 0
        self.top.v_distance = 0
        self.top.vs_offset = 0

        if not self.m_is_even:
            if self.m_first_half:
                self.top.v_distance = self.bottom.v_distance - 2
            else:
                self.top.v_distance = self.bottom.v_distance + 2
            self.top.no_of_compounds = gcd(
                self.n,
                min(self.top.v_distance, self.n - self.top.v_distance),
            )
            self.top.no_of_vs_x_gram = self.n // self.top.no_of_compounds
            self.top.vs_offset = n

        self._pcp_data.half_height = self._get_half_height()
        vertices = self._get_vertices()

        faces, col_i = self.get_n_m_gram(self.bottom, True)

        # Only for odd m there is a top n-gram
        if not self.m_is_even:
            f, c = self.get_n_m_gram(self.top, False)
            faces.extend(f)
            col_i.extend(c)

        f, c = self.get_crossed_rectangles()
        faces.extend(f)
        col_i.extend(c)

        f, c = self.get_slanted_faces()
        faces.extend(f)
        col_i.extend(c)

        super().__init__(
            vertices,
            faces,
            colors=(cols, col_i),
            name=f"{n}/{m} pseudo-cupolaic prismatoid",
        )

        if not self.allow_holes:
            self.use_outlines()

    def _get_half_height(self):
        """Calculate the height of the polyhedron."""
        if self.m_is_even:
            # This gives the following length for the crossed rectangle width (the larger value)
            crossed_diagonal = abs(
                2 * self._pcp_data.radius * sin(self._pcp_data.to_top_offset * pi / self.n)
            )
            if not self.crossed_squares:
                assert geomtypes.FloatHandler.gt(self.bottom.diagonal, crossed_diagonal), (
                    "It is not possible to get equilateral triangles, e.g. the polyhedron might "
                    "become completely flat. Try the option with crossed squares instead"
                )
                half_height = sqrt(self.bottom.diagonal**2 - crossed_diagonal**2) / 2
            else:
                half_height = crossed_diagonal / 2  # half edge length
        else:  # m is odd
            if self.crossed_squares:
                half_height = 1
            else:
                #        ______________
                #       /              \
                #      /                \
                #     /                  \
                #    /                    \
                #   +----------------------+
                # There are two options to make the isosceles trapezoid to become trisosceles
                # 1. Adapt the height so the top edge gets the same length as the sides
                # 2. Adapt the height so the bottom edge gets the same length as the sides
                #
                # The crossed square diagonal should be equal to the diagonal used at the top or
                # bottom (top_distance).
                # This while the crossed square diagonal = 2 * √(1 + half_height**2)
                roots = []
                for distance in (self.bottom.v_distance, self.top.v_distance):
                    half_crossed_square_diagonal = self._pcp_data.radius * sin(
                        distance * pi / self.n
                    )
                    root = half_crossed_square_diagonal**2 - 1
                    if geomtypes.FloatHandler.gt(root, 0):
                        roots.append(sqrt(root))
                assert roots, (
                    "It is not possible to get equilateral triangles, e.g. the polyhedron might "
                    "become completely flat. Try with crossed squares instead"
                )
                # For now use the minimum
                # TODO: decide how to handle (command line parameter?)
                LOGGER.info("Possible heights: %s", roots)
                half_height = min(roots)
        return half_height

    def _get_vertices(self):
        """Return a list of vertices

        The vertices are the vertices of an n-gon with length 2 at height self._pcp_data.half_height
        and at -self._pcp_data.half_height
        """
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
                self._pcp_data.radius * cos(i * TWO_PI / self.n),
                self._pcp_data.radius * sin(i * TWO_PI / self.n),
                self._pcp_data.half_height,
            ) for i in range(self.n)
        ]
        vs.extend([geom.vec(v[0], v[1], -v[2]) for v in vs])

        return vs

    def use_outlines(self):
        """Raplace the n-grams by there outlines to prevent holes."""
        for face_index in range(self.bottom.no_of_compounds):
            self.replace_face_by_outline(face_index, self.exp_tol_eq_float)
        for face_index in range(
            self.bottom.no_of_compounds, self.top.no_of_compounds + self.bottom.no_of_compounds
        ):
            self.replace_face_by_outline(face_index, self.exp_tol_eq_float)

    def get_n_m_gram(self, n_gram, opposite_direction):
        """Get a tuple with face and colour indices for the {n/m}-gram

        If the {n/m} is a compound then there will be more than one. This call with set the
        attributes faces and col_i for the n_gram object.

        n_gram: either self.bottom or self.top attribute
        opposite_direction: set to True to turn in the opposite direction. Top and bottom should use
            opposite direction to get the face normal consistently pointing to the inside or the
            outside.
        """
        inv = -1 if opposite_direction else 1
        if self.n_is_even and n_gram.no_of_vs_x_gram == 2:
            # E.g. this is a {10/5} i.e. 5 digons. For n is even these become edges where the
            # isosceles rectangles meet.
            n_gram.no_of_compounds = 0
            faces = []
            col_i = []
        elif n_gram.no_of_compounds > 1:
            # handle e.g. that {9/3} consists of three triangles
            faces = [
                [
                    (m + inv * i * n_gram.v_distance) % self.n + n_gram.vs_offset
                    for i in range(n_gram.no_of_vs_x_gram)
                ]
                for m in range(n_gram.no_of_compounds)
            ]
            col_i = [self.n_gram_col for m in range(n_gram.no_of_compounds)]
        else:
            faces = [
                [(inv * i * n_gram.v_distance) % self.n + n_gram.vs_offset for i in range(self.n)]
            ]
            col_i = [self.n_gram_col]

        return faces, col_i

    def get_crossed_rectangles(self):
        """Get a tuple with face and colour indices for the vertical crossed rectangles."""
        if self.m_is_even:
            faces = [
                [
                    i,
                    i + self.n,
                    (i + self._pcp_data.to_top_offset) % self.n,
                    (i + self._pcp_data.to_top_offset) % self.n + self.n
                ] for i in range(self.n)
            ]
        else:
            # For odd m the crossed rectangles are along the prism sides
            faces = [
                [
                    i,
                    i + self.n,
                    (i + 1) % self.n,
                    (i + 1) % self.n + self.n
                ] for i in range(self.n)
            ]
        return faces, [self.crossed_col for _ in range(self.n)]

    def get_slanted_faces(self):
        """Get a tuple with face and colour indices for the slanted faces."""
        if self.m_is_even:
            # The equilateral triangles
            faces = [
                [
                    i,
                    (i + self.bottom.v_distance) % self.n,
                    (i + self._pcp_data.to_top_offset) % self.n + self.n,
                ]
                for i in range(self.n)
            ]
        else:
            # Use isosceles (or trisosceles) trapezoids
            faces = [
                [
                    i,
                    (i + 1) % self.n + self.n,
                    (i + self.m - 1) % self.n + self.n,
                    (i + self.m) % self.n,
                ] for i in range(self.n)
            ]
        return faces, [self.slanted_col for _ in range(self.n)]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("n", type=int, help="Number of vertices of the {n/m}-gram.")
    parser.add_argument("m", type=int, help="Number of full rounds the {n/m}-gram makes.")
    parser.add_argument("filename", help="Path to the off file to be save.")
    parser.add_argument(
        "-H", "--allow_holes",
        action="store_true",
        help="If specified the {n/m} polygon will saved using n vertices, which results in holes "
        "in orbitit for the parts that have even coverage due to the stencil buffer.",
    )
    parser.add_argument(
        "-s", "--crossed_squares",
        action="store_true",
        help="If specified the {n/m} polygon will be kept as is, which results in holes in orbitit "
        "for a 3D player using a stencil buffer."
    )
    parser.add_argument(
        "-w", "--overwrite",
        action="store_true",
        help="If specified an existing file will be overwritten without asking. Otherwise the "
        "the script will ask interactively whether to overwrite an existing file."
    )
    parser.add_argument(
        "-x", "--x-rotate",
        metavar="DEG",
        type=float,
        help="Rotate the model a certain amount of degrees around the x-axis."
    )
    args = parser.parse_args()
    shape = PseudoCupolaicPrismatoid(args.n, args.m, args.crossed_squares, args.allow_holes)

    if args.x_rotate:
        shape.transform(
            geomtypes.Rot3(angle=geom.DEG2RAD * args.x_rotate, axis=geomtypes.Vec3([1, 0, 0]))
        )

    filepath = Path(args.filename)
    if not args.overwrite and filepath.is_file():
        yes_or_no = input(f"{filepath} exists. Overwrite? y/N\n")
        if not yes_or_no or yes_or_no.lower()[0] != "y":
            LOGGER.warning("No overwrite requested; bailing out")
            sys.exit(1)

    with open(args.filename, "w") as fd:
        fd.write(shape.to_off())
        LOGGER.info("Written %s", filepath)
