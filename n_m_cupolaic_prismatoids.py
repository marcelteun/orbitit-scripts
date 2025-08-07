"""Script to generate a cupalaic prismatoid based on a {n/m}-gram"""
import argparse
import logging
import os
from pathlib import Path
import sys

from math import cos, gcd, pi, sin, sqrt
from orbitit import geom_3d as geom
from orbitit.colors import STD_COLORS as cols
from orbitit import geomtypes

# TODO: update text below
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

    scaling_slanted = "slanted"
    scaling_vertical = "vertical"

    def __init__(self, n, m, p, allow_holes, scaling=None):
        """Initialise object

        This will create a pseudo-cupoalic prismatoid with bases {n/m} and {n/p}.

        Note that not all combinations of n, m and p will lead to valid PCP.

        n: amount of vertices in the bases
        m: the number m in the primary base {n/m} with m < n/2.
        p: the number p in the secondary base {n/p} with m < n/2. Use p = n to get a pseudo-base.
        scaling: This affects the scaling of the height of the PCP.
            It either evaluates to False or equals to 'slanted' or 'vertical'. For the 'vertical'
            the height will be scaled so that the crossed rectangles so they fit in a square. For
            the 'slanted' it will try to scale in such a way that the slanted faces attached to a
            polygon base through an edge has three edges with the same length. Note that this isn't
            always possible and it might result in an error.
        allow_holes: if set to True then the {n/n}-gram(s) will be replaced by a polygon following
            the outline. If set to False the polygon will follow the n edges, which might not be
            shown well in a 3D player, e.g. holes might appear at parts that have even coverage.
        """
        if n < 3:
            raise ValueError("n must be bigger than 3")

        self.n = n
        self.m = m
        self.p = p
        self.scaling = scaling
        self.allow_holes = allow_holes

        self.pseudo_base = p == n
        m_low = min(m, n - m)
        p_low = min(p, n - p)

        if self.pseudo_base:
            # pseudo-base should have even m
            if m % 2 != 0:
                raise ValueError(f"Expected an even m for a PCP with pseudo-base, got m={m}")
        else:
            p_delta_m_low = p - m_low

            # 2 bases
            if not 0 < m < n:
                raise ValueError(f"Make sure that 0 < m < n, got m = {m} and n = {n}")
            m_low = n - m if m > n / 2 else m
            if not m_low <= p <= n - m_low:
                raise ValueError(f"Make sure that {m_low} <= p <= {n - m_low}, got p={p}")
            if p_delta_m_low % 2 != 0:
                raise ValueError("Invalid difference between m and p")
            if m_low == m and p != m:
                raise ValueError(f"Invalid m and p: did you mean m = {n - m}?")

        # opposite value of m in the modulo-n group
        self.m_opp = n - m

        # any data for this pseudo cupolaic prismatoid
        self._pcp_data = Object()
        self._pcp_data.radius = 1 / sin(pi / n)

        # properties of the primary base polygon {n/m}
        self.bases = [Object()]
        # {n/m} may be a compound, e.g. for {8/2}
        self.bases[0].no_of_compounds = gcd(n, m)
        # The number of vertices per sub-polygon:
        self.bases[0].no_of_vs_x_gram = n // self.bases[0].no_of_compounds
        # always draw in one correction to make it easy to decide the face normal
        self.bases[0].v_distance = m_low

        if self.bases[0].no_of_vs_x_gram == 2:
            if m == p:
                raise ValueError(
                    f"Cannot generate PCP {n}/{m} | {n}/{p} with two digon bases"
                )
            LOGGER.info("The provided values for n (%d) and m (%d) lead to digons", n, m)

        # Calculate the length of the diagonal. This assumes the side of the n-gon has length 2
        self.bases[0].diagonal = 2 * self._pcp_data.radius * sin(self.bases[0].v_distance * pi / n)

        # The crossed rectangles also follow some diagonal in the {n/q}-polygon which results in a
        # {n/to_secondary_offset} polygon, secondary as in secondary base
        if self.pseudo_base:
            self._pcp_data.to_secondary_offset = m // 2
            m_first_half = m <= n / 2
            if not m_first_half:
                self._pcp_data.to_secondary_offset = -self._pcp_data.to_secondary_offset
        else:
            offset = (p - self.m_opp) % n
            # should have been taken care of above, so this shouldn't be a problem
            assert offset % 2 == 0, f"p - (n - m) = {offset} should be even"
            # to_secondary_offset is the δ in pcp_notes.txt
            self._pcp_data.to_secondary_offset = offset // 2

        # The length of the rectangular hull of the crossed rectangle, which is an n/q diagonal
        self._pcp_data.crossed_diagonal = abs(
            2 * self._pcp_data.radius * sin(self._pcp_data.to_secondary_offset * pi / self.n)
        )

        if self.pseudo_base:
            self._pcp_data.half_height = self._get_half_height_pseudo_base()
        else:
            self.bases.append(Object())
            self.bases[1].v_distance = p_low
            self.bases[1].no_of_compounds = gcd(self.n, p)
            self.bases[1].no_of_vs_x_gram = self.n // self.bases[1].no_of_compounds
            self._pcp_data.half_height = self._get_half_height_double_base()

        vertices = self._get_vertices()
        faces, col_i = self.get_n_m_gram(0, False)

        # TODO: use four classes:
        # 1. shared abstract base class
        # 2. pseudo-base class
        # 3. double base class
        # 4. general PCP class
        if self.pseudo_base:
            f, c = self._get_crossed_rectangles_pseudo_base()
            self._pcp_data.crossed_squares_index = len(faces)
            faces.extend(f)
            col_i.extend(c)

            f, c = self._get_slanted_faces_pseudo_base()
            self._pcp_data.slanted_faces_index = len(faces)
            faces.extend(f)
            col_i.extend(c)
        else:
            f, c = self.get_n_m_gram(1, True)
            faces.extend(f)
            col_i.extend(c)

            f, c = self._get_crossed_rectangles_double_base()
            self._pcp_data.crossed_squares_index = len(faces)
            faces.extend(f)
            col_i.extend(c)

            f, c = self._get_slanted_faces_double_base()
            self._pcp_data.slanted_faces_index = len(faces)
            faces.extend(f)
            col_i.extend(c)

        super().__init__(
            vertices,
            faces,
            colors=(cols, col_i),
            name=f"{n}/{m} pseudo-cupolaic prismatoid",
        )

        self.crossed_squares_use_outlines()
        if not self.allow_holes:
            self.use_outlines()

    def _get_half_height_pseudo_base(self):
        """Calculate the height of the polyhedron using a pseudo-base."""
        scaling = self.scaling
        if scaling != self.scaling_vertical:
            cannot_do_slanted = geomtypes.FloatHandler.le(
                self.bases[0].diagonal, self._pcp_data.crossed_diagonal
            )
            if cannot_do_slanted:
                if self.scaling == self.scaling_slanted:
                    raise ValueError(
                        "It is not possible to get equilateral triangles, e.g. the polyhedron "
                        f"might become completely flat. Try the option '{self.scaling_vertical}' "
                        "scaling instead"
                    )
                scaling = self.scaling_vertical

        if scaling == self.scaling_vertical:
            half_height = self._pcp_data.crossed_diagonal / 2  # half edge length
        else:
            half_height = sqrt(self.bases[0].diagonal**2 - self._pcp_data.crossed_diagonal**2) / 2
        return half_height

    def _get_half_height_double_base(self):
        """Calculate the height of the polyhedron with a double base."""
        scaling = self.scaling
        if scaling != self.scaling_vertical:

            #        ______________
            #       /              \
            #      /                \
            #     /                  \
            #    /                    \
            #   +----------------------+
            # There are two options to make the isosceles trapezoid to become trisosceles
            # 1. Adapt the height so the secondary edge gets the same length as the sides
            # 2. Adapt the height so the primary edge gets the same length as the sides
            #
            # The crossed square diagonal should be equal to the diagonal used at the secondary
            # or primary (top_distance).
            # This while the crossed square diagonal = 2 * √(1 + half_height**2)
            roots = []
            for distance in (self.bases[0].v_distance, self.bases[1].v_distance):
                half_crossed_square_diagonal = self._pcp_data.radius * sin(
                    distance * pi / self.n
                )
                root = half_crossed_square_diagonal**2 - 1
                if geomtypes.FloatHandler.gt(root, 0):
                    roots.append(sqrt(root))

            if not roots:
                if self.scaling == self.scaling_slanted:
                    raise ValueError(
                        "It is not possible to get equilateral triangles, e.g. the polyhedron "
                        f"might become completely flat. Try the option '{self.scaling_vertical}' "
                        "scaling instead"
                    )
                scaling = self.scaling_vertical

        if scaling == self.scaling_vertical:
            half_height = self._pcp_data.crossed_diagonal / 2
        else:
            # For now use the minimum
            # TODO: decide how to handle (command line parameter?)
            half_height = min(roots)
            LOGGER.info("Possible heights: %s, using %s", roots, half_height)
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
        """Replace the n-grams by there outlines to prevent holes."""
        face_index = 0
        for base in self.bases:
            for _ in range(base.no_of_compounds):
                self.replace_face_by_outline(face_index, self.exp_tol_eq_float)
                face_index += 1

    def crossed_squares_use_outlines(self):
        """Replace the crossed rectangles by there outlines to flashing."""
        for i in range(self.n):
            self.replace_face_by_outline(
                self._pcp_data.crossed_squares_index + i, self.exp_tol_eq_float
            )
        if self.m == self.p:
            for i in range(self.n):
                face_id = self._pcp_data.slanted_faces_index + i
                self.replace_face_by_outline(face_id, self.exp_tol_eq_float)
                # When the trapezoids become crossed rectangles and the outline is taken then need
                # to be reversed for the normal to point outwards.
                self.reverse_face(face_id)

    def _primary_base_index(self, index):
        "Ensure a vertex index is from the primary base."
        return index % self.n

    def _secondary_base_index(self, index):
        "Ensure a vertex index is from the secondary base."
        return index % self.n + self.n

    def get_n_m_gram(self, offset, opposite_direction):
        """Get a tuple with face and colour indices for the {n/m}-gram

        If the {n/m} is a compound then there will be more than one. This call with set the
        attributes faces and col_i for the n_gram object.

        offset: index in self.bases
        opposite_direction: set to True to turn in the opposite direction. Primary and secondary
            should use opposite direction to get the face normal consistently pointing to the inside
            or the outside.
        """
        n_gram = self.bases[offset]
        inv = -1 if opposite_direction else 1
        v_index = [self._primary_base_index, self._secondary_base_index][offset]
        n_is_even = self.n % 2 == 0
        if n_is_even and n_gram.no_of_vs_x_gram == 2:
            # E.g. this is a {10/5} i.e. 5 digons. For n is even these become edges where the
            # isosceles rectangles meet.
            n_gram.no_of_compounds = 0
            faces = []
            col_i = []
        elif n_gram.no_of_compounds > 1:
            # handle e.g. that {9/3} consists of three triangles
            faces = [
                [
                    v_index(m + inv * i * n_gram.v_distance)
                    for i in range(n_gram.no_of_vs_x_gram)
                ]
                for m in range(n_gram.no_of_compounds)
            ]
            col_i = [self.n_gram_col for m in range(n_gram.no_of_compounds)]
        else:
            faces = [
                [v_index(inv * i * n_gram.v_distance) for i in range(self.n)]
            ]
            col_i = [self.n_gram_col]

        return faces, col_i

    def _get_slanted_faces_pseudo_base(self):
        """Get a tuple with face and colour indices for the slanted faces."""
        # The equilateral triangles
        faces = [
            [
                i,
                (i + self.bases[0].v_distance) % self.n,
                (i + self._pcp_data.to_secondary_offset) % self.n + self.n,
            ]
            for i in range(self.n)
        ]
        return faces, [self.slanted_col for _ in range(self.n)]

    def _get_crossed_rectangles_pseudo_base(self):
        """Get a tuple with face and colour indices for the vertical crossed rectangles."""
        faces = [
            [
                i,
                i + self.n,
                (i + self._pcp_data.to_secondary_offset) % self.n,
                (i + self._pcp_data.to_secondary_offset) % self.n + self.n
            ] for i in range(self.n)
        ]
        return faces, [self.crossed_col for _ in range(self.n)]

    def _get_slanted_faces_double_base(self):
        """Get a tuple with face and colour indices for the slanted faces."""
        # Use isosceles (or trisosceles) trapezoids
        # See file pcp_notes.txt to understand how these values are obtained:
        faces = [
            [
                self._primary_base_index(i),
                self._primary_base_index(i + self.m),
                self._secondary_base_index(i + self.m - self._pcp_data.to_secondary_offset),
                self._secondary_base_index(
                    i + self.m - self._pcp_data.to_secondary_offset + self.p
                ),
            ] for i in range(self.n)
        ]
        return faces, [self.slanted_col for _ in range(self.n)]

    def _get_crossed_rectangles_double_base(self):
        """Get a tuple with face and colour indices for the vertical crossed rectangles."""
        # See file pcp_notes.txt to understand how these values are obtained:
        faces = [
            [
                self._primary_base_index(i + self.m),
                self._secondary_base_index(i + self.m - self._pcp_data.to_secondary_offset),
                self._primary_base_index(i + self.m - self._pcp_data.to_secondary_offset),
                self._secondary_base_index(i + self.m),
            ] for i in range(self.n)
        ]
        return faces, [self.crossed_col for _ in range(self.n)]


if __name__ == "__main__":

    SCALE_1 = PseudoCupolaicPrismatoid.scaling_slanted
    SCALE_2 = PseudoCupolaicPrismatoid.scaling_vertical

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("n", type=int, help="Number of vertices of the {n/m}-polygon.")
    parser.add_argument(
        "-m",
        type=int,
        default=0,
        help="Vertex offset of the primary base polygon {n/m}. Use m > n/2 for "
        "connecting faces that fold inwards. "
        "If not specified then all PCPs for the specified value of 'n' are generated.",
    )
    parser.add_argument(
        "-p",
        type=int,
        default=0,
        help="Vertex offset of the secondary base polygon {n/p}. Use m > n/2 for "
        "connecting faces that fold inwards or p = n for pseudo-base. "
        "If not specified then all PCPs for the specified value of 'n' and 'm' are generated.",
    )
    parser.add_argument(
        "-f", "--file_base_name",
        default="pcp_",
        help="A header to name the file. This will be used as a base for the OFF file. "
        "It will be appended by n_m__n_p.off."
    )
    parser.add_argument(
        "-o", "--out_dir",
        default=".",
        help="path to directory to save the resulting OFF file(s)."
    )
    parser.add_argument(
        "-H", "--allow_holes",
        action="store_true",
        help="If specified the {n/m} polygon will saved using n vertices, which results in holes "
        "in orbitit for the parts that have even coverage due to the stencil buffer.",
    )
    parser.add_argument(
        "-s", "--scaling",
        help=f"This parameter affect the scaling of the height of the PCP. "
        f"If specified this should either be '{SCALE_1}' or '{SCALE_2}'. "
        "For the latter the program fit the crossed rectangles so they fit in a square. "
        "For the former it will try to scale in such a way that the slanted faces attached to a "
        "polygon base through an edge has three edges with the same length. "
        "Note that this isn't always possible and it might result in an error."
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

    if os.path.exists(args.out_dir):
        if not os.path.isdir(args.out_dir):
            raise ValueError(
                f"The output directory {args.out_dir} exists, but isn't a valid directory"
            )
    else:
        os.mkdir(args.out_dir)

    if args.scaling and not args.scaling in [
        PseudoCupolaicPrismatoid.scaling_slanted,
        PseudoCupolaicPrismatoid.scaling_vertical,
    ]:
        raise ValueError(
            f"--scaling should be either 'squared' or 'trisosceles', git {args.scaling}"
        )


    def generate_one_pcp(args, m, p):
        """Generate one PCP for the specified value of p using args and save as OFF file."""
        shape = PseudoCupolaicPrismatoid(args.n, m, p, args.allow_holes, args.scaling)

        if args.x_rotate:
            shape.transform(
                geomtypes.Rot3(angle=geom.DEG2RAD * args.x_rotate, axis=geomtypes.Vec3([1, 0, 0]))
            )

        filepath = Path(args.out_dir) / f"pcp_{args.n}_{m}__{args.n}_{p}.off"
        if not args.overwrite and filepath.is_file():
            yes_or_no = input(f"{filepath} exists. Overwrite? y/N\n")
            if not yes_or_no or yes_or_no.lower()[0] != "y":
                LOGGER.warning("No overwrite requested; bailing out")
                sys.exit(1)

        with open(filepath, "w") as fd:
            minimized_shape = shape.clean_shape(shape.exp_tol_eq_float)
            fd.write(minimized_shape.to_off())
            LOGGER.info("Written %s", filepath)

    if args.m and args.p:
        m_p_values = [(args.m, args.p)]
    else:
        if args.m:
            m_values = [args.m]
        else:
            m_p_values = [
                # two polygon bases:
                (m if m == p else args.n - m, p)
                for m in range(1, args.n // 2 + 1)
                for p in range(m, args.n - m + 1, 2)
                # if not double digons
                if not(args.n / m == 2 and p == m)
            ] + [
                # one polygon base and a pseudo-base
                (m, args.n) for m in range(2, args.n, 2)
            ]

    for m, p in m_p_values:
        if not(args.n / m == 2 and p == m):
            generate_one_pcp(args, m, p)
