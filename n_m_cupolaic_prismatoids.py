"""Script to generate a cupalaic prismatoid based on a {n/m}-gram"""
import abc
import argparse
import logging
import os
from pathlib import Path
import sys

from math import cos, gcd, pi, sin, sqrt
from orbitit.colors import STD_COLORS as cols
from orbitit import geom_3d, geomtypes

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

class PcpAbcMeta(abc.ABC, type(geom_3d.SimpleShape)):
    """Meta class for ABC and SimpleShape."""

class PcpAbc(abc.ABC, metaclass=PcpAbcMeta):
    """Abstract class for ABC and SimpleShape."""

    @abc.abstractmethod
    def _calc_to_secondary_offset(self):
        """Return the vertex index difference from primary base to secondary base.

        The crossed rectangles also follow some diagonal in the {n/q}-polygon which results in a
        {n/to_secondary_offset} polygon, secondary as in secondary base
        """

    @abc.abstractmethod
    def _calc_half_height(self):
        """Calculate the height of the polyhedron using the scaling."""

    @abc.abstractmethod
    def _get_crossed_rectangles(self):
        """Get a tuple with face and colour indices for the vertical crossed rectangles."""

    @abc.abstractmethod
    def _get_slanted_faces(self):
        """Get a tuple with face and colour indices for the slanted faces."""

class PcpBase(PcpAbc, geom_3d.SimpleShape):
    """Abstract base class for pseudo cupolaic prismatoid

    A pseudo cupolaic prismatoid consists of a {n/m} polygon and a {n/p} polygon or a pseudo-base
    consisting only using the n vertices of {n}. Generally it has n vertical crossed rectangles and
    either n (equilateral) triangles or n isosceles trapezoids.
    Note that even crossed rectangles and isosceles triangles can be seen as special kinds of
    trapezoids.
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

    def __init__(self, n, m, p, use_outlines, scaling=None):
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
        use_outlines: if set to True then the {n/m}-gram(s) will be replaced by a polygon following
            the outline. If set to False the polygon will follow the n edges, which might not be
            shown well in a 3D player, e.g. holes might appear at parts that have even coverage.
        """
        PcpAbc.__init__(self)

        if n < 3:
            raise ValueError("n must be bigger than 3")

        self.n = n
        self.m = m
        self.p = p
        self.scaling = scaling
        self.use_outlines = use_outlines

        m_low = min(m, n - m)

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
        self.bases[0].flip_face = False
        self.bases[0].ensure_index = self._primary_base_index

        if self.bases[0].no_of_vs_x_gram == 2:
            LOGGER.info("The provided values for n (%d) and m (%d) lead to digons", n, m)

        # Calculate the length of the diagonal. This assumes the side of the n-gon has length 2
        self.bases[0].diagonal = 2 * self._pcp_data.radius * sin(self.bases[0].v_distance * pi / n)

        self._pcp_data.to_secondary_offset = self._calc_to_secondary_offset()

        # The length of the rectangular hull of the crossed rectangle, which is an n/q diagonal
        self._pcp_data.crossed_diagonal = abs(
            2 * self._pcp_data.radius * sin(self._pcp_data.to_secondary_offset * pi / self.n)
        )

        self._add_2nd_base()
        self._pcp_data.half_height = self._calc_half_height()

        vertices = self._get_vertices()

        faces = []
        col_i = []
        for base in self.bases:
            f, c = self.get_n_m_gram(base)
            faces.extend(f)
            col_i.extend(c)

        f, c = self._get_crossed_rectangles()
        self._pcp_data.crossed_squares_index = len(faces)
        faces.extend(f)
        col_i.extend(c)

        f, c = self._get_slanted_faces()
        self._pcp_data.slanted_faces_index = len(faces)
        faces.extend(f)
        col_i.extend(c)

        super().__init__(
            vertices,
            faces,
            colors=(cols, col_i),
            name=f"{n}/{m} | {n} | {p} pseudo-cupolaic prismatoid",
        )

        if self.use_outlines:
            LOGGER.info("Using outlines, OFF files will not load in Stella!")
            self.crossed_squares_use_outlines()
            self.bases_use_outlines()
        else:
            LOGGER.info("Not using outlines, OFF files might show holes in Orbitit!")

    def _add_2nd_base(self):
        """Possibly add secondary base."""

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
            geom_3d.vec(
                self._pcp_data.radius * cos(i * TWO_PI / self.n),
                self._pcp_data.radius * sin(i * TWO_PI / self.n),
                self._pcp_data.half_height,
            ) for i in range(self.n)
        ]
        vs.extend([geom_3d.vec(v[0], v[1], -v[2]) for v in vs])

        return vs

    def _primary_base_index(self, index):
        "Ensure a vertex index is from the primary base."
        return index % self.n

    def _secondary_base_index(self, index):
        "Ensure a vertex index is from the secondary base."
        return index % self.n + self.n

    def get_n_m_gram(self, n_gram):
        """Get a tuple with face and colour indices for the {n/m}-gram

        If the {n/m} is a compound then there will be more than one. This call with set the
        attributes faces and col_i for the n_gram object.

        offset: index in self.bases
        """
        inv = -1 if n_gram.flip_face else 1
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
                    n_gram.ensure_index(m + inv * i * n_gram.v_distance)
                    for i in range(n_gram.no_of_vs_x_gram)
                ]
                for m in range(n_gram.no_of_compounds)
            ]
            col_i = [self.n_gram_col for m in range(n_gram.no_of_compounds)]
        else:
            faces = [
                [n_gram.ensure_index(inv * i * n_gram.v_distance) for i in range(self.n)]
            ]
            col_i = [self.n_gram_col]

        return faces, col_i

    def bases_use_outlines(self):
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

class PcpPseudoBase(PcpBase):
    """This class of pseudo cupolaic prismatoids (PCPs) consists of one polygon base {n/m} and a
    pseudo-base, generally indicated by {n/n}. Besides that it has vertical crossed rectangles and
    n (equilateral) triangles.
    """

    def __init__(self, n, m, use_outlines, scaling=None):
        """Initialise object

        n: amount of vertices in the bases
        m: the number m in the primary base {n/m} with m < n/2.
        scaling: This affects the scaling of the height of the PCP.
            It either evaluates to False or equals to 'slanted' or 'vertical'. For the 'vertical'
            the height will be scaled so that the crossed rectangles so they fit in a square. For
            the 'slanted' it will try to scale in such a way that the slanted faces attached to a
            polygon base through an edge has three edges with the same length. Note that this isn't
            always possible and it might result in an error.
        use_outlines: if set to True then the {n/n}-gram(s) will be replaced by a polygon following
            the outline. If set to False the polygon will follow the n edges, which might not be
            shown well in a 3D player, e.g. holes might appear at parts that have even coverage.
        """
        # pseudo-base should have even m
        if m % 2 != 0:
            raise ValueError(f"Expected an even m for a PCP with pseudo-base, got m={m}")

        PcpBase.__init__(self, n, m, n, use_outlines, scaling)

    def _calc_to_secondary_offset(self):
        to_secondary_offset = self.m // 2
        m_first_half = self.m <= self.n / 2
        if not m_first_half:
            to_secondary_offset = -to_secondary_offset

        return to_secondary_offset

    def _calc_half_height(self):
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

    def _get_crossed_rectangles(self):
        faces = [
            [
                i,
                i + self.n,
                (i + self._pcp_data.to_secondary_offset) % self.n,
                (i + self._pcp_data.to_secondary_offset) % self.n + self.n
            ] for i in range(self.n)
        ]
        return faces, [self.crossed_col for _ in range(self.n)]

    def _get_slanted_faces(self):
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

class PcpPseudoTwoBases(PcpBase):
    """This class of pseudo cupolaic prismatoids (PCPs) consists of one polygon base {n/m} and a
    another base {n/p} Besides that it has vertical crossed rectangles and
    n trapezoids.
    """

    def __init__(self, n, m, p, use_outlines, scaling=None):
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
        use_outlines: if set to True then the {n/n}-gram(s) will be replaced by a polygon following
            the outline. If set to False the polygon will follow the n edges, which might not be
            shown well in a 3D player, e.g. holes might appear at parts that have even coverage.
        """
        m_low = min(m, n - m)
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

        if m == p and n == 2 * m:
            raise ValueError(f"Cannot generate PCP {n}/{m} | {n}/{p} with two digon bases")

        # make vertical scaling the default
        if not scaling:
            scaling = self.scaling_vertical

        PcpBase.__init__(self, n, m, p, use_outlines, scaling)

    def _calc_to_secondary_offset(self):
        offset = (self.p - self.m_opp) % self.n
        # should have been taken care of above, so this shouldn't be a problem
        assert offset % 2 == 0, f"p - (n - m) = {offset} should be even"
        # to_secondary_offset is the δ in pcp_notes.txt
        return offset // 2

    def _add_2nd_base(self):
        self.bases.append(Object())
        p_low = min(self.p, self.n - self.p)
        self.bases[1].v_distance = p_low
        self.bases[1].no_of_compounds = gcd(self.n, self.p)
        self.bases[1].no_of_vs_x_gram = self.n // self.bases[1].no_of_compounds
        self.bases[1].flip_face = True
        self.bases[1].ensure_index = self._secondary_base_index

    def _calc_half_height(self):
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

    def _get_crossed_rectangles(self):
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

    def _get_slanted_faces(self):
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

def PseudoCupolaicPrismatoid(n, m, p, use_outlines, scaling=None):
    """Return a the correct PCP class.

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
    use_outlines: if set to True then the {n/n}-gram(s) will be replaced by a polygon following
        the outline. If set to False the polygon will follow the n edges, which might not be
        shown well in a 3D player, e.g. holes might appear at parts that have even coverage.
    """
    pseudo_base = p == n

    if pseudo_base:
        pcp = PcpPseudoBase(n, m, use_outlines, scaling)
    else:
        pcp = PcpPseudoTwoBases(n, m, p, use_outlines, scaling)
    return pcp


if __name__ == "__main__":

    SCALE_1 = PcpBase.scaling_slanted
    SCALE_2 = PcpBase.scaling_vertical

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
        "in orbitit for the parts that have even coverage due to the stencil buffer. "
        "If not specified the n/q polygons and the crossed rectangles will be replaced by their "
        "outline, which results in OFF files where edges are broken and hence the will not have "
        "an even amount of faces joining in each edge, which might result in warnings or errors "
        "for some 3D programs.",
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
    ARGS = parser.parse_args()

    if os.path.exists(ARGS.out_dir):
        if not os.path.isdir(ARGS.out_dir):
            raise ValueError(
                f"The output directory {ARGS.out_dir} exists, but isn't a valid directory"
            )
    else:
        os.mkdir(ARGS.out_dir)

    if ARGS.scaling and not ARGS.scaling in [
        PcpBase.scaling_slanted,
        PcpBase.scaling_vertical,
    ]:
        raise ValueError(
            f"--scaling should be either '{SCALE_1}' or '{SCALE_2}', got {ARGS.scaling}"
        )


    def generate_one_pcp(m, p):
        """Generate one PCP for the specified value of p using ARGS and save as OFF file."""
        shape = PseudoCupolaicPrismatoid(ARGS.n, m, p, not ARGS.allow_holes, ARGS.scaling)

        if ARGS.x_rotate:
            shape.transform(
                geomtypes.Rot3(
                    angle=geom_3d.DEG2RAD * ARGS.x_rotate,
                    axis=geomtypes.Vec3([1, 0, 0]),
                )
            )

        filepath = Path(ARGS.out_dir) / f"pcp_{ARGS.n}_{m}__{ARGS.n}_{p}.off"
        if not ARGS.overwrite and filepath.is_file():
            yes_or_no = input(f"{filepath} exists. Overwrite? y/N\n")
            if not yes_or_no or yes_or_no.lower()[0] != "y":
                LOGGER.warning("No overwrite requested; bailing out")
                sys.exit(1)

        with open(filepath, "w") as fd:
            minimized_shape = shape.clean_shape(shape.exp_tol_eq_float)
            fd.write(minimized_shape.to_off())
            LOGGER.info("Written %s", filepath)

    if ARGS.m and ARGS.p:
        m_p_values = [(ARGS.m, ARGS.p)]
    else:
        if ARGS.m:
            m_values = [ARGS.m]
        else:
            m_p_values = [
                # two polygon bases:
                (m if m == p else ARGS.n - m, p)
                for m in range(1, ARGS.n // 2 + 1)
                for p in range(m, ARGS.n - m + 1, 2)
                # if not double digons
                if not(ARGS.n / m == 2 and p == m)
            ] + [
                # one polygon base and a pseudo-base
                (m, ARGS.n) for m in range(2, ARGS.n, 2)
            ]

    for primary_m, secondary_m in m_p_values:
        if not(ARGS.n / primary_m == 2 and secondary_m == primary_m):
            generate_one_pcp(primary_m, secondary_m)
