"""Script to generate a cupalaic prismatoid based on a {n/m}-gram"""
import argparse
import logging
from math import acos, cos, gcd, pi, sin, sqrt
import os
from pathlib import Path
import sys

from scipy.optimize import minimize
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



def get_polygon(n, m, edge_length=2):
    """Return a the base shape.

    The vertices are the vertices of an n-gon z = 0. The edges of the n/m gram will
    scaled to edge_length.

    edge_length: the required edge length

    return: a list with vertices (geomtypes.Vec3)
    """
    # E.g. {7/2}
    #  Y ^
    #    |           v0
    #    |    v4            v3
    #    |
    #    +
    #    |  v1                v6
    #    |
    #    |
    #    |       v5     v2
    #     -----------+-------------> X
    no_of_compounds = gcd(n, m)
    no_of_vs_x_gram = n // no_of_compounds
    if no_of_compounds > 1:
        LOGGER.info("{%i/%i} has a common divider", n, m)
    vs = [
        [
            geom_3d.vec(
                cos((m * i + j) * TWO_PI / n),
                sin((m * i + j) * TWO_PI / n),
                0,
            )
            for i in range(no_of_vs_x_gram)
        ]
        for j in range(no_of_compounds)
    ]

    diagonal = (vs[0][0] - vs[0][1]).norm()
    scale = edge_length / diagonal
    return [[scale * v for v in v_list] for v_list in vs]


class Object:
    """Just a class to be able to add attributes to an object."""

class Cupola(geom_3d.SimpleShape):
    """class for creating cupola object.
    """

    # When taking the outline of a concave polygon, use this exponent to decide whether two floats
    # are equal
    exp_tol_eq_float = 8

    # colour indices from orbitit.colors.STD_COLORS
    base_col = 3
    triangle_col = 1
    side_polygon_col = 2

    edge_length = 2

    def __init__(self, n, m, p, use_outlines, angle_index=0):
        """Initialise object

        This will create a cupola prismatoid with bases {n/m} sides {n/p}.

        n: amount of vertices in the base and the sides
        m: the number m in the base {n/m} with m < n/2.
        p: the number p in {n/p} used for the sides
        use_outlines: if set to True then the {n/m}-gram(s) will be replaced by a polygon following
            the outline. If set to False the polygon will follow the n edges, which might not be
            shown well in a 3D player, e.g. holes might appear at parts that have even coverage.
        angle_index: If more than one triangle angle is found, use the specified index.
        """
        if n < 3:
            raise ValueError("n must be bigger than 3")

        # any data for this pseudo cupolaic prismatoid
        self._cupola_data = Object()
        self._cupola_data.n = n
        self._cupola_data.m = m
        self._cupola_data.p = p
        self._cupola_data.base_at_z = 0
        self._cupola_data.use_index = angle_index
        self._cupola_data.use_outlines = use_outlines

        if use_outlines:
            LOGGER.info("Using outlines, OFF files will not load in Stella!")
        else:
            LOGGER.info("Not using outlines, OFF files might show holes in Orbitit!")

        vertices = []
        faces = []
        col_i = []

        def add_face(vs, col):
            if use_outlines:
                face = geom_3d.Face(vs)
                vs = face.outline.vs
            offset = len(vertices)
            vertices.extend(vs)
            faces.append([offset + i for i in range(len(vs))])
            col_i.append(col)

        def add_face_list(vs, col):
            for list_of_v in vs:
                add_face(list_of_v, col)


        # TODO: perhaps we should make this a compound shape
        base_vs = get_polygon(n, m)
        add_face_list(base_vs, self.base_col)

        side_vs = get_polygon(n, p)
        side_sub = side_vs[0]
        len_sub = len(side_sub)
        assert len_sub > 2, "Digons not supported here"
        delta = (side_sub[2 % len_sub] - side_sub[0]).norm()
        triangle = self._get_triangle(base_vs[0], delta)
        for i in range(n):
            rotate = geomtypes.Rot3(axis=geom_3d.vec(0, 0, 1), angle=i * TWO_PI / n)
            add_face([rotate * v for v in triangle], self.triangle_col)

        # + 1: because the base has index 0
        triangle_next = faces[1 + m]
        triangle_next_top = vertices[triangle_next[2]]

        side_vs = self._attach_side(
            side_vs,
            [0, 1],
            [triangle[1], triangle[2]],
            -1,
            triangle_next_top,
        )

        for i in range(n):
            rotate = geomtypes.Rot3(axis=geom_3d.vec(0, 0, 1), angle=i * TWO_PI / n)
            add_face_list(
                [
                    [rotate * v for v in v_list]
                    for v_list in side_vs
                ],
                self.side_polygon_col,
            )

        super().__init__(
            vertices,
            faces,
            colors=(cols, col_i),
            name=f"Cupola {n}/{m} with {n}/{p}",
        )

    def _get_triangle(self, base_vs, delta):
        """Attach equilateral triangles to each edge of the base.

        These are folded up so that a n/p polygon fits in between.
        base_vs: the vertices of the base polygon is order.
        delta: the required distance between to neighbouring triangles after being folded up.
        """
        v0 = base_vs[0]
        v1 = base_vs[1]
        # The point between these is also in the direction of v3
        centre = (v0 + v1) / 2
        # equilateral triangle
        v2 = centre + sqrt(3) * centre.normalise()
        vs = [v0, v1, v2]

        # the reflection to use to obtain the triangle next to vs[0].
        refl = geomtypes.Refl3(normal=base_vs[1] - base_vs[-1])

        def how_close_when_folded(fold_angle):
            """Return the difference between angle that is left and the angle that it needed.

            Only positive values are returned.
            """
            return abs(self._get_diff_for_fold(vs, fold_angle, refl) - delta)

        def _try_angles(method="Powell"):
            """Try angles"""
            solutions = []
            no_of_guesses = 10
            angle_step = TWO_PI /  no_of_guesses
            for i in range(no_of_guesses):
                initial_guess = i * angle_step
                result = minimize(how_close_when_folded, initial_guess, method=method)
                if result.success:
                    for fold_angle in result.x.tolist():
                        fold_angle = (fold_angle + pi) % TWO_PI - pi
                        with geomtypes.FloatHandler(4):
                            if geomtypes.FloatHandler.eq(fold_angle, 0):
                                continue
                            if geomtypes.FloatHandler.eq(fold_angle, pi):
                                continue
                            if geomtypes.FloatHandler.eq(fold_angle, -pi):
                                continue
                        already_found = False
                        for solution in solutions:
                            with geomtypes.FloatHandler(self.exp_tol_eq_float):
                                if geomtypes.FloatHandler.eq(fold_angle, solution):
                                    already_found = True
                                    break
                        if not already_found:
                            solutions.append(fold_angle)
            return solutions

        # Powell seems to work best, it finds the minimum and the solution has high precision
        # However the vertices don't fit for that one.
        solutions = _try_angles()
        if solutions:
            if len(solutions) > 1:
                LOGGER.warning("Found %d solutions: %s!", len(solutions), solutions)
                for angle in solutions:
                    LOGGER.info("  %.2f degrees", angle * 180 / pi)
            use_angle = solutions[
                self._cupola_data.use_index if self._cupola_data.use_index < len(solutions) else 0
            ]
            LOGGER.info("Using fold angle %.2f degrees", use_angle * 180 / pi)
            return [v0, v1, self._get_top_for_fold(vs, use_angle)]
        raise ValueError("Couldn't find a triangle fold to fit polygon")

    # TODO: make static?
    def _get_top_for_fold(self, vs, angle):
        """Return the top of the equilateral triangle when folded.

        vs: the three vertices of the triangle. It is assumed that the triangle is lying the in z=0
            plane. The triangle will be folded around the edge connecting vs[0] and vs[1].
        angle: angle in radians to fold up the triangle

        return: the coordinate of the new top.
        """
        fold = geomtypes.Rot3(axis=vs[0] - vs[1], angle=angle)
        # Translate so that the centren of the edge vs[0] -- vs[1] ends up at the origin
        centre = (vs[0] + vs[1]) / 2
        vs2 = vs[2] - centre
        # fold and translate back
        return fold * vs2 + centre

    def _get_diff_for_fold(self, vs, angle, refl):
        """Return the distance between the triangles when folding up.

        vs: the three vertices of the triangle. It is assumed that the triangle is lying the in z=0
            plane. The triangle will be folded around the edge connecting vs[0] and vs[1].
        angle: angle in radians to fold up the triangle
        refl: the reflection to use to obtain the triangle next to vs[0].

        return: the difference between the top of two neighbouring triangles.
        """
        new_top = self._get_top_for_fold(vs, angle)
        other_top = refl * new_top
        return (new_top - other_top).norm()

    # TODO: perhaps this method can be static
    # A more general method would make sure that the distance to vertex gets a certain value, but
    # then minimize must be used (and there will be more than one solution)
    def _attach_side(self, side_vs, side_edge, attach_to, side_vertex, vertex, sub_index=0):
        """Attach an edge of side_vs to a specified edge and fold into correct position.

        side_vs: the vertices of the polygon to attach. This is a list of list to support {9/3} e.g.
        side_edge: two indices in side_vs[sub_index] specifying an edge to attach
        attach_to: a list for two coordinates, where side_edge[0] shall be attached to the first
            coordinate in the list and side_edge[1] to the other one.
        side_vertex: an index in side_vs[sub_index] different from the ones in side_edge. The face
            side_vs will be folded over side_edge so that side_vertex will have a certain distance
            to vertex
        vertex: a coordinate to which side_vertex needs to be attaced
        sub_index: the list of vertices in side_vs to use when attaching

        return: a new array of side_vs for the transformed face
        """
        side_sub = side_vs[sub_index]
        i_v0 = side_edge[0]
        i_v1 = side_edge[1]
        translate = attach_to[0] - side_sub[i_v0]
        vs = [[translate + v for v in v_list] for v_list in side_vs]
        # Rotate to get v1 attached to attach_to[1]
        side_sub = vs[sub_index]
        vec0 = side_sub[i_v1] - side_sub[i_v0]
        vec0 = vec0.normalise()
        vec1 = attach_to[1] - attach_to[0]
        vec1 = vec1.normalise()
        angle = acos(vec0 * vec1)
        axis = vec0.cross(vec1)
        rotate = geomtypes.Rot3(axis=axis, angle=angle)
        make_origin = side_sub[i_v0]
        vs = [[v - make_origin for v in v_list] for v_list in vs]
        vs = [[rotate * v for v in v_list] for v_list in vs]
        vs = [[v + make_origin for v in v_list] for v_list in vs]
        # Now rotate around side_edge to attach side_vertex to vertex

        # TODO: break out code and share (the part the translates and rotates)
        side_sub = vs[sub_index]
        edge = geom_3d.Line3D(side_sub[0], p1=side_sub[1])
        make_origin = edge.project(vertex)
        edge_to_side_vertex = edge.to_point(side_sub[side_vertex]).normalise()
        edge_to_goal_vertex = edge.to_point(vertex).normalise()
        angle = acos(edge_to_side_vertex * edge_to_goal_vertex)
        # try the angle:
        axis = edge.v.normalise()
        rotate = geomtypes.Rot3(axis=axis, angle=angle)
        v = side_sub[side_vertex] - make_origin
        v = rotate * v
        v = v + make_origin
        with geomtypes.FloatHandler(self.exp_tol_eq_float):
            if not v == vertex:
                rotate = geomtypes.Rot3(axis=axis, angle=-angle)
                v = side_sub[side_vertex] - make_origin
                v = rotate * v
                v = v + make_origin
                with geomtypes.FloatHandler(self.exp_tol_eq_float):
                    assert v == vertex, "Expected to have mapped v on vertex by now"
        vs = [[v - make_origin for v in v_list] for v_list in vs]
        vs = [[rotate * v for v in v_list] for v_list in vs]
        vs = [[v + make_origin for v in v_list] for v_list in vs]

        return vs

    def base_use_outlines(self):
        """Replace the n-grams by there outlines to prevent holes."""
        face_index = 0
        self.replace_face_by_outline(face_index, self.exp_tol_eq_float)
        return
        for base in self.bases:
            for _ in range(base.no_of_compounds):
                self.replace_face_by_outline(face_index, self.exp_tol_eq_float)
                face_index += 1


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("n", type=int, help="Number of vertices of the {n/m}-polygon.")
    parser.add_argument(
        "m",
        type=int,
        help="Vertex offset of the primary base polygon {n/m}. Use m > n/2."
    )
    parser.add_argument(
        "p",
        type=int,
        help="Vertex offset of the secondary base polygon {n/p}. Use m > n/2."
    )
    parser.add_argument(
        "-i", "--angle_index",
        type=int,
        default=0,
        help="If more than one triangle angle is found, use the specified index."
    )
    parser.add_argument(
        "-b", "--file_base_name",
        default="cupola_",
        help="A header to name of the file. This will be used as a base for the OFF file. "
        "It will be appended by n_m__n_p.off."
    )
    parser.add_argument(
        "-t", "--file_tail_name",
        default="",
        help="A string to append to name of the file. This will be appended to 'n_m__n_p'."
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

    shape = Cupola(ARGS.n, ARGS.m, ARGS.p, not ARGS.allow_holes, angle_index=ARGS.angle_index)

    if ARGS.x_rotate:
        shape.transform(
            geomtypes.Rot3(
                angle=geom_3d.DEG2RAD * ARGS.x_rotate,
                axis=geomtypes.Vec3([1, 0, 0]),
            )
        )

    model = f"{ARGS.n}_{ARGS.m}__{ARGS.n}_{ARGS.p}"
    filepath = Path(ARGS.out_dir) / f"{ARGS.file_base_name}{model}{ARGS.file_tail_name}.off"
    if not ARGS.overwrite and filepath.is_file():
        yes_or_no = input(f"{filepath} exists. Overwrite? y/N\n")
        if not yes_or_no or yes_or_no.lower()[0] != "y":
            LOGGER.warning("No overwrite requested; bailing out")
            sys.exit(1)

    with open(filepath, "w") as fd:
        minimized_shape = shape.clean_shape(shape.exp_tol_eq_float)
        fd.write(minimized_shape.to_off())
        LOGGER.info("Written %s", filepath)
