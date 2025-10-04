"""Script to generate a cupalaic prismatoid based on a {n/m}-gram"""

import argparse
import logging
from math import acos, cos, gcd, pi, sin, sqrt
import os
from pathlib import Path
import sys

from scipy.optimize import fsolve, minimize
from orbitit.colors import STD_COLORS as cols
from orbitit import geom_3d, geomtypes

# TODO: update text below
DESCRIPTION = """Generate an off file for {n/p} base and {m/q} sides cupolaic shapes.

The scipt tries to put together two sides consisting of {m/q} polygons sharing an edge and it will
fit a equilateral triangle on each side of the edge. The result is a shape that might look like a
half-hip roof.

These roofs are attached to an {n/p} polygon where a free triangle edge is attached to an edge of
the {n/p} polygon. The dihedral angle between these faces is adjusted in such a way that the
distance between the next vertices has the same length as one edge. This would mean that another
equilateral triangle would fit.

The hope is that when the half-hip root is rotated n times that the triangles of the half-hip roof
coincide with that extra one and that one side of the roof coincide with the other half of a rotated
one.

Note that not all parameters have solutions.
"""
# TODO: raise ValueError when n > 11
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s - %(message)s",
    datefmt="%m-%d %H:%M",
)
LOGGER = logging.getLogger("pseudo copulaic prismatoid generator")
TWO_PI = 2 * pi


def get_polygon(n, m, edge_length=2, rotate=True):
    """Return a the base shape.

    The vertices are the vertices of an n-gon z = 0. The edges of the n/m gram will
    scaled to edge_length.

    edge_length: the required edge length
    rotate: if this validates to True, then the first edge will be parallel to the Y-axis. Otherwise
        the first vertex will be on the x-axis.

    return: a list with vertices (geomtypes.Vec3)
    """
    # E.g. {7/2}
    #                           ^ X
    #              v4           |
    #       v1            v0    |
    #                           |
    #                           +
    #     v5                v3  |
    #                           |
    #                           |
    #          v2        v6     |
    # Y <-----------+-----------
    no_of_compounds = gcd(n, m)
    no_of_vs_x_gram = n // no_of_compounds
    if no_of_compounds > 1:
        LOGGER.info("{%i/%i} has a common divider", n, m)
    offset = m * pi / n if rotate else 0
    vs = [
        [
            geom_3d.vec(
                cos((m * i + j) * TWO_PI / n - offset),
                sin((m * i + j) * TWO_PI / n - offset),
                0,
            )
            for i in range(no_of_vs_x_gram)
        ]
        for j in range(no_of_compounds)
    ]

    diagonal = (vs[0][0] - vs[0][1]).norm()
    scale = edge_length / diagonal
    return [[scale * v for v in v_list] for v_list in vs]


def translate_list_of_vs(list_of_vs, vector):
    """Translate a list of list of Vec3 points."""
    return [[v + vector for v in vs] for vs in list_of_vs]


def fold_to_fit(
    fold_func,
    distance_func,
    precision=7,
    angle_domain=None,
    angle_step=TWO_PI / 10,
):
    """Fold a face over one of its edges until a condition is met.

    One example can be to fold a face until the distance to a vertex has a certain value.

    fold_func: function accepting the fold angle as parameter and returns the vs for this fold.
    distance_func: a function that returns the lowest value for when the condition is met. The
        function will accept one parameter: angle. For the example, given above one could write a
        function that returns the absolute difference between the actual distance and the required
        distance.
    precision: the amount of decimals for two floats, e.g. angles, to be interpreted as the same,
    angle_domain: a list with minimal and maximal angle in radians. If not specified, then all
        angles between 0 and 2π will be used as initial guess for a fold angle.
    angle_step: the increase in angle (in radians) when trying different angles in the domain as
        initial guess.

    return: a list of solutions. Each solution consists of a dictionary with keys "angle" and "vs",
        where "angle" gives the angle in radians and "vs" contains the adjusted vertices from the
        input parameter "vs".
    """
    solutions = []
    if angle_domain is None:
        angle_domain = [0, TWO_PI]
    initial_guess = angle_domain[0]
    while not initial_guess > angle_domain[1]:
        results = fsolve(distance_func, initial_guess)
        initial_guess += angle_step
        for fold_angle in results:
            # map angle to domain [-π, π)
            fold_angle = (fold_angle + pi) % TWO_PI - pi
            # Don't return angles where faces end up in the same plane
            with geomtypes.FloatHandler(precision):
                if geomtypes.FloatHandler.eq(fold_angle, 0):
                    continue
                if geomtypes.FloatHandler.eq(fold_angle, pi):
                    continue
                if geomtypes.FloatHandler.eq(fold_angle, -pi):
                    continue
                # For some reason sometimes the solutions isn't really a solutions,
                # e.g. for {7/2} I got -0.9112721567972777 with distance ~= 0.28400
                if geomtypes.FloatHandler.ne(distance_func(fold_angle), 0):
                    continue
            # only add new solutions
            already_found = False
            for solution in solutions:
                with geomtypes.FloatHandler(precision):
                    if geomtypes.FloatHandler.eq(fold_angle, solution["angle"]):
                        already_found = True
                        break
            if not already_found:
                LOGGER.debug(
                    "Solution found for angle %0.2f (minimized distance %0.2f)",
                    geom_3d.RAD2DEG * fold_angle,
                    distance_func(fold_angle),
                )
                solutions.append(
                    {
                        "angle": fold_angle,
                        "vs": fold_func(fold_angle),
                    }
                )
    return solutions


class Object:
    """Just a class to be able to add attributes to an object."""


class Polygram:
    """A polygram specification {n/m}."""

    def __init__(self, n, m):
        "An {n/m} polygram."
        self.n = n
        self.m = m


class RoofTop(geom_3d.SimpleShape):
    """class for creating cupolaic objects with roof top shaped sides."""

    # When taking the outline of a concave polygon, use this exponent to decide whether two floats
    # are equal
    exp_tol_eq_float = 7

    # colour indices from orbitit.colors.STD_COLORS
    base_col = 3
    triangle_col = [1, 6, 7]
    side_polygon_col = 2

    edge_length = 2

    def __init__(
        self,
        base: Polygram,
        side: Polygram,
        use_outlines,
        args,
    ):
        """Initialise object

        This will create a cupolaic prismatoid with bases {n/m} sides {n/p}.

        base: the {n/m} polygram used at the base.
        side: the {n/m} polygram used at the sides.
        use_outlines: if set to True then the {n/m}-gram(s) will be replaced by a polygon following
            the outline. If set to False the polygon will follow the n edges, which might not be
            shown well in a 3D player, e.g. holes might appear at parts that have even coverage.
        args: a Namespace object with the following fields:
            angle_index: If more than one triangle angle is found, use the specified index.
            use_whole_roof: If True then the {n/p} polygon that use to construct a half-hip roof
                is kept in the final shape. Otherwise only one side of the roof is kept.
            add_extra_triangles: If True then the extra equilateral triangles are added. These extra
                triangles are from attaching the roofs to the base polygon.
            no_base: if set then the base isn't added, which can make sense when add_extra_triangles
                is set.
        """
        if base.n < 3:
            raise ValueError("n must be bigger than 3")
        if side.n < 3:
            raise ValueError("n must be bigger than 3")

        # any data for this shape
        self._shape_data = Object()
        self._shape_data.base = base
        self._shape_data.side = side
        self._shape_data.base_at_z = 0
        self._shape_data.use_outlines = use_outlines
        self._shape_data.use_index = args.angle_index
        self._shape_data.use_whole_roof = args.use_whole_roof
        self._shape_data.add_extra_triangles = args.add_extra_triangles
        self._shape_data.add_base = not args.no_base

        if use_outlines:
            LOGGER.info("Using outlines, OFF files will not load in Stella!")
        else:
            LOGGER.info("Not using outlines, OFF files might show holes in Orbitit!")

        vertices = []
        faces = []
        col_i = []

        def add_face(vs, col):
            """Add a face to the vertices, faces and col_i variables."""
            if use_outlines and len(vs) > 3:
                face = geom_3d.Face(vs)
                vs = face.outline.vs
            offset = len(vertices)
            vertices.extend(vs)
            faces.append([offset + i for i in range(len(vs))])
            col_i.append(col)

        def add_face_list(vs, col):
            """Add a list of faces to the vertices, faces and col_i variables."""
            for list_of_v in vs:
                add_face(list_of_v, col)

        # TODO: move this into a method:
        #######
        #  1  #
        #######
        # Try to make a "half-hip roof" of two sides, using an equilateral triangle as hip.
        side_vs = get_polygon(self._shape_data.side.n, self._shape_data.side.m)
        side_sub = side_vs[0]
        len_sub = len(side_sub)
        assert len_sub > 2, "Digons not supported here"

        # translate origin to centre of first edge to construct the roof
        make_origin = (side_sub[0] + side_sub[1]) / 2
        side_vs_t = translate_list_of_vs(side_vs, -make_origin)
        fold_axis = side_vs_t[0][1] - side_vs_t[0][0]

        def fold_side(angle):
            """Rotate side_vs_t around first edge by alpha radians and return new vs."""
            rotate = geomtypes.Rot3(axis=fold_axis, angle=angle)
            return [[rotate * v for v in v_list] for v_list in side_vs_t]

        def get_distance(angle):
            """Get how far off the distand is from the required distance of 2.

            The distance that is used is the distance between third vertex and rotated third
            vertex.
            """
            with geomtypes.FloatHandler(4):
                if geomtypes.FloatHandler.ne(angle, 0):
                    pass
            v0 = side_vs_t[0][2]
            v1 = geomtypes.Rot3(axis=fold_axis, angle=angle) * v0
            distance = (v1 - v0).norm()
            abs_diff = abs(distance - 2)
            return abs_diff

        LOGGER.info(
            "Looking for fold angle to create half-hip roof of {%d/%d} (and triangle)",
            self._shape_data.side.n,
            self._shape_data.side.m,
        )
        solutions = fold_to_fit(
            fold_side,
            get_distance,
            precision=self.exp_tol_eq_float,
        )
        if not solutions:
            LOGGER.error("No angle found to create a half-hip roof")
            sys.exit(1)

        index = 0
        if len(solutions) > 1:
            LOGGER.info(
                "Found %d solutions, will use fold angle alpha = %.2f degrees",
                len(solutions),
                solutions[index]["angle"] * 180 / pi,
            )

        # For balancing the two sides: it is better to rotate both: -angle/2 and angle/2
        # The the rotationa axis for the orbit shape will be the z-axis
        half_angle = solutions[index]["angle"] / 2

        # Add translated origin and make the fold incl. the origin
        side_vs_t.append([-make_origin])
        side_vs_t0 = fold_side(-half_angle)
        side_vs_t1 = fold_side(half_angle)

        # remove origin:
        new_origin = -(side_vs_t0[-1][0] + side_vs_t1[-1][0]) / 2
        side_vs_t0 = side_vs_t0[:-1]
        side_vs_t1 = side_vs_t1[:-1]

        # translate back
        both_sides = translate_list_of_vs(side_vs_t0, new_origin)
        next_side_i = len(both_sides)
        both_sides.extend(translate_list_of_vs(side_vs_t1, new_origin))

        # add one edge, that needs to attached to {n/p}
        both_sides.append([both_sides[0][2], both_sides[next_side_i][2]])

        #######
        #  2  #
        #######
        # Attach to the base {n/m}
        top_vs = get_polygon(self._shape_data.base.n, self._shape_data.base.m)
        both_sides = self._attach_edges(both_sides, [0, 1], top_vs[0][:2], sub_index=-1)

        # remove that edge again
        del both_sides[-1]

        # Now fold around that edge to form more triangles.
        LOGGER.info(
            "Fitting half-hip roof to top {%d/%d} polygon",
            self._shape_data.base.n,
            self._shape_data.base.m,
        )
        axis_direction = top_vs[0][0] - top_vs[0][1]
        axis_through = top_vs[0][0]
        fold_vertex = both_sides[0][3]

        def fold_result(alpha):
            transform = geomtypes.Rot3NonCentered(axis_direction, axis_through, alpha)
            new_vec = transform * fold_vertex
            return abs(2 - (top_vs[0][-1] - new_vec).norm())

        solutions = []
        no_of_steps = 20
        for step in range(no_of_steps):
            alpha = -pi + step * TWO_PI / no_of_steps

            result = minimize(fold_result, alpha, method="Powell")
            if result.success:
                if geomtypes.FloatHandler.ne(result.fun, 0):
                    continue
                fold_angle = result.x[0]
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

        if not len(solutions):
            LOGGER.error(
                "No solutions found to fit {%d/%d} roofs to {%d/%d}",
                self._shape_data.side.n,
                self._shape_data.side.m,
                self._shape_data.base.n,
                self._shape_data.base.m,
            )
            raise ValueError("Try with other parameters")

        LOGGER.info("Found %d solutions:", len(solutions))
        for angle in solutions:
            LOGGER.info("  %0.2f degrees", geom_3d.RAD2DEG * angle)
        # TODO: check length
        using = (
            self._shape_data.use_index
            if self._shape_data.use_index < len(solutions)
            else 0
        )
        angle = solutions[using]
        LOGGER.info("Will use %0.2f degrees (index %d)", geom_3d.RAD2DEG * angle, using)
        transform = geomtypes.Rot3NonCentered(axis_direction, axis_through, angle)

        both_sides = [[transform * v for v in face] for face in both_sides]
        extra_triangle = [
            both_sides[0][3],
            top_vs[0][0],
            top_vs[0][-1],
        ]

        #######
        #  3  #
        #######
        # Put these together
        if self._shape_data.add_base:
            for face in top_vs:
                face.reverse()
            add_face_list(top_vs, self.base_col)
        if self._shape_data.use_whole_roof:
            to_orbit = both_sides
        else:
            to_orbit = both_sides[next_side_i:]
        triangles = [
            [both_sides[0][1], both_sides[0][2], both_sides[next_side_i][2]],
            [both_sides[0][0], both_sides[next_side_i][-1], both_sides[0][-1]],
            extra_triangle,
        ]
        # for face in to_orbit:
        #    face.reverse()
        angle_step = TWO_PI / self._shape_data.base.n
        # Note: it is less efficient to have separate loops here, but then the faces are sorted
        for i in range(self._shape_data.base.n):
            transform = geomtypes.Rot3(
                axis=geomtypes.Vec3([0, 0, 1]), angle=i * angle_step
            )
            add_face_list(
                [[transform * v for v in face] for face in to_orbit],
                self.side_polygon_col,
            )
        no_of_triangles = 3 if self._shape_data.add_extra_triangles else 2
        for j in range(no_of_triangles):
            for i in range(self._shape_data.base.n):
                transform = geomtypes.Rot3(
                    axis=geomtypes.Vec3([0, 0, 1]), angle=i * angle_step
                )
                add_face(
                    [transform * v for v in triangles[j]],
                    self.triangle_col[j],
                )

        super().__init__(
            vertices,
            faces,
            colors=(cols, col_i),
            name=(
                f"Roof top polyhedron {self._shape_data.base.n}/{self._shape_data.base.n} with "
                f"{self._shape_data.side.n}/{self._shape_data.side.n}"
            ),
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

        solutions = fold_to_fit(
            lambda a: [v0, v1, self._get_top_for_fold(vs, a)],
            lambda a: abs(self._get_diff_for_fold(vs, a, refl) - delta),
            precision=self.exp_tol_eq_float,
        )

        if solutions:
            if len(solutions) > 1:
                LOGGER.warning("Found %d solutions!", len(solutions))
                for s in solutions:
                    LOGGER.info("  %.2f degrees", s["angle"] * 180 / pi)
            solution = solutions[
                self._shape_data.use_index
                if self._shape_data.use_index < len(solutions)
                else 0
            ]
            LOGGER.info("Using fold angle %.2f degrees", solution["angle"] * 180 / pi)
            return solution["vs"]
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
    # TODO: this method should return a 4D matrix, then you don't need to
    def _attach_edges(self, faces, edge, attach_to, sub_index=0):
        """Attach an edge of faces to a specified edge

        faces: the vertices of the polygon to attach. This is a list of list to support {9/3} e.g.
        edge: two indices in faces[sub_index] specifying an edge to attach
        attach_to: a list for two coordinates, where vertex index edge[0] shall be attached to the
            first coordinate in the list and edge[1] to the other one.
        sub_index: the list of vertices in faces to use when attaching

        return: a new array of faces for the transformed face
        """
        # Move this first part out: attach edge to edge
        side_sub = faces[sub_index]
        i_v0 = edge[0]
        i_v1 = edge[1]
        vs = translate_list_of_vs(faces, attach_to[0] - side_sub[i_v0])
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
        vs = translate_list_of_vs(vs, -make_origin)
        vs = [[rotate * v for v in v_list] for v_list in vs]
        vs = translate_list_of_vs(vs, make_origin)
        return vs

    def fold_until(
        self,
        vs,
        axis: geomtypes.Vec3,
        axis_through: geomtypes.Vec3,
        init_angle,
        to_minimize,
    ):
        """Fold vertices around an axis to minimize a certain value.

        vs: an array with 3D vertices, must be of Vec3 type
        axis: a direction vector.
        axis_through: a point on the axis
        init_angle: start with this fold angle, which should be close to the solution
        to_minimize: a function that accepts the folded vs and returns a values. It is assumed
            that for the fold angle that is being looked for the returned values is the lowest.

        Return: a 4D matrix that expresses the transform using homogeneous coordinates.
        """

    def fold_to_attach_vertices(self, vs, sub_index, side_vertex, vertex):
        """Old code to fold a side along a (shared) edge until vertices coincide."""
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
    parser.add_argument(
        "base_n",
        type=int,
        help="Number of vertices of the base {n/m}-polygon.",
    )
    parser.add_argument(
        "base_m",
        type=int,
        help="Vertex offset of the base {n/m}-polygon. Use m > n/2.",
    )
    parser.add_argument(
        "side_n",
        type=int,
        help="Number of vertices of the side {n/m}-polygon.",
    )
    parser.add_argument(
        "side_m",
        type=int,
        help="Vertex offset of the side {n/m}-polygon. Use m > n/2.",
    )
    parser.add_argument(
        "-i",
        "--angle_index",
        type=int,
        default=0,
        help="If more than one triangle angle is found, use the specified index.",
    )
    parser.add_argument(
        "-b",
        "--file_base_name",
        default="roof_top_",
        help="A header to name of the file. This will be used as a base for the OFF file. "
        "It will be appended by n_m__n_p.off.",
    )
    parser.add_argument(
        "-t",
        "--file_tail_name",
        default="",
        help="A string to append to name of the file. This will be appended to 'n_m__n_p'.",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        default=".",
        help="path to directory to save the resulting OFF file(s).",
    )
    parser.add_argument(
        "-H",
        "--allow_holes",
        action="store_true",
        help="If specified a {n/m} polygon will saved using n vertices and edges, which results "
        "in holes in orbitit for the parts that have even coverage due to the stencil buffer. "
        "If not specified the {n/m} polygons and the crossed rectangles will be replaced by their "
        "outline, which results in OFF files where edges are broken and hence the will not have "
        "an even amount of faces joining in each edge, which might result in warnings or errors "
        "for some 3D programs.",
    )
    parser.add_argument(
        "-w",
        "--overwrite",
        action="store_true",
        help="If specified an existing file will be overwritten without asking. Otherwise the "
        "the script will ask interactively whether to overwrite an existing file.",
    )
    parser.add_argument(
        "--use_whole_roof",
        action="store_true",
        help="If specified the {n/p} polygon that was used to construct a half-hip roof is kept in "
        "the final shape. Otherwise only one side of the roof is kept.",
    )
    parser.add_argument(
        "--add_extra_triangles",
        action="store_true",
        help="When attaching the sides to the base, the program will try to used a dihedral angle "
        "between the triangle of the roof and the base so that the 'next' vertex of the side "
        "polygon is an edge away from the next vertex in of the base. The next vertex of the side "
        "polygon is seen from the top of the roof. This means that equilateral triangles are "
        "formed by the side and base edges to the 'next' vertex and an edge between these 'next' "
        "vertices. If you set this option these triangles are added to the shape. Note that that "
        "means that each edge of the base will join three faces: the triangles from the roof, the "
        "the base and the extra triangles. Therefore it makes sense to remove the base in this "
        "case.",
    )
    parser.add_argument(
        "--no_base",
        action="store_true",
        help="If specified the {n/p} polygon that is the base will not be added to the final "
        "shape. This is useful when --add_extra_triangles is set.",
    )
    parser.add_argument(
        "-x",
        "--x-rotate",
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

    shape = RoofTop(
        Polygram(ARGS.base_n, ARGS.base_m),
        Polygram(ARGS.side_n, ARGS.side_m),
        not ARGS.allow_holes,
        ARGS,
    )
    shape.transform(geomtypes.Roty(angle=-pi / 2))

    sum_of_vs = geomtypes.Vec3([0, 0, 0])
    for v in shape.vs:
        sum_of_vs += v
    shape.translate(-sum_of_vs / len(shape.vs))

    if ARGS.x_rotate:
        shape.transform(
            geomtypes.Rot3(
                angle=geom_3d.DEG2RAD * ARGS.x_rotate,
                axis=geomtypes.Vec3([1, 0, 0]),
            )
        )

    model = (
        f"{ARGS.base_n}_{ARGS.base_m}__{ARGS.side_n}_{ARGS.side_m}_{ARGS.angle_index}"
    )
    filepath = (
        Path(ARGS.out_dir) / f"{ARGS.file_base_name}{model}{ARGS.file_tail_name}.off"
    )
    if not ARGS.overwrite and filepath.is_file():
        yes_or_no = input(f"{filepath} exists. Overwrite? y/N\n")
        if not yes_or_no or yes_or_no.lower()[0] != "y":
            LOGGER.warning("No overwrite requested; bailing out")
            sys.exit(1)

    with open(filepath, "w") as fd:
        minimized_shape = shape.clean_shape(shape.exp_tol_eq_float)
        fd.write(minimized_shape.to_off())
        LOGGER.info("Written %s", filepath)
