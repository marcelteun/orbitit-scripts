"""Script to generate JSON files for different kinds of tetartoid polyhedra

A tetartoid is a dodecahedron of 5 pentagons. The whole model has at least the rotation symmetries
of a tetrahedron.
"""
# pylint: disable=too-many-lines
import argparse
from collections.abc import Callable
from enum import IntEnum, auto
import json
import logging
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

class EdgeCat(IntEnum):
    """Defines the category of edge properties that a tetartoid can have."""
    GENERAL = auto()  # there are three different edge lengths
    E0_EQ_E1_2 = auto()  # edge 0 has the same length as edge 1 and 2
    E0_EQ_E3_4 = auto()  # edge 0 has the same length as edge 3 and 4
    E1_2_EQ_E3_4 = auto()  # edge 1 and 2 have the same length as edge 3 and 4
    ALL_EQ = auto()  # all edges have equal length

class AngleCat(IntEnum):
    """Defines the category of angle properties that a tetartoid can have.

    With angle it is meant the angle in a face between the edges of a vertex.
    """
    GENERAL = auto()  # there a five difference angles
    ONE_EQ_PAIR = auto()  # there is one pair with the same angle
    TWO_EQ_PAIRS = auto()  # there are two pairs with the same angle
    EQ_PAIR_AND_TRIPLE = auto()  # there is one pairs and a triplet with the same angle
    EQ_QUARTET = auto()  # there are four angles the same angle
    ALL_EQ = auto()  # all angles are the same

LOGGER = logging.getLogger("tetartoid")
logging.basicConfig(
    format="%(levelname)s: %(message)s",
    level=logging.INFO,
)
DESCR = """Generate Orbitit JSON files for a tetartoid polyhedron.

A tetartoid is a dodecahedron of 5 pentagons. The whole model has at least the rotation symmetries
of a tetrahedron.

The definition that is used here is from https://en.wikipedia.org/wiki/Dodecahedron#Tetartoid
with the exeption of the first requirement: 0 <= a <= b <= c which is mainly to prevent pentagrams.

The JSON files are in the format of Orbitit describing orbit Shapes. The JSON file defines one
face and the final symmetry (using E as the stabiliser symmetry).
"""

def tetar_n(a, b, c):
    """Get the value of 'n' for a tetartoid defined by a, b, c."""
    return a**2 * c - b * c**2

def tetar_d1(a, b, c):
    """Get the value of 'd1' for a tetartoid defined by a, b, c."""
    return a * (a - b + c) + b * (b - 2 * c)

def tetar_d2(a, b, c):
    """Get the value of 'd2' for a tetartoid defined by a, b, c."""
    return a * (a + b - c) + b * (b - 2 * c)

def tetar_face(a, b, c, e1, e2):
    """From the five input numbers return one face of a tetartoid.

    The values of a, b, and c are the standard values mention on the Wikipedia page.
    e1 = n / d1
    e2 = n / d2
    """
    vertices = [
        [a, b, c],
        [-a, -b, c],
        [-e1, -e1, e1],
        [-c, -a, b],
        [-e2, e2, e2],
    ]
    vertices = np.array(vertices)
    return vertices


def get_edges(face):
    f_len = len(face)
    return [
        face[i] - face[(i + 1) % f_len]
        for i in range(f_len)
    ]


def get_edge_lengths(edges) -> list[float]:
    return np.linalg.norm(edges, axis=1).tolist()


def get_edge_angles(face) -> list[float]:
    """Return angles between sides in a face in radians.

    Note that for edges with a length that is almost 0 the angles doesn't say so much.
    """
    edges = get_edges(face)
    edge_lengths = get_edge_lengths(edges)
    no_of_vs = len(edges)
    result = [
        np.arccos(
            np.clip(
                np.dot(edges[i], edges[(i + 1) % no_of_vs]) / (
                    edge_lengths[i] * edge_lengths[(i + 1) % no_of_vs]
                ), -1, 1
            )
        )
        for i in range(no_of_vs)
    ]
    return result


class Tetartoid():
    # Increased compared to numpy for tetrahedron
    rtol = 1e-4
    atol = 1e-8
    related = None

    def __init__(self, a, b, c, name="A4"):
        """
        name: will be used as Shape name
        """
        # Original requirement: will only allow for convex tetartoids
        # assert a <= b <= c
        self.a = a
        self.b = b
        self.c = c
        self.name = name

        self.shape_props = {
            "class": "orbit_shape_v2",
            "data": {
                "version": 2,
                "name": name if name else "A4",
                "base": {
                    # vs is added in the code below
                    "fs": [[0, 1, 2, 3, 4]],
                },
                "final_sym": {
                    "class": "A4",
                    "data": {
                        "generator": {
                            "o2axis0": [1, 0, 0],
                            "o2axis1": [0, 1, 0],
                        }
                    },
                },
                "stab_sym": {"class": "E", "data": {}},
                "cols": [
                    [64, 104, 224],
                    [153, 204, 49],
                    [220, 159, 220],
                    [254, 214, 0],
                ],
                "no_of_cols": 4,
                "col_sym": "C3",
                "col_alt": 2,  # note 0 based
            }
        }
        assert self.n * self.d1 * self.d2 != 0

    @property
    def abc(self):
        return self.a, self.b, self.c

    @property
    def n(self):
        return tetar_n(self.a, self.b, self.c)

    @property
    def d1(self):
        return tetar_d1(self.a, self.b, self.c)

    @property
    def d2(self):
        return tetar_d2(self.a, self.b, self.c)

    @property
    def face(self):
        """Return one face of a tetartoid."""
        e1 = self.n / self.d1
        e2 = self.n / self.d2
        return tetar_face(self.a, self.b, self.c, e1, e2)

    @property
    def edges(self):
        """Calculate edges vectors."""
        return get_edges(self.face)

    @property
    def edge_lengths(self):
        """Calculate the length of the edges vectors."""
        return get_edge_lengths(self.edges)

    @property
    def angles(self):
        """Calculate the angles (in radians) between sides sharing a vertex.

        Note that for edges with a length that is almost 0 the angles doesn't say so much.
        """
        return get_edge_angles(self.face)

    def set_related(self, relation: Callable[[int, int, int], bool]) -> None:
        """Set a function to check whether another tetartoid is related by using a, b, c.

        relation: a function accepting values for a, b, and c. The function should return True if it
            is related.
        """
        self.related = relation

    def diff(self, t_cmp: "Tetartoid", simple: bool = True) -> list[float]:
        """Compare self with another tetartoid.

        t_cmp: the tetartoid obj to compare with.
        simple: if True then only a, b, and c are compared, otherwise, edge lengths and angles will
        be added as well.

        return: a list of floats with absolute differences.
        """
        if simple:
            return np.abs(np.array(self.abc) - t_cmp.abc)
        return np.abs(
            # Optionally one could check the angles as well, but
            # note that for edge lengths close to 0 the angles don't matter
            # if edge 'i' is close to 0, i.e from vertex 'i' to 'i+1'
            # then you should ignore the angles 'i' and 'i+1'
            np.array(
                list(self.abc) + self.edge_lengths
            ) - np.array(
                list(t_cmp.abc) + t_cmp.edge_lengths
            )
        )

    def __eq__(self, t_cmp: "Tetartoid") -> bool:
        """Check whether two tetartoids are more or less the same.

        One can use self.atol and self.rtol, see numpy.isclose
        """
        return np.allclose(self.abc, t_cmp.abc, self.rtol, self.atol)

    def unify(self):
        """Scale a, b, c, so that a == 1 to get unique form

        This makes it possible to compare."""
        if not np.isclose(self.a, 0):
            self.rescale(1 / self.a)
        elif not np.isclose(self.b, 0):
            self.rescale(1 / self.b)
        elif not np.isclose(self.c, 0):
            self.rescale(1 / self.c)

    def rescale(self, factor):
        """Scale a, b and c."""
        self.a = factor * self.a
        self.b = factor * self.b
        self.c = factor * self.c

    def save_json(self, filename=None):
        vertices = self.face.tolist()
        self.shape_props["data"]["name"] = self.name
        if not filename:
            filename = self.name + ".json"
        self.shape_props["data"]["base"]["vs"] = vertices
        with open(filename, "w") as fd:
            json.dump(self.shape_props, fd)
        LOGGER.info("saved %s", filename)

    def categorize(self):
        """Return a category for the tetartoid

        The category is related to whether the tetartoid has equal edge lengths or whether angles
        between the sides at the vertices are the same.

        The face in a tetartoid there has five edges: 0, 1, 2, 3, 4.
        Edge 1 connects to an edge 2 of another face in the tetartoid and similarly edge 3 connects
        to an edge 4 of another face. This means that edge 1 and 2 always have the same edge length
        and edge 3 and 4 as well. From a perspective of the edges there are three categories:
            EdgeCat.GENERAL: edge 0 has its own edge length, not equal to the others.
            EdgeCat.E0_EQ_E1_2: edge 0 has the same length as edges 1 and 2
            EdgeCat.E0_EQ_E3_4: edge 0 has the same length as edges 3 and 4
            EdgeCat.ALL_EQ: edge 0 has the same length as edges 1, 2, 3 and 4

        The category is a tuple consisting of:
            - first the edge category
            - then the angle category

        """
        edge_lengths = self.edge_lengths
        if np.isclose(edge_lengths[0], edge_lengths[1]):
            if np.isclose(edge_lengths[0], edge_lengths[3]):
                edge_cat = EdgeCat.ALL_EQ
            else:
                edge_cat = EdgeCat.E0_EQ_E1_2
        elif np.isclose(edge_lengths[0], edge_lengths[3]):
            edge_cat = EdgeCat.E0_EQ_E3_4
        elif np.isclose(edge_lengths[1], edge_lengths[3]):
            edge_cat = EdgeCat.E1_2_EQ_E3_4
        else:
            edge_cat = EdgeCat.GENERAL

        angles = self.angles
        eq_pairs = []
        for i, angle in enumerate(angles):
            for j in range(i + 1, len(angles)):
                if np.isclose(angle, angles[j]):
                    eq_pairs.append((i, j))
                    break
        angle_spec = ()
        match len(eq_pairs):
            case 0:
                angle_cat = AngleCat.GENERAL
            case 1:
                angle_cat = AngleCat.ONE_EQ_PAIR
                angle_spec = tuple(eq_pairs)
            case 2:
                i0, j0 = eq_pairs[0]
                if i0 in eq_pairs[1] or j0 in eq_pairs[1]:
                    LOGGER.warning("Edge lengths: %s", edge_lengths)
                    LOGGER.warning("Angles (rad): %s", angles)
                    raise ValueError("Undefined angle category triplets")
                angle_cat = AngleCat.TWO_EQ_PAIRS
                angle_spec = tuple(eq_pairs)
            case 3:
                # There must be a triplet at least:
                combos = []
                # Note that for (i, j) 'i' won't appear in the following tuples because of the break
                # statement above.
                for pair in eq_pairs:
                    # Add to combos
                    for combo in combos:
                        if pair[0] in combo:
                            combo.append(pair[1])
                            break
                    else:
                        # add as new combo
                        combos.append(list(pair))
                if len(combos) == 2:
                    # should be one pair and a triplet:
                    assert len(combos[0]) == 3 or len(combos[1]) == 3, \
                        "Expected one triplet"
                    assert len(combos[0]) == 2 or len(combos[1]) == 2, \
                        "Expected one pair"
                    angle_cat = AngleCat.EQ_PAIR_AND_TRIPLE
                elif len(combos) == 1:
                    angle_cat = AngleCat.EQ_QUARTET
                else:
                    assert False, "Programming error handling three eq_pairs (quartet?)"
                angle_spec = tuple(combos)
            case 4:
                angle_cat = AngleCat.ALL_EQ
                angle_spec = tuple(eq_pairs)
            case _:
                LOGGER.warning("Edge lengths: %s", edge_lengths)
                LOGGER.warning("Angles (rad): %s", angles)
                self.save_json(filename="unknown.json")
                raise ValueError("Undefined angle category")

        return edge_cat, angle_cat, angle_spec

    @property
    def category(self):
        return self.categorize()

    def log_properties(self):
        cat = self.category
        edge_lengths = self.edge_lengths
        angles = np.rad2deg(self.angles)
        title = f"Tetartoid '{self.name}'"
        line = "-" * len(title)
        LOGGER.info(line)
        LOGGER.info(title)
        LOGGER.info(line)
        if np.isnan(np.sum(angles)):
            LOGGER.warning(">>> Improper tetartoid <<<")
        LOGGER.info("A, B, C = %0.15f, %0.15f, %0.15f", self.a, self.b, self.c)
        LOGGER.info("Edge category: %s", cat[0].name)
        LOGGER.info("Angle category: %s", cat[1].name)
        match cat[0]:
            case EdgeCat.E0_EQ_E1_2:
                LOGGER.info("There are two different edge lengths")
                LOGGER.info("Edge lengths 0, 1, 2 are %0.10f", edge_lengths[0])
                LOGGER.info("Edge lengths 3, 4 are %0.10f", edge_lengths[3])
                LOGGER.info(
                    "difference for 0 is %0.1e",
                    np.abs(edge_lengths[0] - edge_lengths[1]),
                )
            case EdgeCat.E0_EQ_E3_4:
                LOGGER.info("There are two different edge lengths")
                LOGGER.info("Edge lengths 1, 2 are %0.10f", edge_lengths[1])
                LOGGER.info("Edge lengths 0, 3, 4 are %0.10f", edge_lengths[3])
                LOGGER.info(
                    "difference for 0 is %0.1e",
                    np.abs(edge_lengths[0] - edge_lengths[3]),
                )
            case EdgeCat.E1_2_EQ_E3_4:
                LOGGER.info("There are two different edge lengths")
                LOGGER.info("Edge length 0 is %0.10f", edge_lengths[0])
                LOGGER.info("Edge lengths 1, 2, 3 and 4 are %0.10f", edge_lengths[1])
                LOGGER.info(
                    "The differences between 1 and 2, and 3 and 4 are %0.1e",
                    np.abs(edge_lengths[1] - edge_lengths[3]),
                )
            case EdgeCat.ALL_EQ:
                LOGGER.info("All edge lengths are %0.10f", edge_lengths[0])
                LOGGER.info(
                    "differences for 0 are %0.1e and %0.1e",
                    np.abs(edge_lengths[0] - edge_lengths[1]),
                    np.abs(edge_lengths[0] - edge_lengths[3]),
                )
            case EdgeCat.GENERAL:
                LOGGER.info("There are three different edge lengths")
                LOGGER.info("Edge length 0 is %0.10f", edge_lengths[0])
                LOGGER.info("Edge lengths 1, 2 are %0.10f", edge_lengths[1])
                LOGGER.info("Edge lengths 3, 4 are %0.10f", edge_lengths[3])
                LOGGER.info(
                    "differences for 0 are %0.1e and %0.1e",
                    np.abs(edge_lengths[0] - edge_lengths[1]),
                    np.abs(edge_lengths[0] - edge_lengths[3]),
                )
            case _:
                raise ValueError("Unhandled edge category")

        if np.any(np.isclose(edge_lengths, 0)):
            LOGGER.info("Note that some edge lengths are (almost) 0 and angle might not make sense")
        for equals in cat[2]:
            if len(equals) == 2:
                LOGGER.info(
                    "Angle no. %i close to %i with difference of %0.1e°",
                    equals[0],
                    equals[1],
                    np.abs(angles[equals[0]] - angles[equals[1]]),
                )
            elif len(equals) == 3:
                LOGGER.info(
                    "Angle no. %i, %i and %i close to each other and difference with first of "
                    "%0.1e° and %0.1e°",
                    equals[0],
                    equals[1],
                    equals[2],
                    np.abs(angles[equals[0]] - angles[equals[1]]),
                    np.abs(angles[equals[0]] - angles[equals[2]]),
                )
            elif len(equals) == 4:
                LOGGER.info(
                    "Angle no. %i, %i, %i and %i close to each other and difference with first of "
                    "%0.1e°, %0.1e° and %0.1e°",
                    equals[0],
                    equals[1],
                    equals[2],
                    equals[3],
                    np.abs(angles[equals[0]] - angles[equals[1]]),
                    np.abs(angles[equals[0]] - angles[equals[2]]),
                    np.abs(angles[equals[0]] - angles[equals[3]]),
                )
            else:
                #breakpoint()
                assert False, "Programming error handling three eq_pairs (quartet?)"
        for i, angle in enumerate(angles):
            LOGGER.info("Angle at vertex no. %d is %0.10f°", i, angle)


class TetartoidEqEdgeLengths(Tetartoid):

    valid_eq_edge_len = (0, 1, 2, 3)
    valid_eq_angle = (0, 1, 2, 3, 4)

    def __init__(self, init_abc, optimize_for, name=""):
        """
        Initialise object

        init_abc: a tuple with the initial values for 'a', 'b' and 'c', see
            https://en.wikipedia.org/wiki/Dodecahedron#Tetartoid
        optimize_for: a dictionary configuring the optimizer. It can have the following fields:
            'optimizer': the name of the optimizer, see scipy minimize
            'opt_i': specify which indices in init_abc should be optimized
            'eq_edge_len': This parameter will instruct the optimizer to strive for edges with
                equal lengths. For a Tetartoid holds that the edges between vertex 1 and 2 and 2 and
                3 always have the same length and the edges between vertices with indices 3 and 4
                and 4 and 0 always have the same length because the different faces share and edge.
                For the edges between vertex 0 and 1 holds that these share an edge at a 2-fold
                axis. You can specify an array of values. The value affects the length of this edge
                in relation to the other edges:
                0: the length of this edge will not be optimized.
                1: the optimizer will try to parameters for a, b, c so that the edge length becomes
                   equal to the edges between vertex 1 and 2.
                2: the optimizer will try to parameters for a, b, c so that the edge length becomes
                   equal to the edges between vertex 3 and 4.
                Later a the following special value was added:
                3: Try make edges 1 and 2 having the same length as edges 3 and 4.
            'eq_angle': This parameter instructs the optimizer to strive for a equal angles between
                the sides at certain vertices. It is a list of two tuple with two vertex indices
                specifying which angles should be equal.
        name: will be used as Shape name
        """
        if "eq_edge_len" not in optimize_for:
            optimize_for["eq_edge_len"] = []
        if "eq_angle" not in optimize_for:
            optimize_for["eq_angle"] = []
        if "opt_i" not in optimize_for:
            optimize_for["opt_i"] = (1, 2)
        if "method" not in optimize_for:
            optimize_for["method"] = "Powell"
        for value in optimize_for["eq_edge_len"]:
            if value not in self.valid_eq_edge_len:
                raise ValueError(
                    f"Expected one of {self.valid_eq_edge_len} for 'eq_edge_len' parameter, got "
                    f"{value}"
                )
        for i0, i1 in optimize_for["eq_angle"]:
            if i0 not in self.valid_eq_angle or i1 not in self.valid_eq_angle:
                raise ValueError(
                    f"Expected one of {self.valid_eq_angle} for eq_angle tuple, got {i0}, {i1}"
                )
        if len(init_abc) != 3:
            raise ValueError(
                f"Expected three values for 'init_abc' got {init_abc}"
            )
        self._try_abc = list(init_abc)
        self.optimize = optimize_for
        for i in optimize_for["opt_i"]:
            self._try_abc[i] = init_abc[i]
        self.calc_x()
        super().__init__(*self._try_abc, name)

    @property
    def _n(self):
        return tetar_n(*self._try_abc)

    @property
    def _d1(self):
        return tetar_d1(*self._try_abc)

    @property
    def _d2(self):
        return tetar_d2(*self._try_abc)

    @property
    def _face(self):
        e1 = self._n / self._d1
        e2 = self._n / self._d2
        return tetar_face(self._try_abc[0], self._try_abc[1], self._try_abc[2], e1, e2)

    def value_to_minimize(self, values):
        for i, opt_i in enumerate(self.optimize["opt_i"]):
            self._try_abc[opt_i] = values[i]

        face = self._face
        delta_edge_len = 0
        for value in self.optimize["eq_edge_len"]:
            edges = get_edge_lengths(get_edges(face))
            i_j = [None, (0, 1), (0, 3), (1, 3)][value]
            if i_j:
                i, j = i_j
                delta_edge_len += np.abs(edges[i] - edges[j])
        # If |e0| == |e1| (== |e2|)
        # Then preferably the angles (e0, e1) == (e1, e2)
        # Similarly:
        # If |e0| == |e3| (== |e4|)
        # Then preferably the angles (e3, e4) == (e4, e0)
        angle_diff = 0
        for i0, i1 in self.optimize["eq_angle"]:
            angles = get_edge_angles(face)
            angle_diff += np.abs(angles[i0] - angles[i1])
        # use a factor to increase the importance of the angle. TODO: use a parameter
        angle_factor = 180 / np.pi
        return delta_edge_len + angle_factor * angle_diff

    def calc_x(self):
        result = minimize(self.value_to_minimize, self._try_abc, method=self.optimize["method"])
        print(result)
        if not result.success:
            LOGGER.warning("Not really converged, error %e", result.fun)
        else:
            LOGGER.info("Converged, minimum delta %e", result.fun)
        # Check whether the result.fun is close to 0:
        # We can have reached a minimum, but we required it to be close to 0
        # FIXME
        success = result.success or result.fun < 1e-8
        assert success, "Couldn't find 'x' try other input values"
        return result.x


def generate_uniform(outdir: Path):
    """Generate tetartoids with faces or sides that give rise to uniform polyhedra.

    Included are polyhedra where a side, e.g. a square, consists of several other faces, e.g.
    two rectangles.

    outdir: path to directotry where to save the JSON files.
    """
    τ = (np.sqrt(5) + 1) / 2
    τ1 = τ + 1
    tetartoids = {
        # "cube": (0, 1, 1),
        # Use the following instead to get the same file as the cube below
        "cube": (0, τ1, τ1),
        "cubic": (1, 1.5, 1.5),
        "tetrahedron": (1, 1, 3),
        "regular_dodecahedron": (0, τ1, τ1**2),
        "great_stellated_dodecahedron": (0, τ1, 1),
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


def generate_pyritohedra(outdir: Path):
    """Generate tetartoids with faces that have bilateral symmetry.

    These all become pyritohedra, which have more symmetry: S4xI.

    Though I am not sure whether there are more of these.
    I.e. are there more tetartoids with faces with bilateral symmetry that aren't pyritohedra?
    Or are there more tetartoids that are pyritohedra?

    outdir: path to directotry where to save the JSON files.
    """
    τ = (np.sqrt(5) + 1) / 2
    τ1 = τ + 1
    δ = [1e-3, 1e-2]
    big = 4000
    tetartoids = {
        "rhombic_dodecahedron": (0, τ1, big),
        "pyritohedron_slim_pentagons": (0, τ1, τ1**2 + 5),
        "regular_dodecahedron": (0, τ1, τ1**2),
        "pyritohedron_wide_pentagons": (0, τ1, 2 * τ + 1),
        "cube": (0, τ1, τ1),
        "pyritohedron_concave_obtuse": (0, τ1, τ + 1 / 2),
        "endododecahedron": (0, τ1, τ),
        "pyritohedron_concave_sharp": (0, τ1, τ - 1 / 5),
        "d0_by_abc_021_00x": (0, τ1, τ1 / 2 + 2 * δ[0]),
        "pyritohedron_pentagrams_short_single_top": (0, τ1, 10 / 9),
        "great_stellated_dodecahedron": (0, τ1, 1),
        "pyritohedron_pentagrams_long_single_top": (0, τ1, τ - 1),
        "n0_by_ac_00_0x": (0, τ1, δ[1]),  # three crossing lines
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


def generate_at_singularity_case_n_1(outdir: Path):
    """Generate the case where n = 0 by a=b=c."""
    # TODO: remove not used ξ and δ values
    ξ = [1e-4, 1e-3, 5e-3, 1e-5, 1e-2]
    δ = [2e-2, 4e-2, 1e-1, 2e-1, 3e-1]
    tetartoids = {
        "n0_by_abc_111_n00": (1 - δ[2], 1, 1),
        "n0_by_abc_111_x00": (1 - ξ[0], 1, 1),
        "n0_by_abc_111_p00": (1 + δ[2], 1, 1),
        "n0_by_abc_111_0n0": (1, 1 - δ[2], 1),
        "n0_by_abc_111_0x0": (1, 1 - ξ[3], 1),
        "n0_by_abc_111_0p0": (1, 1 + δ[2], 1),
        "n0_by_abc_111_00n": (1, 1, 1 - δ[2]),
        "n0_by_abc_111_00x": (1, 1, 1 - ξ[0]),
        "n0_by_abc_111_00p": (1, 1, 1 + δ[2]),
        "n0_by_abc_111_np0": (1 - δ[2], 1 + δ[2], 1),
        "n0_by_abc_111_-+0": (1 - ξ[0], 1 + ξ[0], 1),
        "n0_by_abc_111_pn0": (1 + δ[2], 1 - δ[2], 1),
        "n0_by_abc_111_n0p": (1 - δ[2], 1, 1 + δ[2]),
        "n0_by_abc_111_-0+": (1 - ξ[0], 1, 1 + ξ[0]),
        "n0_by_abc_111_p0n": (1 + δ[2], 1, 1 - δ[2]),
        "n0_by_abc_111_0np": (1, 1 - δ[2], 1 + δ[2]),
        "n0_by_abc_111_0-+": (1, 1 - ξ[0], 1 + ξ[0]),
        "n0_by_abc_111_0pn": (1, 1 + δ[2], 1 - δ[2]),
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


def generate_at_singularity_case_n_2(outdir: Path):
    """Generate the case where n = 0 by bc/aa=1."""
    # TODO: remove not used ξ and δ values
    ξ = [1e-4, 1e-3, 5e-3, 1e-5, 1e-2]
    δ = [2e-2, 4e-2, 1e-1, 2e-1, 3e-1]
    tetartoids = {
        "n0_by_abc_1xx-1_0n0": (1, 5 / 4 - δ[2], 4 / 5),
        "n0_by_abc_1xx-1_0x0": (1, 5 / 4 - ξ[0], 4 / 5),
        "n0_by_abc_1xx-1_0p0": (1, 5 / 4 + δ[2], 4 / 5),
        "n0_by_abc_1x-1x_0n0": (1, 4 / 5 - δ[0], 5 / 4),
        "n0_by_abc_1x-1x_0x0": (1, 4 / 5 - ξ[0], 5 / 4),
        "n0_by_abc_1x-1x_0p0": (1, 4 / 5 + δ[0], 5 / 4),
        "n0_by_abc_1xx-1_00p": (1, 4 / 5, 5 / 4 + 0.16),
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


def generate_at_singularity_case_n_3(outdir: Path):
    """Generate the case where n = 0 by c=0."""
    # TODO: remove not used ξ and δ values
    ξ = [1e-4, 1e-3, 5e-3, 1e-5, 1e-2]
    δ = [2e-2, 4e-2, 1e-1, 2e-1, 3e-1]
    tetartoids = {
        "n0_by_ac_00_0n": (0, 2, 0 - δ[3]),
        "n0_by_ac_00_0x": (0, 2, 0 - ξ[4]),
        "n0_by_ac_00_0p": (0, 2, 0 + δ[3]),
        "n0_by_ac_10_0n": (1, 2, 0 - δ[3]),
        "n0_by_ac_10_0x": (1, 2, 0 - ξ[4]),
        "n0_by_ac_10_0p": (1, 2, 0 + δ[3]),
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


def generate_at_singularity_case_d_1(outdir: Path):
    """Generate the case where n = 0 by a=b=0b=c."""
    # TODO: remove not used ξ and δ values
    ξ = [1e-4, 1e-3, 5e-3, 1e-5, 1e-2]
    δ = [2e-2, 4e-2, 1e-1, 2e-1, 3e-1]
    tetartoids = {
        "d0_by_ab_00_n0": (-δ[0], 0, 1),
        "d0_by_ab_00_x0": (ξ[1], 0, 1),
        "d0_by_ab_00_p0": (δ[0], 0, 1),
        "d0_by_ab_00_0n": (0, -δ[1], 1),
        "d0_by_ab_00_0x": (0, ξ[1], 1),
        "d0_by_ab_00_0p": (0, δ[1], 1),
        "d0_by_ab_00_nn": (-δ[1], -δ[1], 1),
        "d0_by_ab_00_--": (-ξ[1], -ξ[1], 1),
        "d0_by_ab_00_pp": (δ[1], δ[1], 1),
        "d0_by_ab_00_np": (-δ[1], δ[1], 1),
        "d0_by_ab_00_-+": (-ξ[1], ξ[1], 1),
        "d0_by_ab_00_pn": (δ[1], -δ[1], 1),
        # "d0_by_ab_00_+-": (ξ[1], -ξ[1], 1), same as -+
        # "d0_by_ab_00_++": (ξ[1], ξ[1], 1), same as ++
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


def generate_at_singularity_case_d_2(outdir: Path):
    """Generate the case where n = 0 by b=0 and c=-a."""
    # TODO: remove not used ξ and δ values
    ξ = [1e-4, 1e-3, 5e-3, 1e-5, 1e-2]
    δ = [2e-2, 4e-2, 1e-1, 2e-1, 3e-1]
    tetartoids = {
        "d0_by_abc_10-1_n00": (1 - δ[2], 0, -1),
        "d0_by_abc_10-1_x00": (1 - ξ[2], 0, -1),
        "d0_by_abc_10-1_p00": (1 + δ[2], 0, -1),
        "d0_by_abc_10-1_0n0": (1, 0 - δ[2], -1),
        "d0_by_abc_10-1_0x0": (1, 0 + ξ[2], -1),
        "d0_by_abc_10-1_0p0": (1, 0 + δ[2], -1),
        "d0_by_abc_10-1_00n": (1, 0, -1 - δ[2]),
        "d0_by_abc_10-1_00x": (1, 0, -1 + ξ[2]),
        "d0_by_abc_10-1_00p": (1, 0, -1 + δ[2]),
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


def generate_at_singularity_case_d_3(outdir: Path):
    """Generate the case where n = 0 by b=0 and c=a."""
    # TODO: remove not used ξ and δ values
    ξ = [1e-4, 1e-3, 5e-3, 1e-5, 1e-2]
    δ = [2e-2, 4e-2, 1e-1, 2e-1, 3e-1]
    tetartoids = {
        "d0_by_abc_101_n00": (1 - δ[2], 0, 1),
        "d0_by_abc_101_x00": (1 - ξ[2], 0, 1),
        "d0_by_abc_101_p00": (1 + δ[2], 0, 1),
        "d0_by_abc_101_0n0": (1, 0 - δ[2], 1),
        "d0_by_abc_101_0x0": (1, 0 - ξ[2], 1),
        "d0_by_abc_101_0p0": (1, 0 - δ[2], 1),
        "d0_by_abc_101_00n": (1, 0, 1 - δ[2]),
        "d0_by_abc_101_00x": (1, 0, 1 + ξ[2]),
        "d0_by_abc_101_00p": (1, 0, 1 + δ[2]),
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


def generate_at_singularity_case_d_4(outdir: Path):
    """Generate the case where n = 0 by b≠0 and c=f1×b."""
    # TODO: remove not used ξ and δ values
    ξ = [1e-4, 1e-3, 5e-3, 1e-5, 1e-2]
    δ = [2e-2, 4e-2, 1e-1, 2e-1, 3e-1]
    tetartoids = {
        # choose a = 1 and b = 4/5 to be similar to bc/aa=1
        # Now here fa = 4/5 and f1 = 7/4, and c = 7/5
        # Now scale to b = 1
        "d1_0_by_abc_fabf1_00n": (5 / 4, 1, 7 / 4 - δ[1]),
        "d1_0_by_abc_fabf1_00x": (5 / 4, 1, 7 / 4 - ξ[2]),
        "d1_0_by_abc_fabf1_00p": (5 / 4, 1, 7 / 4 + δ[1]),
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


def generate_at_singularity_case_d_5(outdir: Path):
    """Generate the case where n = 0 by b≠0 and c=f2×b."""
    # TODO: remove not used ξ and δ values
    ξ = [1e-4, 1e-3, 5e-3, 1e-5, 1e-2]
    δ = [2e-2, 4e-2, 1e-1, 2e-1, 3e-1]
    tetartoids = {
        # Take b = 1, fa = 0.5, then f2 = 7/10
        "d2_0_by_abc_fabf2_00n": (1 / 2, 1, 7 / 10 - δ[0]),
        "d2_0_by_abc_fabf2_00x": (1 / 2, 1, 7 / 10 + ξ[1]),
        "d2_0_by_abc_fabf2_00p": (1 / 2, 1, 7 / 10 + δ[0]),
    }
    for name, abc in tetartoids.items():
        Tetartoid(*abc).save_json(outdir / (name + ".json"))


NAMED_SET_MAP = {
    "pyritohedra": generate_pyritohedra,
    "uniform polyhedra": generate_uniform,
    # singularities when c = 0
    # ------------------------
    "a=b=c": generate_at_singularity_case_n_1,
    "bc/aa=1": generate_at_singularity_case_n_2,
    "c=0": generate_at_singularity_case_n_3,
    # singularities when d = 0
    # ------------------------
    "a=b=0": generate_at_singularity_case_d_1,
    "b=0 and c=-a": generate_at_singularity_case_d_2,
    "b=0 and c=a": generate_at_singularity_case_d_3,
    "b≠0 and c=f1×b": generate_at_singularity_case_d_4,
    "b≠0 and c=f2×b": generate_at_singularity_case_d_5,
}


if __name__ == "__main__":
    groups = list(NAMED_SET_MAP.keys())
    groups.append("pyritohedra")
    groups.append("uniform polyhedra")
    parser = argparse.ArgumentParser(description=DESCR)
    parser.add_argument(
        "named_set",
        help=f"Named set up tetartoids. Should be one of {groups}",
    )
    parser.add_argument(
        "-o", "--outdir",
        default=".",
        help="Specify the output directory of the JSON files. The directory must exist prior to "
        "the call.",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir) if args.outdir else Path(".")

    set_of_tetartoids = {
    }

    def add_to_set(tetartoid):
        tetartoid.unify()
        set_of_tetartoids[tetartoid.abc] = tetartoid

    def find_in_set(tetartoid):
        tetartoid.unify()
        candidates = {}
        for t in set_of_tetartoids.values():
            if tetartoid == t:
                candidates[t.name] = t
            else:
                # if a relation is defined check it:
                if t.related:
                    if t.related(*tetartoid.abc):
                        candidates[t.name] = t
        if len(candidates) == 0:
            return None
        if len(candidates) == 1:
            for t_d in candidates.values():
                return t_d
        diffs = [
            np.sum(tetartoid.diff(t_d, simple=False))
            for t_d in candidates.values()
        ]
        objs = np.array(list(candidates.values()))
        result = objs[np.array(diffs).argmin()]
        return result

    #t = TetartoidOneReq(a=1, b=1.1)
    #t = TetartoidOneReq(b=1, c=1.1)
    #t = TetartoidOneReq(a=0.4, c=2)

    # It seems that if you start with b=1 then you will get a convex polyhedron.
    # That happened at least with Nelder-Mead and Powell and the default method
    # Powel seems to converge better.
    #t = TetartoidEqEdgeLengths(b=1, method="Powell", eq_edge_len=1)
    #t = TetartoidEqEdgeLengths(a=0, method="Nelder-Mead", eq_edge_len=1)
    #t = TetartoidEqEdgeLengths(b=0, method="BFGS", eq_edge_len=1)

    # Hats...
    #t = TetartoidEqEdgeLengths((0.9, 0.9, 2), keep_index=2, method="Powell", eq_edge_len=2)
    # Not really converged: 3 edge lengths
    # t = TetartoidEqEdgeLengths((0.4, 0.8, 2.0), keep_index=0, method="Powell", eq_edge_len=2)
    # Same result:
    # t = TetartoidEqEdgeLengths((0.4, 0.8, 2.0), keep_index=0, method="Powell", eq_edge_len=2)

    # Self intersecting face and butterfly shaped
    # Not really converged: no equal angles
    #    (0.4, 0.8, 2.0), keep_index=2, method="Powell", eq_edge_len=1, eq_angle=1

    # Regular dodecahedron.
    # t = TetartoidEqEdgeLengths((0.1, 1.0, 2.0), keep_index=1, method="Powell", eq_edge_len=2)
    # t = TetartoidEqEdgeLengths((0.4, 0.8, 2.0), keep_index=1, method="Powell", eq_edge_len=2)
    # t = TetartoidEqEdgeLengths((0.4, 0.8, 2.0), keep_index=1, method="Powell", eq_edge_len=1)
    # (0.4, 0.8, 2.0), keep_index=2, method="Powell", eq_edge_len=1, eq_angle=3

    name = "regular_dodecahedron"

    # I tried combinations of eq_angle. They all lead to "almost" regular dodecahedra and they
    # didn't optimize completely. I.e. would they go the whole way, then it would be a regular
    # regular dodecahedron

    tau = (np.sqrt(5) + 1) / 2
    try_methods = ("Powell", "Nelder-Mead", "COBYQA", "BFGS", "SLSQP")
    opt_setup = {
        # ######################################################
        # Regular shaped
        # ######################################################
        # edge: ALL_EQ
        # angle: ALL_EQ
        "regular_dodecahedron": {
            "abc": (0, 1, tau + 1),
        },
        # edge: ALL_EQ
        # angle: ALL_EQ
        "great_stellated_dodecahedron": {
            "abc": (0, tau + 1, 1),
        },
        # edge: ALL_EQ
        # angle: TWO_EQ_PAIRS
        "tetrahedron": {
            "abc": (1, 1, 3),
        },
        # edge: GENERAL
        # angle: ONE_EQ_PAIR
        # each side consists of two right trapezoids
        "cubic": {
            "abc": (1, 1.5, 1.5),
            "related": lambda a, b, c: np.isclose(a, 1) and np.isclose(b, c) and (b > 1 or b < -1),
        },
        # edge: E1_2_EQ_E3_4
        # angle: EQ_QUARTET
        # each side consists of two rectangles
        "cube": {
            "abc": (0, 1, 1),
            "related": lambda a, b, c: np.isclose(a, 1) and np.isclose(b, c) and b > 100,
        },
        # ######################################################
        # On the limit
        # ######################################################
        # edge: GENERAL
        # angle: ONE_EQ_PAIR
        # each side consists of two triangles
        "almost_cube": {
            "abc": (1 - 1e-10, 1, 1),
            # Same result for
            # "abc": (1 + 1e-10, 1, 1),
        },
        # edge: E1_2_EQ_E3_4
        # angle: EQ_PAIR_AND_TRIPLE
        # should be 1, lim x->0: x, x
        "3_crossing_lines": {
            # degenerate
            "abc": (1, 1e-14, 1e-14),
            # TODO: b can have any value
        },
        # edge: GENERAL
        # angle: TWO_EQ_PAIRS
        # each side is covered 3 times
        "almost_tetrahedron_0": {
            # abc = 1, 1, lim x↑1 x
            "abc": (1, 1, 1 - 1e-10)
            # same result for
            # abc = 1, 1, lim x↓1 x
            # "abc": (1, 1, 1 + 1e-10)
        },
        "almost_tetrahedron_0_alt": {
            # abc = 1, 1, lim x↓1 x
            "abc": (1, 1, 1 + 1e-10)
        },
        # edge: E0_EQ_E3_4
        # angle: ONE_EQ_PAIR
        # Each side consists of triangles meeting in the side centre
        "almost_tetrahedron_1": {
            # abc = 1, 1, lim x↑1 x
            "abc": (1, 1 - 1e-10, 1)
            # same result for
            # abc = 1, lim x↓1 x, 1
            # "abc": (1, 1 + 1e-10, 1)
        },
        "almost_tetrahedron_1_alt": {
            # abc = 1, lim x↓1 x, 1
            "abc": (1, 1 + 1e-10, 1)
        },
        # edge: E0_EQ_E3_4 (δ = 2.7e-11)
        # angle: ONE_EQ_PAIR (δ = 4.5e-10)
        # irregular tetraaugmented tetrahedron
        "four_tetras": {
            "abc": (1, 1 - 1.5e-11, 1 - 1e-11),
            "related": lambda a, b, c: np.isclose(a, 1) and \
                np.isclose((1 - b) / (1 - c), 1.5) and \
                np.isclose(b, 1),
        },
        # edge: E0_EQ_E1_2 (δ = 2.7e-11)
        # angle: ONE_EQ_PAIR (δ = 4.5e-10)
        # Four tetrahedra on one tetrahedron
        "four_tetras_alt": {
            "abc": (1, -1 - 1.5e-11, -1 - 1e-11),
            "related": lambda a, b, c: np.isclose(a, 1) and \
                np.isclose((1 + b) / (1 + c), 1.5) and \
                np.isclose(b, -1),
        },

        # ######################################################
        # Other
        # ######################################################
        # edge: ALL_EQ
        # angle: EQ_PAIR_AND_TRIPLE
        "extended_regular_dodecahedron": {
            "abc": (0, 1, -tau),
        },
        # edge: GENERAL
        # angle: TWO_EQ_PAIRS
        # A whole series for which a=1, 0 < b=c < 1
        "extended_cube": {
            "abc": (1, 0.8, 0.8),
            "related": lambda a, b, c: np.isclose(a, 1) and np.isclose(b, c) and -1 < b < 1,
        },
        "pentaspikes": {
            "start_with": (0.4, 0.8, 2.0),
            "optimize_for": {
                "opt_i": (1, 2),
                "eq_edge_len": [1],
            },
        },
        "almost_bilateral": {
            "start_with": (0.4, 0.8, 2.0),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [1],
            },
        },
        "v_shape_face_self_intersect_butterfly": {
            # Local minimum (delta = 1.9)
            # With method Nelder-Mead / SLSQP a regular dodecahedron is obtained
            # method: "COBYQA" see below
            "start_with": (0.4, 0.8, 2.0),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [1],
                "eq_angle": [(0, 1)],
            },
        },
        # edge: E0_EQ_E1_2
        # angle: TWO_EQ_PAIRS
        # δ = 1.5e-4
        "v_shape_face_butterfly_0": {
            # "abc": (1, 0.210138278109585, 2.484461652493819),
            "start_with": (1., 0.2, 2.5),
            "optimize_for": {
                "opt_i": (0, 1),
                "method": try_methods[1],
                "eq_edge_len": [1],
                "eq_angle": [(0, 4), ],
            },
        },
        # edge: E0_EQ_E3_4
        # angle: TWO_EQ_PAIRS
        # δ = 1.2e-7
        "v_shape_face_butterfly_1": {
            # "abc": (1, 0.484454027157523, 1.274316688874315),
            "start_with": (1., 0.2, 2.5),
            "optimize_for": {
                "opt_i": (0, 1),
                "method": try_methods[1],
                "eq_edge_len": [2],
                "eq_angle": [(0, 4), ],
            },
        },
        # edge: GENERAL
        # angle: TWO_EQ_PAIRS
        # δ = 1.5e-10
        "v_shape_face_butterfly_alt": {
            # "abc": (1, 0.204142371201545, 2.551342236267743),
            # Note
            # that even though we are optimizing for one edge length and one angle, the result has
            # two equal angles and no edges with the same length
            # Changing eq_edge_len to 2 generates a slightly different model, but still having the
            # same category:
            # "abc": (1, 0.204125029072595, 2.551541651277570)
            # and δ = 2.6e-12
            "start_with": (1., 0.2, 2.5),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [1],
                "eq_angle": [(0, 4), ],
            },
        },
        "spiky_butterfly": {
            # local minimum: (minimum = 8°)
            "start_with": (0.2, 0.4, 1),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_angle": [(1, 2), (0, 4)],
            },
        },
        "4_point_star_self_inters": {
            "start_with": (1, 2.3, 1.3),
            "optimize_for": {
                "opt_i": (1, 2),
                #"method": try_methods[2],
                "eq_edge_len": [1],
                "eq_angle": [(3, 4), ],
            },
        },
        # The following ones come in many variations, since only one requirement is met:
        # ------------------------------------------------------------
        # edge: E0_EQ_E1_E2
        # angle: GENERAL
        # δ = 0.0
        "classic_tetartoid": {
            # "abc": (1, -3.671535138208082, -9.339026119964645),
            "start_with": (0.4, 0.8, 2),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [1],
            },
        },
        # edge: E0_EQ_E3_E4
        # angle: GENERAL
        # δ = 3.5e-15
        "classic_tetartoid_1": {
            # "abc": (1, 3.615386818341456, 9.196464579601590),
            "start_with": (0.4, 0.8, 2),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [1],
            },
        },
        "self_inters_face_pyramids_0": {
            "start_with": (1, 0.53, 0.73),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_angle": [(3, 4), ],
            },
        },
        "like_52nd_ico_stellation_0": {
            "start_with": (1, 3.0, 2.7),
            "optimize_for": {
                "opt_i": (1, 2),
                "eq_edge_len": [1],
                "eq_angle": [(3, 4), ],
            },
        },
        # edge: E0_EQ_E3_4
        # angle: GENERAL
        "v_shape_face_pyramids_1": {
            # 1, -0.079390874047626, -2.551894699132569
            "start_with": (1.1, 0.0, 2.0),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [2],
            },
        },
        # edge: E0_EQ_E1_2
        # angle: GENERAL
        "v_shape_face_pyramids_inv": {
            # 1, , 0.004761203017896, 2.563964321119920
            "start_with": (1.1, 0.0, 2.0),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [1],
            },
        },
        # edge: GENERAL
        # angle: ONE_EQ_PAIR
        "v_shape_face_pyramids_small": {
            "start_with": (1.6, -0.1, 1.5),
            "optimize_for": {
                "opt_i": (0, 1),
                "method": try_methods[4],
                "eq_angle": [(3, 4), (0, 1)],
            },
        },
        "twist": {
            # local minimum: delta 1.5
            "start_with": (0.2, 0.4, 1),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_angle": [(3, 4), (0, 2)],
            },
        },
        "v_shape_face_butterfly_0": {
            "start_with": (0.4, 0.8, 2.0),
            "optimize_for": {
                "opt_i": (0, 1),
                "method": "COBYQA",
                "eq_edge_len": [1],
                "eq_angle": [(0, 1)],
            },
        },
    }

    def find_tetartoid(name):
        """Find a pre-defined tetartoid in opt_setup

        name: the name used in set_of_tetartoids

        return: the Tetartoid object
        """
        setup = opt_setup[name]
        if "abc" in setup:
            tetartoid = Tetartoid(*setup["abc"], name=name)
            if "related" in setup:
                tetartoid.set_related(setup["related"])
            return tetartoid
        else:
            return TetartoidEqEdgeLengths(setup["start_with"], setup["optimize_for"], name=name)

    for name in opt_setup.keys():
        t = find_tetartoid(name)
        # always add, don't check:
        # if find_in_set(t) is None:
        # since some abc are very close but with very different results
        add_to_set(t)
        t.save_json()

    for t in set_of_tetartoids.values():
        t.log_properties()

    # TODO: You can still require that lengths 1,2 == 3, 4
    # This is what you want for
    # start_with = (1., x, x)

    # Try:
    #tau = (np.sqrt(5) + 1) / 2
    #start_with = (0., 1.0, -tau)
    start_with = (1, 2.3, 1.3)
    start_with = (1, 3.0, 2.7)
    start_with = (1.1, 0., 2.0)
    start_with = (1.0, 1., .9)
    optimize_for = {
        "opt_i": (1, 2),
        "method": try_methods[1],
        #"eq_edge_len": [3],
        "eq_angle": [(3, 4), (0, 1)],
    }
    # Set if you want to test a, b, c directly
    abc = ()
    # extended tetrahedron: abc = 1.0, 1-2e-11, 1-1e-1
    # delta = -1e-1
    # delta = 1e-1
    delta = 2e-1
    factor = 1.5  # differences for 0 are 2.8e+00 and 1.2e-04
    abc = 1.0, -1 + factor * delta, -1 + delta

    # Investigate approaching limits (e.g. d1 = 0)
    f_a = 3
    f_c = (f_a**2 - f_a + 1) / (2 - f_a)
    δ = -1e-1
    #abc = f_a + δ, 1, f_c
    #abc = f_a, 1 + δ, f_c
    abc = f_a, 1, f_c + δ

    # try out when the tetartoid is convex / concave
    # I think requirement 1 is related to that.
    # we have
    # 1   3    2.7
    # 1+ξ 1-ξ  1
    # With a and b symmetrical: c < a or c < b?
    # 1, 3, 3 is square with right trapezoids
    # Tried: 1.65, 3, 2.2
    # If you increase 'a' to 1.6 then the point between the tops gets
    # so low that the face self-intersects.
    #      => At least b < c < a or a < c < b though this isn't enough
    # Lowering the value of 'a' will lower the top
    # around abc = 0, 3, 2.2
    # the tops are (sort of) equal
    # around abc = 0.0, 3, 2.2
    # flatter: abc = 0.0, 3, 2.5
    # very spiky: abc = 0.0, 3, 1.6
    # At 0, 3, 1.5 d1 = 0 -> singular point
    # leads to almost Great stellated dodecahedron: find out when that is happening..
    abc = 0, tau + 1, 1e-4

    # d1 = 0
    abc = -2e-2, 0, 1

    if args.named_set:
        NAMED_SET_MAP[args.named_set](outdir)

    #if abc:
    #    t = Tetartoid(*abc)
    #    t.save_json()
    #else:
    #    t = TetartoidEqEdgeLengths(start_with, optimize_for, name="test")
    #LOGGER.info("=================================")
    #t1 = find_in_set(t)

    #if t1 is None:
    #    LOGGER.info("*** New one found!")
    #    t.log_properties()
    #    t.save_json()
    #else:
    #    LOGGER.info("*** Exists as %s", t1.name)
    #    LOGGER.info("===============OLD===============")
    #    t1.log_properties()
    #    LOGGER.info("===============NEW===============")
    #    t.log_properties()
