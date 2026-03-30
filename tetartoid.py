from enum import IntEnum, auto
import json
import logging

import numpy as np
from scipy.optimize import minimize

class EdgeCat(IntEnum):
    GENERAL = auto()  # there are three different edge lengths
    E0_EQ_E1_2 = auto()  # edge 0 has the same length as edge 1 and 2
    E0_EQ_E3_4 = auto()  # edge 0 has the same length as edge 3 and 4
    E1_2_EQ_E3_4 = auto()  # edge 1 and 2 have the same length as edge 3 and 4
    ALL_EQ = auto()  # all edges have equal length

class AngleCat(IntEnum):
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

def tetar_n(a, b, c):
    return a**2 * c - b * c**2

def tetar_d1(a, b, c):
    return a * (a - b + c) + b * (b - 2 * c)

def tetar_d2(a, b, c):
    return a * (a + b - c) + b * (b - 2 * c)

def tetar_face(a, b, c, e1, e2):
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


def get_edge_lengths(edges):
    return np.linalg.norm(edges, axis=1)


def get_edge_angles(face):
    edges = get_edges(face)
    edge_lengths = get_edge_lengths(edges)
    no_of_vs = len(edges)
    return [
        np.arccos(
            np.dot(edges[i], edges[(i + 1) % no_of_vs]) / (
                edge_lengths[i] * edge_lengths[(i + 1) % no_of_vs]
            )
        )
        for i in range(no_of_vs)
    ]


class Tetartoid():
    rtol = 1e-5
    atol = 1e-8

    def __init__(self, a, b, c, name="A4"):
        """
        name: will be used as Shape name
        """
        #assert a <= b <= c
        self.a = a
        self.b = b
        self.c = c
        self.name = name

        self.shape_props = {
            "class": "orbit_shape",
            "data": {
                "cols": [
                    [64, 104, 224],
                    [64, 104, 224],
                    [153, 204, 49],
                    [153, 204, 49],
                    [220, 159, 220],
                    [220, 159, 220],
                    [254, 214, 0],
                    [254, 214, 0],
                    [153, 204, 49],
                    [64, 104, 224],
                    [220, 159, 220],
                ],
                "final_sym": {
                    "class": "A4",
                    "data": {
                        "generator": {
                            "o2axis0": [1, 0, 0],
                            "o2axis1": [0, 1, 0],
                        }
                    },
                },
                "fs": [[0, 1, 2, 3, 4]],
                "name": name if name else "A4",
                "no_of_cols": 12,
                "stab_sym": {"class": "E", "data": {}},
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

    def __eq__(self, t):
        """Check whether two tetartoids are more or less the same.

        One can use self.atol and self.rtol, see numpy.isclose
        """
        return np.allclose(self.abc, t.abc, self.rtol, self.atol)

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
        self.shape_props["data"]["vs"] = vertices
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
        edges = get_edges(self.face)
        edge_lengths = get_edge_lengths(edges)
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

        angles = get_edge_angles(self.face)
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
                t.save_json(filename="unknown.json")
                raise ValueError("Undefined angle category")

        return edge_cat, angle_cat, angle_spec

    @property
    def category(self):
        return self.categorize()

    def log_properties(self):
        cat = self.category
        edges = get_edges(self.face)
        edge_lengths = get_edge_lengths(edges)
        angles = get_edge_angles(self.face)
        title = f"Tetartoid '{self.name}'"
        line = "-" * len(title)
        LOGGER.info(line)
        LOGGER.info(title)
        LOGGER.info(line)
        if np.isnan(np.sum(angles)):
            LOGGER.warning(">>> Inproper tetartoid <<<")
        LOGGER.info("A, B, C = %f, %f, %f", self.a, self.b, self.c)
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
                LOGGER.info("Edge lengths 1, 2 are %0.10f", edge_lengths[0])
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

        for equals in cat[2]:
            if len(equals) == 2:
                LOGGER.info(
                    "Angle no. %i close to %i with difference of %0.1e°",
                    equals[0],
                    equals[1],
                    np.rad2deg(np.abs(angles[equals[0]] - angles[equals[1]])),
                )
            elif len(equals) == 3:
                LOGGER.info(
                    "Angle no. %i, %i and %i close to each other and difference with first of "
                    "%0.1e° and %0.1e°",
                    equals[0],
                    equals[1],
                    equals[2],
                    np.rad2deg(np.abs(angles[equals[0]] - angles[equals[1]])),
                    np.rad2deg(np.abs(angles[equals[0]] - angles[equals[2]])),
                )
            elif len(equals) == 4:
                LOGGER.info(
                    "Angle no. %i, %i, %i and %i close to each other and difference with first of "
                    "%0.1e°, %0.1e° and %0.1e°",
                    equals[0],
                    equals[1],
                    equals[2],
                    equals[3],
                    np.rad2deg(np.abs(angles[equals[0]] - angles[equals[1]])),
                    np.rad2deg(np.abs(angles[equals[0]] - angles[equals[2]])),
                    np.rad2deg(np.abs(angles[equals[0]] - angles[equals[3]])),
                )
            else:
                breakpoint()
                assert False, "Programming error handling three eq_pairs (quartet?)"
        for i, angle in enumerate(angles):
            LOGGER.info("Angle at vertex no. %d is %0.10f°", i, np.rad2deg(angle))


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
            if value < 3:
                # 1. edges v[1] - v[2] and v[2] - v[3] are shared: i.e. they have the same length
                # 2. edges v[3] - v[4] and v[4] - v[0] are also shared.
                # edge length v[0] - v[1] should be equal to either of these
                # They meet each other at a 2-fold axis
                d_edge = edges[1] if value else edges[3]
                delta_edge_len += np.abs(edges[0] - d_edge)
            else:
                delta_edge_len += np.abs(edges[1] - edges[3])
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


if __name__ == "__main__":
    set_of_tetartoids = {
    }

    def add_to_set(tetartoid):
        tetartoid.unify()
        set_of_tetartoids[tetartoid.abc] = tetartoid

    def find_in_set(tetartoid):
        tetartoid.unify()
        result = None
        for t in set_of_tetartoids.values():
            if tetartoid == t:
                result = t
                break
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

    name = "almost_tetrahedron"
    name = "regular_dodecahedron"

    # I tried combinations of eq_angle. They all lead to "almost" regular dodecahedra and they
    # didn't optimize completely. I.e. would they go the whole way, then it would be a regular
    # regular dodecahedron

    tau = (np.sqrt(5) + 1) / 2
    try_methods = ("Powell", "Nelder-Mead", "COBYQA", "BFGS", "SLSQP")
    opt_setup = {
        # edge: ALL_EQ, angle: ALL_EQ (1.5e-4)
        # TODO: Use the real values, don't optimize
        "regular_dodecahedron": {
            "start_with": (0.0, 0.8, 2.0),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_angle": [(0, 3), (1, 2)],
            },
        },
        # edge: ALL_EQ, angle: EQ_PAIR_AND_TRIPLE
        # TODO: don't optimize
        "extended_regular_dodecahedron": {
            "start_with": (1., 0.2, 2.5),
            # Leads to a divide by 0 in method _face for e2:
            # "start_with": (0.0, 1.0, -tau),
            "optimize_for": {
                "opt_i": (0, 2),
                "eq_edge_len": [2],
                "eq_angle": [(0, 4), (1, 3)],
            },
        },
        # TODO: specify
        "cube": {
            "abc": (0, 1, 1),
        },
        # edge: GENERAL, angle: TWO_EQ_PAIRS
        # A whole series for which a=1, 0 < b=c < 1
        # TODO: need a way to check which one we found
        "extended_cube": {
            "abc": (1, 0.8, 0.8),
        },
        "classic_tetartoid": {
            "start_with": (0.4, 0.8, 2),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [2],
            },
        },
        "four_tetras": {
            "start_with": (1.6, -0.1, 1.5),
            "optimize_for": {
                "opt_i": (0, 1),
                "method": try_methods[1],
                "eq_angle": [(3, 4), (0, 1)],
            },
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
        # edge: E0_EQ_E1_2, angle: TWO_EQ_PAIRS (1.5e-4)
        "v_shape_face_butterfly_alt": {
            "start_with": (1., 0.2, 2.5),
            "optimize_for": {
                "opt_i": (0, 1),
                "method": try_methods[1],
                "eq_edge_len": [2],
                "eq_angle": [(0, 4), ],
            },
        },
        # edge: GENERAL, angle: TWO_EQ_PAIRS (1.5e-10)
        "v_shape_face_butterfly": {
            "start_with": (1., 0.2, 2.5),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [2],
                "eq_angle": [(0, 4), ],
            },
        },
        "v_shape_face_pyramids_small": {
            "start_with": (1.6, -0.1, 1.5),
            "optimize_for": {
                "opt_i": (0, 1),
                "method": try_methods[4],
                "eq_angle": [(3, 4), (0, 1)],
            },
        },
        "almost_tetrahedron": {
            # doesn't optimize well for angles (mininum = 7.3e-6)
            "start_with": (-0.67, 0.67, 2),
            "optimize_for": {
                "opt_i": (0, 1),
                "method": "Nelder-Mead",
                "eq_angle": [(1, 2), (0, 3)],
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
        "v_shape_face_pyramids": {
            "start_with": (1.1, 0.0, 2.0),
            "optimize_for": {
                "opt_i": (0, 1),
                "eq_edge_len": [2],
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
        setup = opt_setup[name]
        if "abc" in setup:
            return Tetartoid(*setup["abc"], name=name)
        else:
            return TetartoidEqEdgeLengths(setup["start_with"], setup["optimize_for"], name=name)

    for name in opt_setup.keys():
        t = find_tetartoid(name)
        if find_in_set(t) is None:
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
    start_with = (1., 0.2, 2.5)
    #start_with = (1., 1.5, 1.5)
    optimize_for = {
        "opt_i": (0, 2),
        #"method": try_methods[1],
        "eq_edge_len": [3],
        #"eq_angle": [(0, 4), ],
    }
    t = TetartoidEqEdgeLengths(start_with, optimize_for, name="test")
    LOGGER.info("=================================")
    t1 = find_in_set(t)

    if t1 is None:
        LOGGER.info("*** New one found!")
        t.log_properties()
        t.save_json()
    else:
        LOGGER.info("*** Exists as %s", t1.name)
        LOGGER.info("===============OLD===============")
        t1.log_properties()
        LOGGER.info("===============NEW===============")
        t.log_properties()

    # Check directly
    # Cube:
    #a, b, c = 1e-9, 1, 1
    a, b, c = 1, 1e-14, 1e-14
    t = Tetartoid(a, b, c)
    t.unify()
    t.log_properties()
    print(t.face)
    t.save_json("checking.json")

    # TODO: handle 1, 1, 1 (tetrahedron)
    # a, b, c = 1, 1, 1 - 1e-12
    # t = Tetartoid(a, b, c)
    # t.log_properties()
    # t.save_json()

