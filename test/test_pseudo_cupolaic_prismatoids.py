"""Test pseudo cupolaic prismatoid generartor

Generate models and compare with existing one
"""

import os
import pathlib
import unittest

from n_m_cupolaic_prismatoids import PseudoCupolaicPrismatoid as PCP

class TestPCP(unittest.TestCase):
    """Test the class PseudoCupolaicPrismatoid."""
    out_dir = pathlib.Path("tmp")
    exp_dir = pathlib.Path("expected")
    filename = "cupolaic_prismatoid"

    def __init__(self, *args, **opts):
        super().__init__(*args, **opts)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

    # values of n and m that are excluded
    excluded = (
    )

    # values of n and m that are must use crossed squares
    only_crossed = (
        (5, 4),
        (6, 4),
        (7, 6),
        (8, 6),
        (9, 6),
        (9, 8),
        (10, 8),
    )

    # combinations of n and m without double coverage
    no_double_coverage = (
        (6, 2),
        (8, 2),
        (8, 6),
        (9, 6),
        (10, 2),
        (10, 8),
    )

    # combinations of n and m without double coverage
    value_error = (
        (8, 4),
    )

    def test_with_holes(self):
        "Test PCP with some that allow holes andhave double coverage"""

    def test_no_holes(self):
        """Test PCP that where n-grams will not show any holes."""
        for n in range(5, 11):
            for m in range(1, n):
                if (n, m) in self.value_error:
                    continue
                crossed_squares = m % 2 == 1
                if (n, m) in self.only_crossed:
                    crossed_squares = True

                tail= "_crossed_squares" if crossed_squares else ""
                allow_holes = (n, m) in self.no_double_coverage
                if m in (1, n - 1):
                    allow_holes = True
                try:
                    shape = PCP(n, m, crossed_squares, allow_holes)
                except AssertionError as ae:
                    raise ValueError(f"The values n = {n} and m = {m} lead to an assertion") from ae

                # for now compare strings
                filename = f"{n}_{m}_{self.filename}{tail}.off"
                off_str = shape.to_off()
                exp_file = self.exp_dir / filename
                with open(exp_file) as fd:
                    org_str = fd.read()
                if org_str != off_str:
                    # save off file for analysis
                    out_file = self.out_dir / filename
                    with open(out_file, "w") as fd:
                        fd.write(off_str)
                        self.fail(f"PCP({n}, {m}) failed: {out_file} not same as {exp_file}")

if __name__ == "__main__":
    unittest.main()
