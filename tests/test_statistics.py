import unittest

import pandas as pd

from adipose_senescence_ifn.statistics import (
    benjamini_hochberg,
    permutation_slope,
)


class StatisticsTests(unittest.TestCase):
    def test_animal_permutation_is_exact_and_detects_direction(self):
        frame = pd.DataFrame(
            {
                "animal_id": [f"a{i}" for i in range(6)],
                "age_months": [4, 4, 4, 24, 24, 24],
                "score": [0.0, 0.1, -0.1, 1.0, 1.1, 0.9],
            }
        )
        result = permutation_slope(frame, "score")
        self.assertGreater(result["slope_per_month"], 0)
        self.assertIn("exact", result["permutation_method"])
        self.assertGreater(result["permutation_p"], 0)

    def test_bh_is_monotone_and_bounded(self):
        adjusted = benjamini_hochberg(pd.Series([0.01, 0.04, 0.03, 0.9]))
        self.assertTrue(((adjusted >= 0) & (adjusted <= 1)).all())
        ordered = adjusted.iloc[[0, 2, 1, 3]].to_numpy()
        self.assertTrue(all(a <= b for a, b in zip(ordered, ordered[1:])))


if __name__ == "__main__":
    unittest.main()

