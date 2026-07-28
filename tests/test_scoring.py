import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from adipose_senescence_ifn.scoring import (
    aggregate_animal_scores,
    score_programs,
)


class ScoringTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(7)
        self.genes = [f"G{i}" for i in range(700)]
        self.cells = [f"c{i}" for i in range(18)]
        self.expression = pd.DataFrame(
            rng.poisson(3, size=(len(self.cells), len(self.genes))),
            index=self.cells,
            columns=self.genes,
        )
        full = {gene: (-0.01 if i % 2 else 0.02) for i, gene in enumerate(self.genes[:650])}
        self.signature = {
            "sencat_variants": {
                "sencat_full": full,
                "sencat_residual": {k: v for k, v in full.items() if k not in self.genes[:20]},
            },
            "binary_axes": {
                "ifng_response": self.genes[:20],
                "mhcii_antigen_presentation": self.genes[20:40],
                "cell_cycle_arrest": self.genes[40:60],
            },
        }

    def test_scores_and_animal_aggregation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "signatures.json"
            path.write_text(json.dumps(self.signature), encoding="utf-8")
            scores, coverage = score_programs(self.expression, path)
            self.assertEqual(scores.shape, (18, 5))
            self.assertTrue((coverage["coverage_fraction"] > 0).all())

            metadata = pd.DataFrame(
                {
                    "cell_id": self.cells,
                    "animal_id": [f"a{i // 3}" for i in range(18)],
                    "age_months": [4 if i < 9 else 24 for i in range(18)],
                    "depot": ["SAT"] * 18,
                    "sex": ["male"] * 18,
                    "cell_type": ["Endothelial"] * 18,
                    "subtype": ["EC1"] * 18,
                }
            ).set_index("cell_id")
            aggregate = aggregate_animal_scores(scores, metadata)
            self.assertEqual(len(aggregate), 6)
            self.assertTrue((aggregate["n_cells"] == 3).all())


if __name__ == "__main__":
    unittest.main()

