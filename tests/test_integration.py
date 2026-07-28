import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from adipose_senescence_ifn.analysis import run_analysis


class SyntheticIntegrationTests(unittest.TestCase):
    def test_end_to_end_uses_animals_and_writes_outputs(self):
        rng = np.random.default_rng(19)
        genes = [f"G{i}" for i in range(700)]
        rows = []
        cell_ids = []
        expression_rows = []
        for age in (4, 24):
            for animal_index in range(3):
                animal = f"a{age}-{animal_index}"
                for cell_index in range(4):
                    cell = f"{animal}-c{cell_index}"
                    cell_ids.append(cell)
                    values = rng.poisson(3, size=len(genes)).astype(float)
                    if age == 24:
                        values[:20] += 2
                    expression_rows.append(values)
                    rows.append(
                        {
                            "cell_id": cell,
                            "animal_id": animal,
                            "age_months": age,
                            "depot": "SAT",
                            "sex": "male",
                            "cell_type": "Endothelial",
                            "subtype": "EC1",
                            "sample_id": f"{animal}-SAT",
                            "qc_pass": True,
                        }
                    )
        metadata = pd.DataFrame(rows)
        expression = pd.DataFrame(expression_rows, index=cell_ids, columns=genes)
        full = {gene: (0.02 if i % 2 == 0 else -0.01) for i, gene in enumerate(genes[:650])}
        signature = {
            "sencat_variants": {
                "sencat_full": full,
                "sencat_minus_ifng": {k: v for k, v in full.items() if k not in genes[:20]},
                "sencat_minus_mhcii": full,
                "sencat_minus_arrest": full,
                "sencat_residual": {k: v for k, v in full.items() if k not in genes[:60]},
            },
            "binary_axes": {
                "ifng_response": genes[:20],
                "mhcii_antigen_presentation": genes[20:40],
                "cell_cycle_arrest": genes[40:60],
            },
        }

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata_path = root / "metadata.tsv"
            expression_path = root / "expression.tsv"
            signature_path = root / "signatures.json"
            output = root / "results"
            metadata.to_csv(metadata_path, sep="\t", index=False)
            expression.to_csv(expression_path, sep="\t")
            signature_path.write_text(json.dumps(signature), encoding="utf-8")

            summary = run_analysis(
                metadata_path,
                expression_path,
                signature_path,
                output,
                n_permutations=99,
                n_bootstrap=30,
            )
            self.assertTrue(summary["gate_passed"])
            self.assertEqual(summary["animals"], 6)
            self.assertTrue((output / "inference.tsv").exists())
            self.assertTrue((output / "figures/animal-sencat_full.svg").exists())


if __name__ == "__main__":
    unittest.main()

