import tempfile
import unittest
from pathlib import Path

import pandas as pd

from adipose_senescence_ifn.gate import audit_metadata


def valid_metadata(animals_per_age=3):
    rows = []
    for age in (4, 24):
        for index in range(animals_per_age):
            animal = f"a{age}-{index}"
            for depot in ("SAT", "VAT"):
                for cell in range(2):
                    rows.append(
                        {
                            "cell_id": f"{animal}-{depot}-{cell}",
                            "animal_id": animal,
                            "age_months": age,
                            "depot": depot,
                            "sex": "male",
                            "cell_type": "Endothelial",
                            "subtype": "Vascular EC",
                            "sample_id": f"{animal}-{depot}",
                            "qc_pass": True,
                        }
                    )
    return pd.DataFrame(rows)


class GateTests(unittest.TestCase):
    def test_valid_animal_resolved_metadata_passes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "metadata.tsv"
            valid_metadata().to_csv(path, sep="\t", index=False)
            report, metadata = audit_metadata(path)
            self.assertTrue(report.passed)
            self.assertEqual(metadata["animal_id"].nunique(), 6)

    def test_under_replicated_design_stops(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "metadata.tsv"
            valid_metadata(animals_per_age=2).to_csv(path, sep="\t", index=False)
            report, _ = audit_metadata(path)
            self.assertFalse(report.passed)
            failed = {x["name"] for x in report.checks if not x["passed"]}
            self.assertIn("animal_level_replication", failed)

    def test_missing_file_stops_without_exception(self):
        report, metadata = audit_metadata(Path("/nonexistent/metadata.tsv"))
        self.assertFalse(report.passed)
        self.assertIsNone(metadata)


if __name__ == "__main__":
    unittest.main()

