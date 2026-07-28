import csv
import json
import tempfile
import unittest
from pathlib import Path

from adipose_senescence_ifn.signatures import (
    _reactome_symbols,
    read_mgi_one_to_one,
    read_sencat,
)


class SignatureTests(unittest.TestCase):
    def test_sencat_parser_enforces_shape_and_signs(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "markers.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["converted_alias", "gene", "coef"])
                for index in range(5000):
                    writer.writerow(
                        [f"ENSG{index:011d}", f"GENE{index}", -0.1 if index % 2 else 0.1]
                    )
            weights = read_sencat(path)
            self.assertEqual(len(weights), 5000)
            self.assertEqual(weights["GENE1"], -0.1)

    def test_mgi_retains_only_one_to_one_classes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "homology.tsv"
            columns = ["DB Class Key", "NCBI Taxon ID", "Symbol"]
            rows = [
                ["1", "9606", "CDKN1A"],
                ["1", "10090", "Cdkn1a"],
                ["2", "9606", "GENEA"],
                ["2", "10090", "Genea"],
                ["2", "10090", "Genea2"],
            ]
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(columns)
                writer.writerows(rows)
            mapping, stats = read_mgi_one_to_one(path)
            self.assertEqual(mapping, {"CDKN1A": "Cdkn1a"})
            self.assertEqual(stats["ambiguous_cross_species_classes"], 1)

    def test_reactome_parser_filters_non_gene_entities(self):
        payload = [
            {
                "refEntities": [
                    {
                        "schemaClass": "ReferenceGeneProduct",
                        "displayName": "UniProt:P01579 IFNG",
                    },
                    {
                        "schemaClass": "ReferenceMolecule",
                        "displayName": "ChEBI:123 ATP",
                    },
                ]
            }
        ]
        self.assertEqual(_reactome_symbols(payload), {"IFNG"})


if __name__ == "__main__":
    unittest.main()

