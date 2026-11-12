"""Regression tests for retaining the final junction during merging."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from joinONTLRs import JoinLR


def junction(position):
    return ("chr1", position, "+", "D", "chr1", position + 200, "+", "D")


def evidence(name):
    return [[name], [{}, {}, {}, {}, {}], [60], 0, ["300"], [600], ["0"],
            ["start"], ["1"], ["AAA"], ["CCC"]]


class FinalJunctionTests(unittest.TestCase):
    def merge(self, positions):
        caller = JoinLR.__new__(JoinLR)
        caller.maxDistToJoin = 30
        caller.allFusions = {"chr1": {junction(pos): evidence(str(pos)) for pos in positions}}
        _, merged = caller.joinChromSimilarLRs("chr1")
        return merged

    def test_single_junction_is_retained(self):
        merged = self.merge([100])
        self.assertEqual(list(merged), [junction(100)])
        self.assertEqual(merged[junction(100)][0], ["100"])

    def test_distinct_terminal_junction_is_retained(self):
        self.assertEqual(set(self.merge([100, 1000])), {junction(100), junction(1000)})

    def test_similar_terminal_junction_is_merged_once(self):
        merged = self.merge([100, 110])
        self.assertEqual(list(merged), [junction(100)])
        self.assertEqual(merged[junction(100)][0], ["100", "110"])

    def test_empty_junctions_write_header_without_division_by_zero(self):
        caller = JoinLR.__new__(JoinLR)
        caller.allFusionsJoined1 = {"chr1": {}}
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "junctions.tsv"
            caller.writeFusionToOutput(str(output))
            lines = output.read_text().splitlines()
            self.assertEqual(len(lines), 1)
            self.assertTrue(lines[0].startswith("Chrom1\tPos1\t"))


if __name__ == "__main__":
    unittest.main()
