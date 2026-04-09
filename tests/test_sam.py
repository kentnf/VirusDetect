import tempfile
import unittest
from pathlib import Path

from virusdetect.sam import expand_xa_hits, filter_sam_best_hits, reverse_complement


class SamUtilityTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ACGTacgt"), "acgtACGT")

    def test_expand_xa_hits_adds_secondary_records(self):
        sam_content = (
            "read1\t0\tref1\t1\t37\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0\tXA:Z:ref2,+5,4M,1;\n"
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            sam_path = Path(tmp_dir) / "test.sam"
            sam_path.write_text(sam_content, encoding="utf-8")
            expand_xa_hits(sam_path)
            lines = sam_path.read_text(encoding="utf-8").strip().splitlines()

        self.assertEqual(len(lines), 2)
        self.assertIn("\tref2\t5\t0\t4M\t", lines[1])
        self.assertTrue(lines[1].endswith("NM:i:1"))

    def test_filter_sam_best_hits_keeps_lowest_nm_only(self):
        sam_content = "\n".join(
            [
                "read1\t0\tref1\t1\t37\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:1",
                "read1\t0\tref2\t2\t37\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0",
                "read2\t4\t*\t0\t0\t*\t*\t0\t0\tTGCA\tIIII",
                "read3\t0\tref3\t3\t37\t4M\t*\t0\t0\tTGCA\tIIII\tNM:i:2",
                "",
            ]
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            sam_path = Path(tmp_dir) / "test.sam"
            sam_path.write_text(sam_content, encoding="utf-8")
            mapped_num = filter_sam_best_hits(sam_path)
            lines = sam_path.read_text(encoding="utf-8").strip().splitlines()

        self.assertEqual(mapped_num, 1)
        self.assertEqual(len(lines), 1)
        self.assertIn("\tref2\t2\t37\t4M\t", lines[0])


if __name__ == "__main__":
    unittest.main()
