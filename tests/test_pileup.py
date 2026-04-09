import tempfile
import unittest
from pathlib import Path

from virusdetect.pileup import count_fasta_sequences, filter_pileup_by_coverage, pileup_to_contig
from virusdetect.stages import parse_kmer_values


class PileupTests(unittest.TestCase):
    def test_parse_kmer_values_returns_odd_values(self):
        self.assertEqual(parse_kmer_values("9-15"), [9, 11, 13, 15])

    def test_filter_pileup_by_coverage_keeps_only_supported_reference(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            pileup = tmp / "input.pileup"
            seq_info = tmp / "seq.info"
            output = tmp / "output.pileup"
            pileup.write_text(
                "\n".join(
                    [
                        "ref1\t1\tA\t1\t.\tI",
                        "ref1\t2\tA\t1\t.\tI",
                        "ref1\t3\tA\t1\t.\tI",
                        "ref1\t4\tA\t1\t.\tI",
                        "ref2\t1\tA\t1\t.\tI",
                        ""
                    ]
                ),
                encoding="utf-8",
            )
            seq_info.write_text(
                "ref1\t10\tg\td\t1\th\nref2\t10\tg\td\t1\th\n",
                encoding="utf-8",
            )

            filter_pileup_by_coverage(pileup, seq_info, 0.3, output)
            content = output.read_text(encoding="utf-8")

        self.assertIn("ref1\t1\tA\t1\t.\tI", content)
        self.assertNotIn("ref2\t1\tA\t1\t.\tI", content)

    def test_pileup_to_contig_writes_consensus_fasta(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            reference = tmp / "reference.fa"
            pileup = tmp / "input.pileup"
            output = tmp / "aligned.fa"
            reference.write_text(">ref1\nAAAAA\n", encoding="utf-8")
            pileup.write_text(
                "\n".join(
                    [
                        "ref1\t1\tA\t3\t...\tIII",
                        "ref1\t2\tA\t3\t...\tIII",
                        "ref1\t3\tA\t3\t...\tIII",
                        ""
                    ]
                ),
                encoding="utf-8",
            )

            contig_num = pileup_to_contig(reference, pileup, output, 2, 0, "ALIGNED")
            content = output.read_text(encoding="utf-8")
            fasta_count = count_fasta_sequences(output)

        self.assertEqual(contig_num, 1)
        self.assertIn(">ALIGNED1 ref1 3", content)
        self.assertIn("AAA", content)
        self.assertEqual(fasta_count, 1)


if __name__ == "__main__":
    unittest.main()
