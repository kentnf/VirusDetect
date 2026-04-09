import gzip
import tempfile
import unittest
from pathlib import Path

from virusdetect.identify import (
    BlastHit,
    ContigDepth,
    apply_summary_filter_flags,
    build_reference_detail_html,
    build_reference_report_html,
    build_group_support,
    build_undetermined_report_html,
    filter_reference_records,
    load_database_annotations,
    prune_redundant_group_ids,
    select_best_hits,
    select_best_raw_hits,
    summarize_hits,
    write_reference_detail_pages,
)
from virusdetect.fasta import FastaRecord


class IdentifyTests(unittest.TestCase):
    def test_load_database_annotations_maps_protein_ids(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            with gzip.open(tmp / "vrl_genbank.info.gz", "wt", encoding="utf-8") as handle:
                handle.write("AB123\t1000\tvirusgenus\tExample virus genome\n")
            with gzip.open(tmp / "vrl_idmapping.gz", "wt", encoding="utf-8") as handle:
                handle.write("AB123\tPROT1\n")
            (tmp / "vrl_plant").write_text(">AB123\nAAAA\n", encoding="utf-8")
            (tmp / "vrl_plant_prot").write_text(">PROT1\nMMMM\n", encoding="utf-8")

            annotations, protein_to_reference = load_database_annotations(str(tmp))

        self.assertEqual(annotations["AB123"].genus, "virusgenus")
        self.assertEqual(protein_to_reference["PROT1"], "AB123")

    def test_select_best_hits_keeps_best_hit_per_contig(self):
        hits = [
            BlastHit("blastn", "c1", 100, "h1", "h1", 1000, "g1", "d1", 80, 95.0, 1e-20, 200.0, 1, 80, 1, 80, 80.0),
            BlastHit("blastn", "c1", 100, "h2", "h2", 1000, "g2", "d2", 82, 95.0, 1e-30, 180.0, 1, 82, 1, 82, 82.0),
            BlastHit("blastn", "c2", 100, "h3", "h3", 1000, "g3", "d3", 40, 99.0, 1e-50, 300.0, 1, 40, 1, 40, 40.0),
        ]

        selected = select_best_hits(hits, min_identity=60.0, min_query_coverage=50.0)

        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0].contig_id, "c1")
        self.assertEqual(selected[0].hit_id, "h2")

    def test_summarize_hits_groups_by_analysis_and_reference(self):
        hits = [
            BlastHit("blastn", "c1", 100, "h1", "ref1", 1000, "g1", "d1", 80, 95.0, 1e-20, 200.0, 1, 80, 1, 80, 80.0),
            BlastHit("blastn", "c2", 120, "h2", "ref1", 1000, "g1", "d1", 90, 96.0, 1e-25, 210.0, 1, 90, 1, 90, 75.0),
            BlastHit("blastx", "c3", 200, "p1", "ref2", 200, "g2", "d2", 120, 45.0, 1e-6, 150.0, 5, 124, 10, 49, 60.0),
        ]
        contig_depths = {
            "c1": ContigDepth("c1", 100, 80, 400, 4.0, 4.0),
            "c2": ContigDepth("c2", 120, 75, 480, 4.0, 4.0),
            "c3": ContigDepth("c3", 200, 60, 100, 0.5, 0.5),
        }

        summary = summarize_hits(
            hits,
            contig_depths=contig_depths,
            library_size=1_000_000,
            coverage_cutoff=0.1,
            depth_cutoff=5.0,
            norm_depth_cutoff=5.0,
        )

        self.assertEqual(len(summary), 2)
        self.assertEqual(summary[0]["analysis"], "blastn")
        self.assertEqual(summary[0]["contig_count"], 2)
        self.assertFalse(summary[0]["passes_coverage_depth"])
        self.assertFalse(summary[0]["passes_filters"])
        self.assertEqual(summary[1]["analysis"], "blastx")
        self.assertEqual(summary[1]["group_id"], "p1")
        self.assertEqual(summary[1]["reference_id"], "ref2")

    def test_redundancy_pruning_keeps_only_stronger_reference(self):
        hits = [
            BlastHit("blastn", "c1", 100, "refA", "refA", 1000, "g1", "d1", 100, 95.0, 1e-20, 200.0, 1, 100, 1, 100, 100.0),
            BlastHit("blastn", "c2", 100, "refA", "refA", 1000, "g1", "d1", 100, 95.0, 1e-20, 200.0, 1, 100, 201, 300, 100.0),
            BlastHit("blastn", "c3", 100, "refA", "refA", 1000, "g1", "d1", 100, 95.0, 1e-20, 200.0, 1, 100, 401, 500, 100.0),
            BlastHit("blastn", "c4", 100, "refA", "refA", 1000, "g1", "d1", 100, 95.0, 1e-20, 200.0, 1, 100, 601, 700, 100.0),
            BlastHit("blastn", "c1", 100, "refB", "refB", 1000, "g1", "d1", 100, 95.0, 1e-20, 200.0, 1, 100, 1, 100, 100.0),
            BlastHit("blastn", "c2", 100, "refB", "refB", 1000, "g1", "d1", 100, 95.0, 1e-20, 200.0, 1, 100, 201, 300, 100.0),
            BlastHit("blastn", "c3", 100, "refB", "refB", 1000, "g1", "d1", 100, 95.0, 1e-20, 200.0, 1, 100, 401, 500, 100.0),
            BlastHit("blastn", "c5", 100, "refB", "refB", 1000, "g1", "d1", 40, 95.0, 1e-20, 180.0, 1, 40, 901, 940, 40.0),
        ]
        contig_depths = {
            contig_id: ContigDepth(contig_id, 100, 100, 1000, 10.0, 10.0)
            for contig_id in ("c1", "c2", "c3", "c4", "c5")
        }

        summary = summarize_hits(
            hits,
            contig_depths=contig_depths,
            library_size=1_000_000,
            coverage_cutoff=0.1,
            depth_cutoff=5.0,
            norm_depth_cutoff=5.0,
        )
        kept_group_ids, redundant_group_ids = prune_redundant_group_ids(summary, build_group_support(hits))
        summary = apply_summary_filter_flags(summary, kept_group_ids)
        summary_by_group = {row["group_id"]: row for row in summary}

        self.assertEqual(kept_group_ids, {"refA"})
        self.assertEqual(redundant_group_ids, {"refB"})
        self.assertTrue(summary_by_group["refA"]["passes_filters"])
        self.assertTrue(summary_by_group["refB"]["passes_coverage_depth"])
        self.assertFalse(summary_by_group["refB"]["passes_redundancy"])
        self.assertFalse(summary_by_group["refB"]["passes_filters"])

    def test_filter_reference_records_matches_normalized_accessions(self):
        records = [
            FastaRecord("lcl|ref|AB123.1| Example nucleotide", "AAAA"),
            FastaRecord("sp|P99999| Example protein", "MMMM"),
            FastaRecord("plain_id", "TTTT"),
        ]

        selected = filter_reference_records(records, {"AB123.1", "plain_id"})

        self.assertEqual([record.seq_id for record in selected], ["lcl|ref|AB123.1|", "plain_id"])

    def test_build_reference_report_html_contains_kept_reference(self):
        summary_rows = [
            {
                "analysis": "blastn",
                "group_id": "refA",
                "reference_id": "refA",
                "hit_length": 1000,
                "genus": "g1",
                "description": "desc",
                "contig_count": 2,
                "supporting_contigs": "c1,c2",
                "total_contig_length": 220,
                "covered_bases": 200,
                "reference_coverage": 0.2,
                "reference_depth": 8.0,
                "normalized_depth": 80.0,
                "mean_percent_identity": 95.5,
                "best_evalue": 1e-20,
                "max_query_coverage": 100.0,
                "passes_coverage_depth": True,
                "passes_redundancy": True,
                "passes_filters": True,
            }
        ]
        final_hits = [
            BlastHit("blastn", "c1", 100, "refA", "refA", 1000, "g1", "desc", 80, 95.0, 1e-20, 200.0, 1, 80, 1, 80, 80.0),
            BlastHit("blastn", "c2", 120, "refA", "refA", 1000, "g1", "desc", 90, 96.0, 1e-25, 210.0, 1, 90, 91, 180, 75.0),
        ]

        html_text = build_reference_report_html("blastn", summary_rows, final_hits, "blastn.reference.fa", "blastn_references")

        self.assertIn("BLASTN Reference Summary", html_text)
        self.assertIn("refA", html_text)
        self.assertIn("blastn.reference.fa", html_text)
        self.assertIn("96.00", html_text)
        self.assertIn("95.00", html_text)
        self.assertIn('href="blastn_references/refA.html"', html_text)
        self.assertIn(">Map<", html_text)
        self.assertIn("coverage-track coverage-track-compact", html_text)
        self.assertIn("coverage-block forward", html_text)

    def test_build_reference_detail_html_lists_reference_hits(self):
        hits = [
            BlastHit("blastn", "c1", 100, "refA", "refA", 1000, "g1", "desc", 80, 95.0, 1e-20, 200.0, 1, 80, 5, 84, 80.0),
            BlastHit("blastn", "c2", 120, "refA", "refA", 1000, "g1", "desc", 90, 96.0, 1e-25, 210.0, 1, 90, 91, 180, 75.0),
        ]
        summary_row = {
            "group_id": "refA",
            "reference_id": "refA",
            "hit_length": 1000,
            "covered_bases": 170,
            "reference_coverage": 0.17,
            "contig_count": 2,
            "reference_depth": 8.0,
            "normalized_depth": 80.0,
        }

        html_text = build_reference_detail_html("blastn", "refA", "refA", hits, summary_row, "blastn.html")

        self.assertIn("BLASTN Reference Detail: refA", html_text)
        self.assertIn("c1", html_text)
        self.assertIn("c2", html_text)
        self.assertIn("95.00", html_text)
        self.assertIn("1e-20", html_text)
        self.assertIn('href="../blastn.html"', html_text)
        self.assertIn("Reference Coverage Map", html_text)
        self.assertIn("Reference Length: 1000", html_text)
        self.assertIn("5-84", html_text)
        self.assertIn(">+<", html_text)

    def test_write_reference_detail_pages_writes_only_kept_groups(self):
        summary_rows = [
            {
                "group_id": "refA",
                "reference_id": "refA",
                "hit_length": 1000,
                "covered_bases": 170,
                "reference_coverage": 0.17,
                "contig_count": 1,
                "reference_depth": 8.0,
                "normalized_depth": 80.0,
                "passes_filters": True,
            },
            {
                "group_id": "refB",
                "reference_id": "refB",
                "hit_length": 1000,
                "covered_bases": 120,
                "reference_coverage": 0.12,
                "contig_count": 1,
                "reference_depth": 7.0,
                "normalized_depth": 70.0,
                "passes_filters": False,
            },
        ]
        hits = [
            BlastHit("blastn", "c1", 100, "refA", "refA", 1000, "g1", "desc", 80, 95.0, 1e-20, 200.0, 1, 80, 5, 84, 80.0),
            BlastHit("blastn", "c2", 120, "refB", "refB", 1000, "g1", "desc", 90, 96.0, 1e-25, 210.0, 1, 90, 91, 180, 75.0),
        ]

        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir)

            detail_dir_name = write_reference_detail_pages(output_dir, "blastn", summary_rows, hits)

            detail_dir = output_dir / detail_dir_name
            self.assertEqual(detail_dir_name, "blastn_references")
            self.assertTrue((detail_dir / "refA.html").exists())
            self.assertFalse((detail_dir / "refB.html").exists())
            detail_html = (detail_dir / "refA.html").read_text(encoding="utf-8")
            self.assertIn('href="../blastn.html"', detail_html)
            self.assertIn("Reference Coverage Map", detail_html)
            self.assertIn("c1", detail_html)

    def test_build_undetermined_report_html_uses_best_raw_hit(self):
        undetermined_records = [FastaRecord("c1", "A" * 120), FastaRecord("c2", "C" * 90)]
        contig_depths = {
            "c1": ContigDepth("c1", 120, 120, 600, 5.0, 50.0),
            "c2": ContigDepth("c2", 90, 90, 90, 1.0, 10.0),
        }
        raw_hits = [
            BlastHit("blastn", "c1", 120, "hit1", "ref1", 1000, "g1", "desc1", 60, 90.0, 1e-10, 100.0, 1, 60, 1, 60, 50.0),
            BlastHit("blastx", "c1", 120, "hit2", "ref2", 300, "g2", "desc2", 70, 40.0, 1e-20, 110.0, 1, 70, 1, 70, 58.0),
        ]

        html_text = build_undetermined_report_html(
            undetermined_records,
            contig_depths,
            select_best_raw_hits(raw_hits),
            hits_only=False,
        )

        self.assertIn("Undetermined Contigs", html_text)
        self.assertIn("c1", html_text)
        self.assertIn("hit2", html_text)
        self.assertIn("blastx", html_text)
        self.assertIn("c2", html_text)
        self.assertIn("Total undetermined contigs: 2", html_text)
        self.assertIn("With Virus-Database Hits", html_text)
        self.assertIn("Without Virus-Database Hits", html_text)
        self.assertIn("Contig Coverage", html_text)
        self.assertIn("100.0%", html_text)


if __name__ == "__main__":
    unittest.main()
