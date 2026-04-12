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
    cleanup_intermediate_files,
    compute_alignment_identity_text,
    filter_reference_records,
    load_database_annotations,
    prune_redundant_group_ids,
    select_best_hits,
    select_best_raw_hits,
    summarize_hits,
    write_legacy_alignment_sam,
    write_legacy_alignment_table,
    write_reference_detail_pages,
)
from virusdetect.fasta import FastaRecord


class IdentifyTests(unittest.TestCase):
    def test_compute_alignment_identity_text_formats_legacy_style_summary(self):
        hit = BlastHit(
            "blastn",
            "c1",
            100,
            "refA",
            "refA",
            1000,
            "g1",
            "desc",
            7,
            95.0,
            1e-20,
            200.0,
            2,
            8,
            11,
            5,
            70.0,
            "CGT-AAA",
            "CGTTAAA",
        )

        self.assertEqual(compute_alignment_identity_text(hit), "6/7(86%)")

    def test_write_legacy_alignment_outputs(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            xls_path = tmp / "blastn.xls"
            sam_path = tmp / "blastn.sam"
            hits = [
                BlastHit(
                    "blastn",
                    "c1",
                    100,
                    "refA",
                    "refA",
                    1000,
                    "g1",
                    "desc",
                    7,
                    95.0,
                    1e-20,
                    200.0,
                    2,
                    8,
                    11,
                    5,
                    70.0,
                    "CGT-AAA",
                    "CGTTAAA",
                )
            ]
            records_by_id = {"c1": FastaRecord("c1", "AACCGGTT")}

            write_legacy_alignment_table(xls_path, hits, records_by_id)
            write_legacy_alignment_sam(sam_path, hits)

            xls_text = xls_path.read_text(encoding="utf-8")
            sam_text = sam_path.read_text(encoding="utf-8")
            self.assertIn("#Contig_ID\tContig_Seq\tContig_Len\tHit_ID", xls_text)
            self.assertIn("c1\tAACCGGTT\t100\trefA\t1000\tg1\tdesc\t2\t8\t11\t5\t6/7(86%)\t1e-20\t-1", xls_text)
            self.assertIn("@SQ\tSN:refA\tLN:1000", sam_text)
            self.assertIn("c1\t16\trefA\t5\t255\t92H3M1D3M1H\t*\t0\t0\tTTTACG\t*\tAS:i:200\tEV:Z:1e-20", sam_text)

    def test_cleanup_intermediate_files_removes_only_requested_outputs(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            keep_path = tmp / "keep.tsv"
            remove_a = tmp / "temp.a"
            remove_b = tmp / "temp.b"
            keep_path.write_text("depth\n", encoding="utf-8")
            remove_a.write_text("a\n", encoding="utf-8")
            remove_b.write_text("b\n", encoding="utf-8")

            cleanup_intermediate_files([remove_a, remove_b])

            self.assertTrue(keep_path.exists())
            self.assertFalse(remove_a.exists())
            self.assertFalse(remove_b.exists())

    def test_cleanup_intermediate_files_preserves_outputs_in_debug_mode(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            keep_a = tmp / "temp.a"
            keep_b = tmp / "temp.b"
            keep_a.write_text("a\n", encoding="utf-8")
            keep_b.write_text("b\n", encoding="utf-8")

            cleanup_intermediate_files([keep_a, keep_b], keep_files=True)

            self.assertTrue(keep_a.exists())
            self.assertTrue(keep_b.exists())

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
        self.assertIn('id="detail-top"', html_text)
        self.assertIn("Summary report with grouped coverage and identity metrics", html_text)
        self.assertIn("refA", html_text)
        self.assertIn("blastn.reference.fa", html_text)
        self.assertIn('href="#hit-table"', html_text)
        self.assertIn('href="#detail-top"', html_text)
        self.assertIn('<p class="section-nav"><a href="#hit-table">Reference Table</a> | <a href="#detail-top">Top</a></p>', html_text)
        self.assertIn('id="hit-table"', html_text)
        self.assertIn("<h2>Reference Table</h2>", html_text)
        self.assertIn('id="summary-row-refA"', html_text)
        self.assertIn("96.00", html_text)
        self.assertIn("95.00", html_text)
        self.assertIn('href="blastn_references/refA.html"', html_text)
        self.assertIn('href="blastn_references/refA.html"><code>refA</code></a>', html_text)
        self.assertIn('href="blastn_references/refA.html">200 (20.0%)</a>', html_text)
        self.assertIn('href="blastn_references/refA.html">2</a>', html_text)
        self.assertIn('href="blastn_references/refA.html">desc</a>', html_text)
        self.assertIn('colspan="2">Reference</th>', html_text)
        self.assertIn('colspan="6">Coverage</th>', html_text)
        self.assertIn('colspan="3">Identity</th>', html_text)
        self.assertIn('title="Reference accession or grouped identifier for this report row."', html_text)
        self.assertIn(">Map<", html_text)
        self.assertIn("coverage-track coverage-track-compact", html_text)
        self.assertIn("coverage-block forward", html_text)

    def test_build_reference_detail_html_lists_reference_hits(self):
        hits = [
            BlastHit("blastn", "c1", 100, "refA", "refA", 1000, "g1", "desc", 80, 95.0, 1e-20, 200.0, 1, 80, 5, 84, 80.0, "ACGT-AC", "ACGTAAC"),
            BlastHit("blastn", "c2", 120, "refA", "refA", 1000, "g1", "desc", 90, 96.0, 1e-25, 210.0, 1, 90, 91, 180, 75.0, "TTGGCC", "TTGGCA"),
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
        self.assertIn('href="../blastn.html#summary-row-refA"', html_text)
        self.assertIn("Reference Coverage Map", html_text)
        self.assertIn("Reference Length: 1000", html_text)
        self.assertIn('href="#coverage-map"', html_text)
        self.assertIn('href="#hit-table"', html_text)
        self.assertIn('href="#alignment-index"', html_text)
        self.assertIn('href="#alignment-blocks"', html_text)
        self.assertIn('id="detail-top"', html_text)
        self.assertIn('id="coverage-map"', html_text)
        self.assertIn('id="hit-table"', html_text)
        self.assertIn('id="alignment-index"', html_text)
        self.assertIn('id="hit-row-1"', html_text)
        self.assertIn('<td><a href="#alignment-hit-1">1</a></td>', html_text)
        self.assertIn('colspan="3">Contig</th>', html_text)
        self.assertIn('colspan="5">Alignment Coordinates</th>', html_text)
        self.assertIn('title="Start coordinate of the aligned segment on the reference sequence."', html_text)
        self.assertIn("Reference Start", html_text)
        self.assertIn("Reference End", html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">100</a></td><td><a href="#alignment-hit-1-text">1</a></td><td><a href="#alignment-hit-1-text">80</a></td><td><a href="#alignment-hit-1-text">5</a></td><td><a href="#alignment-hit-1-text">84</a></td>', html_text)
        self.assertIn('href="#alignment-hit-1-text">+</a>', html_text)
        self.assertIn("<svg", html_text)
        self.assertIn("reference-panel", html_text)
        self.assertIn("https://www.ncbi.nlm.nih.gov/nuccore/refA", html_text)
        self.assertIn('href="#alignment-hit-1"', html_text)
        self.assertIn('href="#alignment-hit-1-text">5</a>', html_text)
        self.assertIn('href="#alignment-hit-1-text">84</a>', html_text)
        self.assertIn('href="#alignment-hit-1-text">100</a>', html_text)
        self.assertIn('href="#alignment-hit-1-text">+</a>', html_text)
        self.assertIn('href="#alignment-hit-1-text">95.00</a>', html_text)
        self.assertIn('href="#alignment-hit-1-text">1e-20</a>', html_text)
        self.assertIn('href="#alignment-hit-1-text">80.00</a>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1">desc</a></td>', html_text)
        self.assertIn('href="#hit-row-1"', html_text)
        self.assertIn('href="#detail-top"', html_text)
        self.assertIn("Alignment Index", html_text)
        self.assertIn("Quick links into the contig alignments below.", html_text)
        self.assertIn(">Order<", html_text)
        self.assertIn(">Identity Range<", html_text)
        self.assertIn('colspan="1" title="Anchor for each contig block listed in this index.">Contig Block</th>', html_text)
        self.assertIn('colspan="3" title="Alignment count, stable order range, and covered reference span for this contig block.">Block Summary</th>', html_text)
        self.assertIn('colspan="1" title="Percent-identity range across the alignments in this contig block.">Identity</th>', html_text)
        self.assertIn('colspan="1" title="Direct jumps to this contig block, its alignments, and the matching contig row.">Navigation</th>', html_text)
        self.assertIn('title="Number of alignments grouped into this contig block."', html_text)
        self.assertIn('title="Stable order range covered by this contig block."', html_text)
        self.assertIn('title="Reference span covered by this contig block."', html_text)
        self.assertIn('title="Quick jumps to this contig block, its first or last alignment, and the matching contig row."', html_text)
        self.assertIn('href="#alignment-group-1"', html_text)
        self.assertIn('href="#alignment-hit-1-text"', html_text)
        self.assertIn('id="alignment-index-row-1"', html_text)
        self.assertIn('<td><a href="#alignment-group-1">1</a></td>', html_text)
        self.assertIn('<td><a href="#hit-row-1">1</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-group-1">5-84</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">95.00%</a></td>', html_text)
        self.assertIn(">Contig Block<", html_text)
        self.assertIn(">Alignment<", html_text)
        self.assertIn(">Contig Row<", html_text)
        self.assertIn("Detailed Alignments", html_text)
        self.assertIn('id="alignment-blocks"', html_text)
        self.assertIn('id="alignment-group-1"', html_text)
        self.assertIn('id="alignment-hit-1"', html_text)
        self.assertIn('<strong><a href="#alignment-hit-1-text">c1</a></strong> <span class="muted">(<a href="#coverage-map">Coverage Map</a> | <a href="#hit-table">Contig Hits</a> | <a href="#alignment-index-row-1">Index Entry</a> | <a href="#alignment-index">Alignment Index</a> | <a href="#alignment-group-2">Next Contig Block</a> | <a href="#detail-top">Top</a>)</span>', html_text)
        self.assertIn('<strong><a href="#alignment-hit-2-text">c2</a></strong> <span class="muted">(<a href="#coverage-map">Coverage Map</a> | <a href="#hit-table">Contig Hits</a> | <a href="#alignment-index-row-2">Index Entry</a> | <a href="#alignment-index">Alignment Index</a> | <a href="#alignment-group-1">Previous Contig Block</a> | <a href="#detail-top">Top</a>)</span>', html_text)
        self.assertIn('class="alignment-summary"', html_text)
        self.assertIn("Summary:", html_text)
        self.assertIn('<strong>Contig Length:</strong> 100', html_text)
        self.assertIn('<strong>Order Range:</strong> <a href="#hit-row-1">1</a>', html_text)
        self.assertIn("Retained alignments for this contig.", html_text)
        self.assertIn('class="alignment-segment"', html_text)
        self.assertIn('colspan="1">Navigation</th>', html_text)
        self.assertIn('colspan="3">Query</th>', html_text)
        self.assertIn('colspan="2">Reference</th>', html_text)
        self.assertIn('colspan="6">Metrics</th>', html_text)
        self.assertIn('title="Contig identifier for this alignment row."', html_text)
        self.assertIn('title="Stable order of alignments within this contig block."', html_text)
        self.assertIn('title="Aligned length for this contig/reference segment."', html_text)
        self.assertIn('title="Bit score for this alignment."', html_text)
        self.assertIn('title="Percent of the contig query covered by this alignment."', html_text)
        self.assertIn('<strong>Reference Span:</strong> 5-84', html_text)
        self.assertIn('<strong>Identity Range:</strong> 95.00-95.00%', html_text)
        self.assertIn('<strong>Query Coverage Range:</strong> 80.00-80.00%', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text"><code>c1</code></a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">1</a></td><td><a href="#alignment-hit-1-text">80</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">5</a></td><td><a href="#alignment-hit-1-text">84</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">95.00</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">1e-20</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">80</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">200.0</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">80.00</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">+</a></td>', html_text)
        self.assertIn("<strong>Alignment:</strong>", html_text)
        self.assertNotIn("<strong>Alignment</strong></p>", html_text)
        self.assertIn('class="alignment-header"', html_text)
        self.assertIn('class="alignment-metrics muted"', html_text)
        self.assertIn('class="alignment-links muted"', html_text)
        self.assertIn('class="alignment-text"', html_text)
        self.assertIn('class="alignment-line alignment-line-query"', html_text)
        self.assertIn('class="alignment-line alignment-line-match"', html_text)
        self.assertIn('class="alignment-line alignment-line-hit"', html_text)
        self.assertIn("Query 1-80 | Reference 5-84 | Strand +", html_text)
        self.assertIn("<strong>Metrics:</strong>", html_text)
        self.assertIn("<strong>Identity:</strong> 95.00%", html_text)
        self.assertIn("<strong>E-value:</strong> 1e-20", html_text)
        self.assertIn("<strong>Bit Score:</strong> 200.0", html_text)
        self.assertIn("<strong>Alignment Length:</strong> 80", html_text)
        self.assertIn('<p class="alignment-links muted"><a href="#hit-row-1">Contig Hit Row</a> | <a href="#alignment-hit-1">Block Row</a> | <a href="#alignment-index-row-1">Index Entry</a> | <a href="#alignment-index">Alignment Index</a> | <a href="#coverage-map">Coverage Map</a> | <a href="#detail-top">Top</a></p>', html_text)
        self.assertIn("Query  1 ACGT-AC  6", html_text)
        self.assertIn("Match    |||| ||", html_text)
        self.assertIn("Hit    5 ACGTAAC 11", html_text)

    def test_build_reference_detail_html_wraps_alignment_blocks(self):
        long_query = "A" * 100 + "C" * 20
        long_hit = "A" * 100 + "G" * 20
        hits = [
            BlastHit("blastn", "c1", 140, "refA", "refA", 1000, "g1", "desc", 120, 98.0, 1e-30, 250.0, 1, 120, 20, 139, 85.0, long_query, long_hit),
        ]
        summary_row = {
            "group_id": "refA",
            "reference_id": "refA",
            "hit_length": 1000,
            "covered_bases": 120,
            "reference_coverage": 0.12,
            "contig_count": 1,
            "reference_depth": 5.0,
            "normalized_depth": 50.0,
        }

        html_text = build_reference_detail_html("blastn", "refA", "refA", hits, summary_row, "blastn.html")

        self.assertGreaterEqual(html_text.count("Query  "), 2)
        self.assertIn("CCCCCCCCCCCCCCCCCCCC", html_text)
        self.assertIn("GGGGGGGGGGGGGGGGGGGG", html_text)

    def test_build_reference_detail_html_renders_reverse_hit_coordinates_in_alignment_text(self):
        hits = [
            BlastHit("blastn", "c1", 120, "refA", "refA", 1000, "g1", "desc", 7, 95.0, 1e-20, 200.0, 2, 8, 11, 5, 70.0, "CGT-AAA", "CGTTAAA"),
        ]
        summary_row = {
            "group_id": "refA",
            "reference_id": "refA",
            "hit_length": 1000,
            "covered_bases": 7,
            "reference_coverage": 0.007,
            "contig_count": 1,
            "reference_depth": 2.0,
            "normalized_depth": 20.0,
        }

        html_text = build_reference_detail_html("blastn", "refA", "refA", hits, summary_row, "blastn.html")

        self.assertIn("Query 2-8 | Reference 5-11 | Strand -", html_text)
        self.assertIn("<strong>Identity:</strong> 95.00%", html_text)
        self.assertIn("<strong>E-value:</strong> 1e-20", html_text)
        self.assertIn("<strong>Bit Score:</strong> 200.0", html_text)
        self.assertIn("<strong>Alignment Length:</strong> 7", html_text)
        self.assertIn('<p class="alignment-links muted"><a href="#hit-row-1">Contig Hit Row</a> | <a href="#alignment-hit-1">Block Row</a> | <a href="#alignment-index-row-1">Index Entry</a> | <a href="#alignment-index">Alignment Index</a> | <a href="#coverage-map">Coverage Map</a> | <a href="#detail-top">Top</a></p>', html_text)
        self.assertIn("Query  2 CGT-AAA  7", html_text)
        self.assertIn("Match    ||| |||", html_text)
        self.assertIn("Hit   11 CGTTAAA  5", html_text)

    def test_build_reference_detail_html_groups_panel_tracks_by_contig(self):
        hits = [
            BlastHit("blastn", "c1", 120, "refA", "refA", 1000, "g1", "desc", 40, 97.0, 1e-20, 150.0, 1, 40, 10, 49, 33.3, "A" * 40, "A" * 40),
            BlastHit("blastn", "c1", 120, "refA", "refA", 1000, "g1", "desc", 35, 95.0, 1e-10, 120.0, 50, 84, 120, 154, 29.2, "C" * 35, "C" * 35),
            BlastHit("blastn", "c2", 90, "refA", "refA", 1000, "g1", "desc", 30, 99.0, 1e-25, 180.0, 5, 34, 300, 329, 33.3, "G" * 30, "G" * 30),
        ]
        summary_row = {
            "group_id": "refA",
            "reference_id": "refA",
            "hit_length": 1000,
            "covered_bases": 105,
            "reference_coverage": 0.105,
            "contig_count": 2,
            "reference_depth": 6.0,
            "normalized_depth": 60.0,
        }

        html_text = build_reference_detail_html("blastn", "refA", "refA", hits, summary_row, "blastn.html")

        self.assertEqual(html_text.count('class="panel-label">c1</text>'), 1)
        self.assertIn("2 Alignments /", html_text)
        self.assertGreaterEqual(html_text.count('class="panel-segment"'), 3)
        self.assertIn("Alignment Index", html_text)
        self.assertIn('href="#alignment-group-1"', html_text)
        self.assertIn('href="#alignment-hit-1-text"', html_text)
        self.assertIn('href="#hit-row-1"', html_text)
        self.assertIn('id="alignment-index-row-1"', html_text)
        self.assertIn('id="alignment-index-row-2"', html_text)
        self.assertIn('title="Percent-identity range across the alignments in this contig block."', html_text)
        self.assertIn('<td><a href="#alignment-group-1">2</a></td>', html_text)
        self.assertIn('<td><a href="#hit-row-1">1</a>-<a href="#hit-row-2">2</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-group-1">10-154</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">95.00-97.00%</a></td>', html_text)
        self.assertIn(">95.00-97.00%<", html_text)
        self.assertIn(">First Alignment<", html_text)
        self.assertIn(">Last Alignment<", html_text)
        self.assertIn('href="#alignment-hit-2-text"', html_text)
        self.assertEqual(html_text.count('id="alignment-group-1"'), 1)
        self.assertIn('<strong><a href="#alignment-hit-1-text">c1</a></strong> <span class="muted">(<a href="#coverage-map">Coverage Map</a> | <a href="#hit-table">Contig Hits</a> | <a href="#alignment-index-row-1">Index Entry</a> | <a href="#alignment-index">Alignment Index</a> | <a href="#alignment-group-2">Next Contig Block</a> | <a href="#detail-top">Top</a>)</span>', html_text)
        self.assertIn('class="alignment-summary"', html_text)
        self.assertIn("Summary:", html_text)
        self.assertIn("Alignment Segment 1", html_text)
        self.assertIn("Alignment Segment 2", html_text)
        self.assertIn('(Order <a href="#hit-row-1">1</a>)', html_text)
        self.assertIn('(Order <a href="#hit-row-2">2</a>)', html_text)
        self.assertIn('<strong>Order Range:</strong> <a href="#hit-row-1">1</a>-<a href="#hit-row-2">2</a>', html_text)
        self.assertIn("<strong>Alignment:</strong>", html_text)
        self.assertEqual(html_text.count("<strong>Alignment Segment"), 2)
        self.assertGreaterEqual(html_text.count('class="alignment-header"'), 3)
        self.assertGreaterEqual(html_text.count('class="alignment-metrics muted"'), 3)
        self.assertGreaterEqual(html_text.count('class="alignment-links muted"'), 3)
        self.assertGreaterEqual(html_text.count('class="alignment-segment"'), 3)
        self.assertGreaterEqual(html_text.count('class="alignment-text"'), 3)
        self.assertIn('<strong>Contig Length:</strong> 120', html_text)
        self.assertIn('<strong>Alignments:</strong> 2', html_text)
        self.assertIn("Retained alignments for this contig.", html_text)
        self.assertIn('<strong>Reference Span:</strong> 10-154', html_text)
        self.assertIn('<strong>Identity Range:</strong> 95.00-97.00%', html_text)
        self.assertIn('<strong>Query Coverage Range:</strong> 29.20-33.30%', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text"><code>c1</code></a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">120</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">10</a></td><td><a href="#alignment-hit-1-text">49</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-2-text">120</a></td><td><a href="#alignment-hit-2-text">154</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">97.00</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-2-text">95.00</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">1e-20</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-2-text">1e-10</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">40</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-2-text">35</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">150.0</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-2-text">120.0</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">33.30</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-2-text">29.20</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-1-text">+</a></td>', html_text)
        self.assertIn('<td><a href="#alignment-hit-2-text">+</a></td>', html_text)
        self.assertIn("Query 1-40 | Reference 10-49 | Strand +", html_text)
        self.assertIn("Query 50-84 | Reference 120-154 | Strand +", html_text)
        self.assertIn("<strong>Identity:</strong> 97.00%", html_text)
        self.assertIn("<strong>Identity:</strong> 95.00%", html_text)
        self.assertIn("<strong>Bit Score:</strong> 150.0", html_text)
        self.assertIn("<strong>Bit Score:</strong> 120.0", html_text)
        self.assertIn("<strong>Alignment Length:</strong> 40", html_text)
        self.assertIn("<strong>Alignment Length:</strong> 35", html_text)
        self.assertIn('<p class="alignment-links muted"><a href="#hit-row-1">Contig Hit Row</a> | <a href="#alignment-hit-1">Block Row</a> | <a href="#alignment-hit-2-text">Next Segment</a> | <a href="#alignment-index-row-1">Index Entry</a> | <a href="#alignment-index">Alignment Index</a> | <a href="#coverage-map">Coverage Map</a> | <a href="#detail-top">Top</a></p>', html_text)
        self.assertIn('<p class="alignment-links muted"><a href="#hit-row-2">Contig Hit Row</a> | <a href="#alignment-hit-2">Block Row</a> | <a href="#alignment-hit-1-text">Previous Segment</a> | <a href="#alignment-index-row-1">Index Entry</a> | <a href="#alignment-index">Alignment Index</a> | <a href="#coverage-map">Coverage Map</a> | <a href="#detail-top">Top</a></p>', html_text)
        self.assertIn(">Block Row<", html_text)
        self.assertIn(">Previous Contig Block<", html_text)
        self.assertIn(">Next Contig Block<", html_text)
        self.assertIn(">Previous Segment<", html_text)
        self.assertIn(">Next Segment<", html_text)
        self.assertGreaterEqual(html_text.count('href="#alignment-hit-1"'), 2)
        self.assertGreaterEqual(html_text.count('href="#alignment-hit-2"'), 2)
        self.assertGreaterEqual(html_text.count('href="#hit-row-1"'), 3)
        self.assertGreaterEqual(html_text.count('href="#hit-row-2"'), 3)
        self.assertIn('id="alignment-hit-1"', html_text)
        self.assertIn('id="alignment-hit-2"', html_text)

    def test_build_blastx_reference_detail_html_uses_protein_link(self):
        hits = [
            BlastHit("blastx", "c1", 120, "PROT1", "REFSEQ1", 265, "g1", "desc", 62, 40.3, 1.74e-4, 90.0, 1, 62, 6, 67, 51.0, "EVSD---L", "EVADAAAL"),
        ]
        summary_row = {
            "group_id": "PROT1",
            "reference_id": "REFSEQ1",
            "hit_length": 265,
            "covered_bases": 62,
            "reference_coverage": 0.234,
            "contig_count": 1,
            "reference_depth": 29.5,
            "normalized_depth": 295.4,
        }

        html_text = build_reference_detail_html("blastx", "PROT1", "REFSEQ1", hits, summary_row, "blastx.html")

        self.assertIn("Protein ID: <code>PROT1</code>", html_text)
        self.assertIn("Reference ID: <code>REFSEQ1</code>", html_text)
        self.assertIn("https://www.ncbi.nlm.nih.gov/protein/PROT1", html_text)

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
            self.assertIn('href="../blastn.html#summary-row-refA"', detail_html)
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
        self.assertIn('id="detail-top"', html_text)
        self.assertIn("c2", html_text)
        self.assertIn("Total undetermined contigs: 2", html_text)
        self.assertIn('With virus-database matches: <a href="undetermined_blast.html#hit-table">1</a>', html_text)
        self.assertIn('Without virus-database matches: <a href="#no-match-table">1</a>', html_text)
        self.assertIn("Without Virus-Database Matches", html_text)
        self.assertIn('href="#no-match-table"', html_text)
        self.assertIn('href="undetermined_blast.html#hit-table"', html_text)
        self.assertIn('href="#detail-top"', html_text)
        self.assertIn('<p class="section-nav"><a href="#no-match-table">No-Match Table</a> | <a href="undetermined_blast.html#hit-table">Match Table</a> | <a href="#detail-top">Top</a></p>', html_text)
        self.assertIn('id="no-match-table"', html_text)
        self.assertIn("Contig Coverage", html_text)
        self.assertIn('colspan="2" title="Contig identifier and assembled length for this undetermined sequence.">Contig</th>', html_text)
        self.assertIn('colspan="2" title="Read-supported coverage metrics for the undetermined contig.">Coverage</th>', html_text)
        self.assertIn('colspan="2" title="Read-depth metrics for the assembled undetermined contig.">Depth</th>', html_text)
        self.assertIn('title="Contig identifier reported in this analysis."', html_text)
        self.assertIn('title="Average read depth across the covered contig span."', html_text)
        self.assertIn("undetermined_blast.html", html_text)
        self.assertIn('id="undetermined-no-match-row-c2"', html_text)
        self.assertIn('href="#undetermined-no-match-row-c2"><code>c2</code></a>', html_text)
        self.assertNotIn("hit2", html_text)
        self.assertNotIn("With Virus-Database Matches", html_text)
        self.assertIn("100.0%", html_text)

    def test_build_undetermined_blast_report_uses_hit_specific_columns(self):
        undetermined_records = [FastaRecord("c1", "A" * 120)]
        contig_depths = {
            "c1": ContigDepth("c1", 120, 100, 600, 5.0, 50.0),
        }
        raw_hits = [
            BlastHit("blastx", "c1", 120, "hit2", "ref2", 300, "g2", "desc2", 70, 40.0, 1e-20, 110.0, 1, 70, 1, 70, 58.0),
        ]

        html_text = build_undetermined_report_html(
            undetermined_records,
            contig_depths,
            select_best_raw_hits(raw_hits),
            hits_only=True,
            reference_targets={
                ("blastx", "hit2"): {
                    "summary_href": "blastx.html#summary-row-hit2",
                    "detail_href": "blastx_references/hit2.html",
                }
            },
        )

        self.assertIn('Undetermined contigs with virus-database matches: <a href="#hit-table">1</a>', html_text)
        self.assertIn("Virus-database match", html_text)
        self.assertIn('href="#hit-table"', html_text)
        self.assertIn('href="undetermined.html#no-match-table"', html_text)
        self.assertIn('href="#detail-top"', html_text)
        self.assertIn('<p class="section-nav"><a href="#hit-table">Match Table</a> | <a href="undetermined.html#no-match-table">No-Match View</a> | <a href="#detail-top">Top</a></p>', html_text)
        self.assertIn('id="hit-table"', html_text)
        self.assertIn("Genus", html_text)
        self.assertIn('colspan="5" title="Summary columns for the virus-database match on this undetermined contig.">Virus-database match</th>', html_text)
        self.assertIn('colspan="3" title="Additional identifiers and alignment metrics for this virus-database match.">Match metrics</th>', html_text)
        self.assertIn('title="Reference genus annotation from the database metadata."', html_text)
        self.assertIn(">Acc#<", html_text)
        self.assertIn(">E value<", html_text)
        self.assertIn(">Match ID<", html_text)
        self.assertLess(html_text.index(">Acc#<"), html_text.index(">Genus<"))
        self.assertLess(html_text.index(">Genus<"), html_text.index(">Description<"))
        self.assertLess(html_text.index(">Description<"), html_text.index(">E value<"))
        self.assertLess(html_text.index(">E value<"), html_text.index(">BLAST<"))
        self.assertLess(html_text.index(">BLAST<"), html_text.index(">Match ID<"))
        self.assertIn("g2", html_text)
        self.assertIn('id="undetermined-hit-row-c1"', html_text)
        self.assertIn('href="#undetermined-hit-row-c1"><code>c1</code></a>', html_text)
        self.assertIn('href="blastx_references/hit2.html"><code>hit2</code></a>', html_text)
        self.assertIn('href="blastx_references/hit2.html"><code>ref2</code></a>', html_text)
        self.assertIn('href="blastx.html#summary-row-hit2">blastx</a>', html_text)
        self.assertIn('href="blastx_references/hit2.html">desc2</a>', html_text)
        self.assertNotIn("Without Virus-Database Matches", html_text)

    def test_build_undetermined_report_html_renders_sirna_distribution(self):
        undetermined_records = [FastaRecord("c1", "A" * 120)]
        contig_depths = {
            "c1": ContigDepth("c1", 120, 120, 600, 5.0, 50.0),
        }
        html_text = build_undetermined_report_html(
            undetermined_records,
            contig_depths,
            {},
            read_length_stats={"c1": {21: 6, 22: 4, 24: 2}},
            sirna_percent=0.5,
            hits_only=False,
        )

        self.assertIn(">18<", html_text)
        self.assertIn(">33<", html_text)
        self.assertIn("Candidate contigs are highlighted in green and marked in the Candidate column.", html_text)
        self.assertIn('colspan="2" title="Contig identifier and assembled length for this undetermined sequence.">Contig</th>', html_text)
        self.assertIn('colspan="2" title="Read-supported coverage metrics for the undetermined contig.">Coverage</th>', html_text)
        self.assertIn('colspan="17" title="Per-contig read counts across the historical 18-33 nt small-RNA size bins, plus 21-22 nt enrichment.">siRNA size distribution</th>', html_text)
        self.assertIn('colspan="1" title="Candidate flag driven by the configured 21-22 nt enrichment threshold.">Candidate</th>', html_text)
        self.assertIn('colspan="2" title="Read-depth metrics for the assembled undetermined contig.">Depth</th>', html_text)
        self.assertIn("Candidate threshold: &gt;= 50.00%.", html_text)
        self.assertIn('21-22 nt candidate rows: <a href="#candidate-table">1</a>', html_text)
        self.assertIn('href="#candidate-table"', html_text)
        self.assertIn('href="#detail-top"', html_text)
        self.assertIn('<p class="section-nav"><a href="#candidate-table">Candidate Table</a> | <a href="#no-match-table">No-Match Table</a> | <a href="#detail-top">Top</a></p>', html_text)
        self.assertIn('id="candidate-table"', html_text)
        self.assertIn("Candidate Contigs", html_text)
        self.assertIn("Subset of undetermined contigs in this no-match view that pass the configured 21-22 nt enrichment threshold. Contig IDs jump to the matching row in the full no-match table below.", html_text)
        self.assertIn("21-22 (%)", html_text)
        self.assertIn("Candidate", html_text)
        self.assertIn('title="Number of small-RNA reads of length 21 nt assigned to this undetermined contig."', html_text)
        self.assertIn('title="Percent of 18-33 nt reads assigned to this contig that fall in the 21-22 nt range."', html_text)
        self.assertIn('title="Candidate flag based on the configured 21-22 nt enrichment threshold."', html_text)
        self.assertIn("83.33%", html_text)
        self.assertIn('class="candidate-flag candidate-yes"', html_text)
        self.assertIn(">Yes<", html_text)
        self.assertIn('class="candidate-row"', html_text)
        self.assertIn('id="undetermined-candidate-no-match-row-c1"', html_text)
        self.assertIn('id="undetermined-no-match-row-c1"', html_text)
        self.assertIn('href="#undetermined-no-match-row-c1"><code>c1</code></a>', html_text)

    def test_build_undetermined_blast_report_renders_candidate_focus_section(self):
        undetermined_records = [FastaRecord("c1", "A" * 120)]
        contig_depths = {
            "c1": ContigDepth("c1", 120, 100, 600, 5.0, 50.0),
        }
        raw_hits = [
            BlastHit("blastn", "c1", 120, "hit2", "ref2", 300, "g2", "desc2", 70, 98.0, 1e-20, 110.0, 1, 70, 1, 70, 58.0),
        ]

        html_text = build_undetermined_report_html(
            undetermined_records,
            contig_depths,
            select_best_raw_hits(raw_hits),
            read_length_stats={"c1": {21: 6, 22: 4, 24: 2}},
            sirna_percent=0.5,
            hits_only=True,
        )

        self.assertIn('href="#candidate-table"', html_text)
        self.assertIn('href="#detail-top"', html_text)
        self.assertIn('<p class="section-nav"><a href="#candidate-table">Candidate Table</a> | <a href="#hit-table">Match Table</a> | <a href="undetermined.html#no-match-table">No-Match View</a> | <a href="#detail-top">Top</a></p>', html_text)
        self.assertIn('id="candidate-table"', html_text)
        self.assertIn("Candidate Contigs", html_text)
        self.assertIn("Subset of undetermined contigs in this hit view that pass the configured 21-22 nt enrichment threshold. Contig IDs jump to the matching row in the full match table below.", html_text)
        self.assertIn('21-22 nt candidate rows in this view: <a href="#candidate-table">1</a>', html_text)
        self.assertIn('class="candidate-flag candidate-yes"', html_text)
        self.assertIn('id="undetermined-candidate-hit-row-c1"', html_text)
        self.assertIn('id="undetermined-hit-row-c1"', html_text)
        self.assertIn('href="#undetermined-hit-row-c1"><code>c1</code></a>', html_text)


if __name__ == "__main__":
    unittest.main()
