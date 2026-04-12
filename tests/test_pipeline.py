import tempfile
import unittest
import json
import io
from pathlib import Path
from types import SimpleNamespace
from contextlib import redirect_stdout
from unittest import mock

import virusdetect.cli as cli
import virusdetect.stages as stages
from virusdetect.pipeline import cleanup_plan_file, cleanup_temp_dir, combine_contigs
from virusdetect.stages import remove_redundant_contigs, run_velvet_assembly
from virusdetect.tools import required_runtime_tools
from virusdetect.sample import remove_duplicate_reads


class PipelineStageTests(unittest.TestCase):
    def test_cleanup_plan_file_removes_plan_by_default(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            plan_path = Path(tmp_dir) / "reads.fa.analysis_plan.json"
            plan_path.write_text("{}\n", encoding="utf-8")

            cleanup_plan_file(plan_path)

            self.assertFalse(plan_path.exists())

    def test_cleanup_plan_file_preserves_plan_when_requested(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            plan_path = Path(tmp_dir) / "reads.fa.analysis_plan.json"
            plan_path.write_text("{}\n", encoding="utf-8")

            cleanup_plan_file(plan_path, keep_file=True)

            self.assertTrue(plan_path.exists())

    def test_cleanup_temp_dir_removes_temp_tree_by_default(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            temp_dir = Path(tmp_dir) / "reads.fa_temp"
            temp_dir.mkdir()
            (temp_dir / "temp.txt").write_text("x\n", encoding="utf-8")

            cleanup_temp_dir(temp_dir)

            self.assertFalse(temp_dir.exists())

    def test_cleanup_temp_dir_preserves_temp_tree_when_requested(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            temp_dir = Path(tmp_dir) / "reads.fa_temp"
            temp_dir.mkdir()
            (temp_dir / "temp.txt").write_text("x\n", encoding="utf-8")

            cleanup_temp_dir(temp_dir, keep_temp=True)

            self.assertTrue(temp_dir.exists())

    def test_combine_contigs_concatenates_aligned_and_assembled_records(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            aligned = tmp / "aligned.fa"
            assembled = tmp / "assembled.fa"
            combined = tmp / "combined.fa"
            aligned.write_text(">ALIGNED1\nAAAA\n", encoding="utf-8")
            assembled.write_text(">NODE_1\nCCCC\n>NODE_2\nGGGG\n", encoding="utf-8")

            contig_count = combine_contigs(aligned, assembled, combined)

            self.assertEqual(contig_count, 3)
            content = combined.read_text(encoding="utf-8")
            self.assertIn(">ALIGNED1", content)
            self.assertIn(">NODE_1", content)
            self.assertIn(">NODE_2", content)

    def test_remove_redundant_contigs_drops_duplicates_and_contained_sequences(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            input_fasta = tmp / "combined.fa"
            output_fasta = tmp / "combined.nr.fa"
            input_fasta.write_text(
                (
                    ">ALIGNED1\nATGCAATGCA\n"
                    ">NODE_1\nATGCAATGCA\n"
                    ">NODE_2\nGCAATG\n"
                    ">NODE_3\nTGCATTGCAT\n"
                    ">NODE_4\nCCCCGGGG\n"
                ),
                encoding="utf-8",
            )

            before_count, after_count = remove_redundant_contigs(input_fasta, output_fasta)

            self.assertEqual(before_count, 5)
            self.assertEqual(after_count, 2)
            content = output_fasta.read_text(encoding="utf-8")
            self.assertIn(">ALIGNED1", content)
            self.assertIn(">NODE_4", content)
            self.assertNotIn(">NODE_1", content)
            self.assertNotIn(">NODE_2", content)
            self.assertNotIn(">NODE_3", content)

    def test_run_pipeline_uses_python_identifier_backend(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            input_path = tmp / "reads.fa"
            input_path.write_text(">r1\nAAAA\n", encoding="utf-8")

            args = cli.build_parser().parse_args(
                [
                    "run",
                    str(input_path),
                    "--backend",
                    "python",
                    "--db-path",
                    str(tmp / "db"),
                    "--output",
                    str(tmp / "out"),
                ]
            )

            with mock.patch("virusdetect.pipeline.verify_database_files", return_value={}):
                with mock.patch("virusdetect.pipeline.check_tools", return_value=[]):
                    with mock.patch("virusdetect.pipeline.missing_required_tools", return_value=[]):
                        with mock.patch("virusdetect.pipeline.run_virus_alignment", return_value=(tmp / "virus.sam", 1)):
                            with mock.patch("virusdetect.pipeline.generate_aligned_contigs") as mock_aligned:
                                with mock.patch("virusdetect.pipeline.run_host_subtraction", return_value=(input_path, 0)):
                                    with mock.patch("virusdetect.pipeline.run_denovo_assembly") as mock_assembly:
                                        with mock.patch("virusdetect.pipeline.remove_redundant_contigs", return_value=(2, 1)):
                                            with mock.patch("virusdetect.pipeline.run_python_identifier") as mock_identify:
                                                aligned_path = tmp / "aligned.fa"
                                                assembled_path = tmp / "assembled.fa"
                                                aligned_path.write_text(">ALIGNED1\nAAAA\n", encoding="utf-8")
                                                assembled_path.write_text(">NODE_1\nCCCC\n", encoding="utf-8")
                                                mock_aligned.return_value = (aligned_path, 1)
                                                mock_assembly.return_value = (assembled_path, 1)
                                                mock_identify.return_value = mock.Mock(
                                                    result_dir=str(tmp / "out" / "result_reads.fa"),
                                                    known_contig_count=1,
                                                    novel_contig_count=0,
                                                    undetermined_contig_count=0,
                                                    summary_tsv=str(tmp / "out" / "result_reads.fa" / "summary.tsv"),
                                                )

                                                exit_code = args.handler(args)

            self.assertEqual(exit_code, 0)
            mock_identify.assert_called_once()
            self.assertFalse((tmp / "out" / "reads.fa_temp").exists())
            self.assertFalse((tmp / "out" / "reads.fa.analysis_plan.json").exists())

    def test_run_pipeline_keep_temp_preserves_sample_temp_dir(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            input_path = tmp / "reads.fa"
            input_path.write_text(">r1\nAAAA\n", encoding="utf-8")

            args = cli.build_parser().parse_args(
                [
                    "run",
                    str(input_path),
                    "--backend",
                    "python",
                    "--db-path",
                    str(tmp / "db"),
                    "--output",
                    str(tmp / "out"),
                    "--keep-temp",
                ]
            )

            with mock.patch("virusdetect.pipeline.verify_database_files", return_value={}):
                with mock.patch("virusdetect.pipeline.check_tools", return_value=[]):
                    with mock.patch("virusdetect.pipeline.missing_required_tools", return_value=[]):
                        with mock.patch("virusdetect.pipeline.run_virus_alignment", return_value=(tmp / "virus.sam", 1)):
                            with mock.patch("virusdetect.pipeline.generate_aligned_contigs") as mock_aligned:
                                with mock.patch("virusdetect.pipeline.run_host_subtraction", return_value=(input_path, 0)):
                                    with mock.patch("virusdetect.pipeline.run_denovo_assembly") as mock_assembly:
                                        with mock.patch("virusdetect.pipeline.remove_redundant_contigs", return_value=(2, 1)):
                                            with mock.patch("virusdetect.pipeline.run_python_identifier") as mock_identify:
                                                aligned_path = tmp / "aligned.fa"
                                                assembled_path = tmp / "assembled.fa"
                                                aligned_path.write_text(">ALIGNED1\nAAAA\n", encoding="utf-8")
                                                assembled_path.write_text(">NODE_1\nCCCC\n", encoding="utf-8")
                                                mock_aligned.return_value = (aligned_path, 1)
                                                mock_assembly.return_value = (assembled_path, 1)
                                                mock_identify.return_value = mock.Mock(
                                                    result_dir=str(tmp / "out" / "result_reads.fa"),
                                                    known_contig_count=1,
                                                    novel_contig_count=0,
                                                    undetermined_contig_count=0,
                                                    summary_tsv=str(tmp / "out" / "result_reads.fa" / "summary.tsv"),
                                                )

                                                exit_code = args.handler(args)

            self.assertEqual(exit_code, 0)
            self.assertTrue((tmp / "out" / "reads.fa_temp").exists())
            self.assertFalse((tmp / "out" / "reads.fa.analysis_plan.json").exists())

    def test_run_pipeline_warns_when_velvet_runs_without_rm_dup(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            input_path = tmp / "reads.fa"
            input_path.write_text(">r1\nAAAA\n", encoding="utf-8")

            args = cli.build_parser().parse_args(
                [
                    "run",
                    str(input_path),
                    "--db-path",
                    str(tmp / "db"),
                    "--output",
                    str(tmp / "out"),
                    "--assembler",
                    "velvet",
                ]
            )

            output = io.StringIO()
            with mock.patch("virusdetect.pipeline.verify_database_files", return_value={}):
                with mock.patch("virusdetect.pipeline.check_tools", return_value=[]):
                    with mock.patch("virusdetect.pipeline.missing_required_tools", return_value=[]):
                        with mock.patch("virusdetect.pipeline.run_virus_alignment", return_value=(tmp / "virus.sam", 1)):
                            with mock.patch("virusdetect.pipeline.generate_aligned_contigs") as mock_aligned:
                                with mock.patch("virusdetect.pipeline.run_host_subtraction", return_value=(input_path, 0)):
                                    with mock.patch("virusdetect.pipeline.run_denovo_assembly") as mock_assembly:
                                        with mock.patch("virusdetect.pipeline.remove_redundant_contigs", return_value=(2, 1)):
                                            with mock.patch("virusdetect.pipeline.run_python_identifier") as mock_identify:
                                                aligned_path = tmp / "aligned.fa"
                                                assembled_path = tmp / "assembled.fa"
                                                aligned_path.write_text(">ALIGNED1\nAAAA\n", encoding="utf-8")
                                                assembled_path.write_text(">NODE_1\nCCCC\n", encoding="utf-8")
                                                mock_aligned.return_value = (aligned_path, 1)
                                                mock_assembly.return_value = (assembled_path, 1)
                                                mock_identify.return_value = mock.Mock(
                                                    result_dir=str(tmp / "out" / "result_reads.fa"),
                                                    known_contig_count=1,
                                                    novel_contig_count=0,
                                                    undetermined_contig_count=0,
                                                    summary_tsv=str(tmp / "out" / "result_reads.fa" / "summary.tsv"),
                                                )

                                                with redirect_stdout(output):
                                                    exit_code = args.handler(args)

            self.assertEqual(exit_code, 0)
            self.assertIn("--assembler velvet", output.getvalue())
            self.assertIn("--rm-dup", output.getvalue())

    def test_run_pipeline_python_backend_plan_omits_legacy_modules(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            input_path = tmp / "reads.fa"
            input_path.write_text(">r1\nAAAA\n", encoding="utf-8")

            args = cli.build_parser().parse_args(
                [
                    "run",
                    str(input_path),
                    "--backend",
                    "python",
                    "--db-path",
                    str(tmp / "db"),
                    "--output",
                    str(tmp / "out"),
                    "--check-only",
                ]
            )

            with mock.patch("virusdetect.pipeline.verify_database_files", return_value={}):
                with mock.patch("virusdetect.pipeline.check_tools", return_value=[]):
                    with mock.patch("virusdetect.pipeline.missing_required_tools", return_value=[]):
                        exit_code = args.handler(args)

            plan_path = tmp / "out" / "reads.fa.analysis_plan.json"
            plan = json.loads(plan_path.read_text(encoding="utf-8"))

        self.assertEqual(exit_code, 0)
        self.assertNotIn("missing_legacy_perl_modules", plan)

    def test_required_runtime_tools_match_spades_default(self):
        required = required_runtime_tools(data_type="sRNA", host_reference=None, assembler="spades")
        self.assertIn("spades.py", required)
        self.assertNotIn("velveth", required)
        self.assertNotIn("hisat2", required)

    def test_required_runtime_tools_include_velvet_backend(self):
        required = required_runtime_tools(data_type="sRNA", host_reference=None, assembler="velvet")
        self.assertIn("velveth", required)
        self.assertIn("velvetg", required)
        self.assertNotIn("spades.py", required)

    def test_required_runtime_tools_include_hisat2_for_mrna_host_subtraction(self):
        required = required_runtime_tools(data_type="mRNA", host_reference="host.fa", assembler="spades")
        self.assertIn("hisat2", required)
        self.assertIn("hisat2-build", required)

    def test_remove_duplicate_reads_keeps_first_fasta_record_per_sequence(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            input_fasta = tmp / "reads.fa"
            output_fasta = tmp / "reads.fa.uniq"
            input_fasta.write_text(">r1\nAAAA\n>r2\nAAAA\n>r3\nCCCC\n", encoding="utf-8")

            kept = remove_duplicate_reads(input_fasta, output_fasta)

            self.assertEqual(kept, 2)
            self.assertEqual(output_fasta.read_text(encoding="utf-8"), ">r1\nAAAA\n>r3\nCCCC\n")

    def test_remove_duplicate_reads_keeps_first_fastq_record_per_sequence(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            input_fastq = tmp / "reads.fq"
            output_fastq = tmp / "reads.fq.uniq"
            input_fastq.write_text("@r1\nAAAA\n+\n####\n@r2\nAAAA\n+\n!!!!\n@r3\nCCCC\n+\n$$$$\n", encoding="utf-8")

            kept = remove_duplicate_reads(input_fastq, output_fastq)

            self.assertEqual(kept, 2)
            self.assertEqual(output_fastq.read_text(encoding="utf-8"), "@r1\nAAAA\n+\n####\n@r3\nCCCC\n+\n$$$$\n")

    def test_run_velvet_assembly_keeps_first_maxlen_match_like_v1(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            sample_path = tmp / "reads.fa"
            sample_path.write_text(">r1\nAAAA\n", encoding="utf-8")
            args = SimpleNamespace(kmer_range="9-13", debug=False, threads=2, rm_dup=False)
            calls = []

            def fake_run_velvet(assembly_input_path, candidate_dir, kmer_length, cov_cutoff, tool_paths, debug):
                calls.append((kmer_length, cov_cutoff, candidate_dir.name))
                candidate_dir.mkdir(parents=True, exist_ok=True)
                contigs_path = candidate_dir / "contigs.fa"
                contigs_path.write_text(f">c_{kmer_length}_{cov_cutoff}\nAAAA\n", encoding="utf-8")
                return contigs_path

            def fake_summarize(contigs_path):
                name = contigs_path.parent.name
                if name.endswith(".9.5"):
                    return {"maxLen": 100.0}
                if name.endswith(".11.5"):
                    return {"maxLen": 100.0}
                if name.endswith(".13.5"):
                    return {"maxLen": 90.0}
                if name.endswith(".9.7"):
                    return {"maxLen": 99.0}
                return {"maxLen": 80.0}

            with mock.patch("virusdetect.stages.tool_path_map", return_value={"velveth": "/bin/velveth", "velvetg": "/bin/velvetg"}):
                with mock.patch("virusdetect.stages.prepare_assembly_input", return_value=sample_path):
                    with mock.patch("virusdetect.stages.run_velvet", side_effect=fake_run_velvet):
                        with mock.patch("virusdetect.stages.summarize_contigs", side_effect=fake_summarize):
                            with mock.patch("virusdetect.stages.count_fasta_sequences", return_value=1):
                                assembled_path, assembled_count = run_velvet_assembly(str(sample_path), args, tmp, "sRNA")

            self.assertEqual(assembled_count, 1)
            self.assertTrue(assembled_path.exists())
            self.assertEqual(calls[0:3], [(9, 5, "reads.fa.velvet.9.5"), (11, 5, "reads.fa.velvet.11.5"), (13, 5, "reads.fa.velvet.13.5")])
            self.assertEqual(calls[3][0:2], (9, 7))
            self.assertEqual(calls[-1][0:2], (9, 5))

    def test_run_velvet_assembly_uses_deduplicated_input_when_requested(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            sample_path = tmp / "reads.fa"
            sample_path.write_text(">r1\nAAAA\n>r2\nAAAA\n", encoding="utf-8")
            args = SimpleNamespace(kmer_range="9-9", debug=False, threads=2, rm_dup=True)
            seen_inputs = []

            def fake_run_velvet(assembly_input_path, candidate_dir, kmer_length, cov_cutoff, tool_paths, debug):
                seen_inputs.append(Path(assembly_input_path))
                candidate_dir.mkdir(parents=True, exist_ok=True)
                contigs_path = candidate_dir / "contigs.fa"
                contigs_path.write_text(">c1\nAAAA\n", encoding="utf-8")
                return contigs_path

            with mock.patch("virusdetect.stages.tool_path_map", return_value={"velveth": "/bin/velveth", "velvetg": "/bin/velvetg"}):
                with mock.patch("virusdetect.stages.prepare_assembly_input", return_value=sample_path):
                    with mock.patch("virusdetect.stages.run_velvet", side_effect=fake_run_velvet):
                        with mock.patch("virusdetect.stages.summarize_contigs", return_value={"maxLen": 100.0}):
                            with mock.patch("virusdetect.stages.count_fasta_sequences", return_value=1):
                                run_velvet_assembly(str(sample_path), args, tmp, "sRNA")

            self.assertTrue(seen_inputs)
            self.assertTrue(seen_inputs[0].name.endswith(".uniq"))
            self.assertEqual(seen_inputs[0].read_text(encoding="utf-8"), ">r1\nAAAA\n")


if __name__ == "__main__":
    unittest.main()
