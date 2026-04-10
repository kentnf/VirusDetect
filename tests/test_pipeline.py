import tempfile
import unittest
from contextlib import redirect_stdout
import io
from pathlib import Path
from unittest import mock

import virusdetect.cli as cli
from virusdetect.pipeline import combine_contigs
from virusdetect.stages import remove_redundant_contigs


class PipelineStageTests(unittest.TestCase):
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
                                                    summary_tsv=str(tmp / "out" / "result_reads.fa" / "reads.fa.summary.tsv"),
                                                )

                                                exit_code = args.handler(args)

            self.assertEqual(exit_code, 0)
            mock_identify.assert_called_once()

    def test_run_pipeline_warns_when_using_legacy_backend(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            input_path = tmp / "reads.fa"
            input_path.write_text(">r1\nAAAA\n", encoding="utf-8")

            args = cli.build_parser().parse_args(
                [
                    "run",
                    str(input_path),
                    "--backend",
                    "legacy",
                    "--db-path",
                    str(tmp / "db"),
                ]
            )

            output = io.StringIO()
            with mock.patch("virusdetect.pipeline.verify_database_files", return_value={}):
                with mock.patch("virusdetect.pipeline.check_tools", return_value=[]):
                    with mock.patch("virusdetect.pipeline.missing_required_tools", return_value=[]):
                        with mock.patch("virusdetect.pipeline.missing_legacy_perl_modules", return_value=[]):
                            with mock.patch("virusdetect.pipeline.run_legacy_pipeline", return_value=0) as mock_legacy:
                                with redirect_stdout(output):
                                    exit_code = args.handler(args)

        self.assertEqual(exit_code, 0)
        mock_legacy.assert_called_once()
        self.assertIn("Legacy backend is deprecated on `main`", output.getvalue())
        self.assertIn("use the `v1` branch", output.getvalue())


if __name__ == "__main__":
    unittest.main()
