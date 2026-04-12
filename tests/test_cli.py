import argparse
import io
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

import virusdetect.cli as cli
from virusdetect.db import LEGACY_MARKER_FILES, LEGACY_REFERENCE_FILES


class VirusDetectCliTests(unittest.TestCase):
    def test_normalize_argv_strips_leading_separator(self):
        self.assertEqual(cli.normalize_argv(["--", "db", "path"]), ["db", "path"])
        self.assertEqual(cli.normalize_argv(["run", "reads.fa"]), ["run", "reads.fa"])

    def test_db_path_uses_resolved_location(self):
        args = cli.build_parser().parse_args(["db", "path"])

        with mock.patch("virusdetect.cli.resolve_database_location", return_value=mock.Mock(path="/tmp/db", source="env", kind="legacy")):
            output = io.StringIO()
            with redirect_stdout(output):
                exit_code = cli.handle_db(args)

        self.assertEqual(exit_code, 0)
        self.assertEqual(output.getvalue().strip(), "/tmp/db")

    def test_db_verify_accepts_legacy_layout(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            for file_name in LEGACY_MARKER_FILES + LEGACY_REFERENCE_FILES:
                (tmp_path / file_name).write_text("x", encoding="utf-8")

            args = cli.build_parser().parse_args(["db", "verify", "--path", tmp_dir])
            output = io.StringIO()
            with redirect_stdout(output):
                exit_code = cli.handle_db(args)

        self.assertEqual(exit_code, 0)
        self.assertIn("Database verification passed", output.getvalue())

    def test_db_bundle_uses_resolved_database_and_default_output(self):
        args = cli.build_parser().parse_args(["db", "bundle"])

        with mock.patch("virusdetect.cli.resolve_database_location", return_value=mock.Mock(path="/tmp/db", source="cwd", kind="legacy")):
            with mock.patch(
                "virusdetect.cli.bundle_database_archive",
                return_value=mock.Mock(
                    archive_path="/tmp/dist/virusdetect-db-2.0.0a0.tar.gz",
                    sha256_path="/tmp/dist/virusdetect-db-2.0.0a0.tar.gz.sha256",
                ),
            ) as mock_bundle:
                output = io.StringIO()
                with redirect_stdout(output):
                    exit_code = cli.handle_db(args)

        self.assertEqual(exit_code, 0)
        self.assertEqual(mock_bundle.call_args.kwargs["source_dir"], "/tmp/db")
        self.assertTrue(mock_bundle.call_args.kwargs["destination_archive"].endswith("dist/virusdetect-db-2.0.0a0.tar.gz"))
        self.assertEqual(mock_bundle.call_args.kwargs["db_version"], "2.0.0a0")
        self.assertIn("Database bundle created", output.getvalue())
        self.assertIn("SHA256 file written", output.getvalue())

    def test_db_download_uses_default_release_source(self):
        args = cli.build_parser().parse_args(["db", "download"])

        with mock.patch("virusdetect.cli.resolve_database_target_dir", return_value="/tmp/install") as mock_target:
            with mock.patch(
                "virusdetect.cli.resolve_database_download_source",
                return_value=mock.Mock(
                    url="https://github.com/kentnf/VirusDetect/releases/download/v2.0.0a0/virusdetect-db-2026.04.tar.gz",
                    sha256=None,
                    sha256_url="https://github.com/kentnf/VirusDetect/releases/download/v2.0.0a0/virusdetect-db-2026.04.tar.gz.sha256",
                ),
            ) as mock_source:
                with mock.patch("virusdetect.cli.install_database_archive", return_value="/tmp/install") as mock_install:
                    output = io.StringIO()
                    with redirect_stdout(output):
                        exit_code = cli.handle_db(args)

        self.assertEqual(exit_code, 0)
        mock_target.assert_called_once_with(None)
        self.assertEqual(mock_source.call_args.kwargs["db_version"], "2026.04")
        self.assertEqual(mock_install.call_args.kwargs["destination_dir"], "/tmp/install")
        self.assertIn("Database source:", output.getvalue())
        self.assertIn("Database installed to: /tmp/install", output.getvalue())

    def test_tools_check_json(self):
        args = cli.build_parser().parse_args(["tools", "check", "--json"])
        output = io.StringIO()
        with mock.patch("virusdetect.cli.check_tools", return_value=[]):
            with redirect_stdout(output):
                exit_code = cli.handle_tools(args)

        self.assertEqual(exit_code, 0)
        self.assertEqual(output.getvalue().strip(), "[]")

    def test_tools_install_hint_outputs_conda_guidance(self):
        args = cli.build_parser().parse_args(["tools", "install-hint"])
        output = io.StringIO()
        with redirect_stdout(output):
            exit_code = cli.handle_tools_install_hint(args)

        self.assertEqual(exit_code, 0)
        self.assertIn("Recommended pixi setup for the Python backend:", output.getvalue())
        self.assertIn("pixi add bwa samtools blast hisat2 spades velvet", output.getvalue())
        self.assertIn("mamba install -c conda-forge -c bioconda", output.getvalue())
        self.assertIn("alternate assembler via `--assembler velvet`", output.getvalue())
        self.assertIn("default `main` workflow no longer needs Perl packages", output.getvalue())
        self.assertIn("historical Perl workflow, use the `v1` branch", output.getvalue())
        self.assertNotIn("--legacy", output.getvalue())

    def test_run_parser_accepts_multiple_inputs_and_prepare_flag(self):
        args = cli.build_parser().parse_args(["run", "a.fa", "b.fa", "--prepare-only"])
        self.assertEqual(args.input_paths, ["a.fa", "b.fa"])
        self.assertTrue(args.prepare_only)
        self.assertEqual(args.hisat_dist, 5)
        self.assertEqual(args.diff_ratio, 0.25)
        self.assertEqual(args.diff_contig_cover, 0.5)
        self.assertEqual(args.diff_contig_length, 100)

    def test_run_parser_defaults_to_python_backend(self):
        args = cli.build_parser().parse_args(["run", "a.fa"])
        self.assertEqual(args.backend, "python")
        self.assertEqual(args.assembler, "spades")

    def test_run_parser_accepts_velvet_assembler(self):
        args = cli.build_parser().parse_args(["run", "a.fa", "--assembler", "velvet"])
        self.assertEqual(args.assembler, "velvet")

    def test_run_parser_accepts_rm_dup_aliases(self):
        args_dash = cli.build_parser().parse_args(["run", "a.fa", "--rm-dup"])
        args_underscore = cli.build_parser().parse_args(["run", "a.fa", "--rm_dup"])
        self.assertTrue(args_dash.rm_dup)
        self.assertTrue(args_underscore.rm_dup)

    def test_run_parser_rejects_legacy_backend(self):
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            with self.assertRaises(SystemExit):
                cli.build_parser().parse_args(["run", "a.fa", "--backend", "legacy"])
        self.assertIn("Legacy backend has been removed from `main`", stderr.getvalue())

    def test_run_parser_rejects_unknown_backend(self):
        stderr = io.StringIO()
        with redirect_stderr(stderr):
            with self.assertRaises(SystemExit):
                cli.build_parser().parse_args(["run", "a.fa", "--backend", "perl"])
        self.assertIn("Unsupported backend: perl. Use `python`.", stderr.getvalue())

    def test_run_help_shows_python_only_interface(self):
        parser = cli.build_parser()
        run_parser = next(
            action.choices["run"]
            for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)
        )
        help_text = run_parser.format_help()
        self.assertNotIn("--backend", help_text)
        self.assertNotIn("{legacy,python}", help_text)

if __name__ == "__main__":
    unittest.main()
