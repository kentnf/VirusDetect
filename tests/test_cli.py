import io
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

import virusdetect.cli as cli
from virusdetect.db import LEGACY_MARKER_FILES, LEGACY_REFERENCE_FILES
from virusdetect.legacy import (
    build_legacy_command,
    build_legacy_identifier_command,
    create_legacy_tool_bridge,
    legacy_execution_env,
)


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
        self.assertIn("pixi add bwa samtools blast", output.getvalue())
        self.assertIn("mamba install -c conda-forge -c bioconda", output.getvalue())
        self.assertIn("host subtraction and de novo assembly", output.getvalue())

    def test_run_parser_accepts_multiple_inputs_and_prepare_flag(self):
        args = cli.build_parser().parse_args(["run", "a.fa", "b.fa", "--prepare-only"])
        self.assertEqual(args.input_paths, ["a.fa", "b.fa"])
        self.assertTrue(args.prepare_only)
        self.assertEqual(args.hisat_dist, 5)
        self.assertEqual(args.diff_ratio, 0.25)
        self.assertEqual(args.diff_contig_cover, 0.5)
        self.assertEqual(args.diff_contig_length, 100)

    def test_build_legacy_command_maps_core_args(self):
        args = cli.build_parser().parse_args(
            [
                "run",
                "a.fa",
                "--host-reference",
                "host.fa",
                "--read-length",
                "21-23",
                "--debug",
            ]
        )
        command = build_legacy_command(args, [Path("a.fa")])
        self.assertIn("--reference", command)
        self.assertIn("vrl_plant", command)
        self.assertIn("--host_reference", command)
        self.assertIn("host.fa", command)
        self.assertIn("--read_length", command)
        self.assertIn("21-23", command)
        self.assertIn("--norm_depth_cutoff", command)
        self.assertIn("5.0", command)
        self.assertIn("--diff-ratio", command)
        self.assertIn("0.25", command)
        self.assertIn("--diff-contig-cover", command)
        self.assertIn("0.5", command)
        self.assertIn("--diff-contig-length", command)
        self.assertIn("100", command)
        self.assertIn("--debug", command)

    def test_run_parser_defaults_to_legacy_backend(self):
        args = cli.build_parser().parse_args(["run", "a.fa"])
        self.assertEqual(args.backend, "legacy")

    def test_build_legacy_identifier_command_uses_resolved_reference(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            db = tmp / "db"
            db.mkdir()
            (db / "vrl_plant").write_text(">ref\nAAAA\n", encoding="utf-8")
            args = cli.build_parser().parse_args(["run", "reads.fa"])
            command = build_legacy_identifier_command(args, str(db), Path("reads.fa"), Path("contigs.fa"))

        self.assertIn("--reference", command)
        self.assertIn(str((db / "vrl_plant").resolve()), command)
        self.assertIn("--norm-depth-cutoff", command)
        self.assertIn("--diff-ratio", command)
        self.assertIn("--diff-contig-cover", command)
        self.assertIn("--diff-contig-length", command)
        self.assertTrue(command[-2].endswith("reads.fa"))
        self.assertTrue(command[-1].endswith("contigs.fa"))

    def test_create_legacy_tool_bridge_writes_expected_wrappers(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            with mock.patch(
                "virusdetect.legacy.resolve_tool",
                side_effect=[
                    ("/opt/bin/bwa", "PATH"),
                    ("/opt/bin/samtools", "PATH"),
                    ("/opt/bin/blastn", "PATH"),
                    ("/opt/bin/blastx", "PATH"),
                    ("/opt/bin/makeblastdb", "PATH"),
                ],
            ):
                bridge_dir = create_legacy_tool_bridge(tmp)
            self.assertTrue((bridge_dir / "bwa").exists())
            self.assertTrue((bridge_dir / "samtools").exists())
            self.assertTrue((bridge_dir / "megablast").exists())
            self.assertTrue((bridge_dir / "blastall").exists())
            self.assertTrue((bridge_dir / "formatdb").exists())

    def test_legacy_execution_env_prepends_bridge_dir(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            with mock.patch("virusdetect.legacy.create_legacy_tool_bridge", return_value=tmp / ".virusdetect-legacy-tools"):
                with mock.patch.dict("os.environ", {"PATH": "/usr/bin"}, clear=True):
                    env = legacy_execution_env(tmp)

        self.assertEqual(env["VIRUSDETECT_TOOL_DIR"], str(tmp / ".virusdetect-legacy-tools"))
        self.assertEqual(env["PATH"], f"{tmp / '.virusdetect-legacy-tools'}:/usr/bin")


if __name__ == "__main__":
    unittest.main()
