import unittest
from unittest import mock

from virusdetect.runtime import ResolvedTool, missing_required_tools, resolve_tool


class RuntimeTests(unittest.TestCase):
    def test_resolve_tool_uses_path_only(self):
        with mock.patch("virusdetect.runtime.shutil.which", return_value=None):
            path, source = resolve_tool("bwa")

        self.assertIsNone(path)
        self.assertIsNone(source)

    def test_missing_required_tools_uses_default_required_flags(self):
        statuses = [
            ResolvedTool(name="bwa", required=True, path=None, source=None),
            ResolvedTool(name="spades.py", required=False, path=None, source=None),
        ]

        self.assertEqual(missing_required_tools(statuses), ["bwa"])

    def test_missing_required_tools_accepts_runtime_specific_required_names(self):
        statuses = [
            ResolvedTool(name="bwa", required=True, path="/bin/bwa", source="PATH"),
            ResolvedTool(name="spades.py", required=False, path=None, source=None),
            ResolvedTool(name="velveth", required=False, path=None, source=None),
        ]

        self.assertEqual(missing_required_tools(statuses, required_names=("bwa", "spades.py")), ["spades.py"])


if __name__ == "__main__":
    unittest.main()
