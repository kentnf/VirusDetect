import unittest
from unittest import mock

from virusdetect.runtime import resolve_tool


class RuntimeTests(unittest.TestCase):
    def test_resolve_tool_uses_path_only(self):
        with mock.patch("virusdetect.runtime.shutil.which", return_value=None):
            path, source = resolve_tool("bwa")

        self.assertIsNone(path)
        self.assertIsNone(source)


if __name__ == "__main__":
    unittest.main()
