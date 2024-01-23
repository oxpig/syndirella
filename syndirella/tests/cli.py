import unittest
import subprocess

class TestCLI(unittest.TestCase):
    def test_cli(self):
        subprocess.check_call(['syndirella', '--help'])