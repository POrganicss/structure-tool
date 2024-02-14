import unittest

class TestLog(unittest.TestCase):

    def test_getxyzs(self):
        log = """
        Some log content...
        Standard orientation:
        Atom 1  Atom 2  Atom 3  Atom 4  Atom 5  Atom 6
        C       0.0000  0.0000  0.0000
        H       1.0000  0.0000  0.0000
        H       0.0000  1.0000  0.0000
        H       0.0000  0.0000  1.0000
        C       1.0000  1.0000  1.0000
        H       2.0000  1.0000  1.0000
        H       1.0000  2.0000  1.0000
        H       1.0000  1.0000  2.0000
        """
        expected_xyzs = [
            "C 0.0000 0.0000 0.0000\nH 1.0000 0.0000 0.0000\nH 0.0000 1.0000 0.0000\nH 0.0000 0.0000 1.0000\n",
            "C 1.0000 1.0000 1.0000\nH 2.0000 1.0000 1.0000\nH 1.0000 2.0000 1.0000\nH 1.0000 1.0000 2.0000\n"
        ]

        result = Log.getxyzs(log)

        self.assertEqual(len(result), 2)
        self.assertEqual(result, expected_xyzs)

    def test_getxyzs_no_coordinates(self):
        log = """
        Some log content...
        No coordinates found.
        """

        result = Log.getxyzs(log)

        self.assertEqual(len(result), 0)
        self.assertEqual(result, [])

if __name__ == '__main__':
    unittest.main()