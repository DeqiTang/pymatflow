import unittest

from pymatflow.cpp.cube_handle import cube_diff_1d


class TestCppCucbeHandle(unittest.TestCase):
    """Test cpp.cube_handle """

    @classmethod
    def setUpClass(cls):
        print("A class method called before tests in an individual class are run.\n")

    @classmethod
    def tearDownClass(cls):
        print("A class method called after tests in an individual class have run.\n")

    def setUp(self):
        print("Method called to prepare the test fixture.\n")

    def tearDown(self):
        print("Method called immediately after the test method has been called and the result recorded.\n")

    def test_cube_diff_1d(self):
        """Test method cube_diff_1d"""
        status = cube_diff_1d(["x", "xx", "xxx"], "output.txt", ["c"])
        self.assertTrue(status)

if __name__ == '__main__':

    unittest.main(verbosity=1)