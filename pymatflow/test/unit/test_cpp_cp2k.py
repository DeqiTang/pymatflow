import unittest

from pymatflow.cpp.cp2k import Cp2k

cp2k = Cp2k()

class TestCppCp2k(unittest.TestCase):
    """Test cpp.cp2k """

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

    def test_cp2k(self):
        """Test method add(a, b)"""
        tmp_str = cp2k.to_string()
        #print(tmp_str)
        #self.assertTrue(len(cp2k.to_string()) > 0)
        self.assertTrue(len(tmp_str) > 0)

if __name__ == '__main__':

    unittest.main(verbosity=1)