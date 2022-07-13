import unittest
import utils

class TestUtils(unittest.TestCase):
    def test_smiles_concat_1(self):
        self.assertEqual(utils.smiles_concat(["CO"], ["C"]), ["CO.C"])
    def test_smiles_concat_2(self):
        self.assertEqual(utils.smiles_concat(["CO", "C"], ["C", "CO"]), ["CO.C", "C.CO"])
    def test_smiles_concat_3(self):
        self.assertNotEqual(utils.smiles_concat([""], ["C"]), [".C"])
    def test_smiles_concat_4(self):
        self.assertNotEqual(utils.smiles_concat(["asd"], ["C"]), ["asd.C"])

if __name__ == '__main__':
    unittest.main()