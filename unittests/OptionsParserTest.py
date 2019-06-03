# from pathlib import Path
# import unittest
# from TMBoundrary import OptionsParser as OP


# class TestOptionParser(unittest.TestCase):
#     def testCorrect(self):
#         corr_args = ['-i', './test_data/5y6p_bL.develop201.blast_summ.xml']
#         expect_val = Path(corr_args[1])
#         options = OP(corr_args)
#         self.assertIsInstance(options.input_xml_filename, Path)
#         self.assertEqual(str(options.input_xml_filename),
#                          str(expect_val.absolute()))
    
#     def testEmpty(self):
#        pass 


# if __name__ == "__main__":
#     unittest.main()
