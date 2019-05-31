from optparse import OptionParser, OptionValueError

import xml.etree.ElementTree as ET
from pathlib import Path
import sys

from TMBoundrary import XMLParser


def _check_inputFile(option, opt_str, value, parser):
    f_path = Path(value)
    if not f_path.is_file():
        raise OptionValueError(f"Cannot read input file {f_path.absolute()}")
    setattr(parser.values, option.dest, Path(value))
    parser.values.saved_infile = True

# if __name__ == "__main__":
args = ['-i', './test_data/5y6p_bL.develop201.blast_summ.xml']
options_parser = OptionParser()
options_parser.add_option("-i", "--input", dest="input_xml_filepath", type='str',
                            help="input blast_summ.xml FILE", metavar="FILE",
                            action='callback', callback=_check_inputFile)
(options, args) = options_parser.parse_args(args)
if options.input_xml_filepath is None:
    options_parser.print_help()
    sys.exit()

XML_Info = XMLParser(options.input_xml_filepath)
