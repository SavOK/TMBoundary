from optparse import OptionParser, OptionValueError
import xml.etree.ElementTree as ET
from pathlib import Path
import sys


from TMBoundrary import XMLParser
from TMBoundrary import RowSQL


def _check_inputFile(option, opt_str, value, parser):
    f_path = Path(value)
    if not f_path.is_file():
        raise OptionValueError(f"Cannot read input file {f_path.absolute()}")
    setattr(parser.values, option.dest, Path(value))
    parser.values.saved_infile = True


def get_path_to_domain(domain_uid: str, domain_root: Path = None):
    def _fill_uid(uid: str):
        return '0'*(9-len(uid))+uid

    if domain_root is None:
        domain_root = Path("/data/ecod/domain_data")
    domain_uid = str(domain_uid)
    str_id = _fill_uid(domain_uid)
    dirpath = domain_root/str_id[2:7]/str_id
    return dirpath



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
test1 = XML_Info.hh_run['hits'][5]['domain_id']
test2 = XML_Info.hh_run['hits'][3]['domain_id']

sql = RowSQL()
row1 = sql.get_domain_row(test1)
row2 = sql.get_domain_row(test2)

print(row1)
domain_path1 = get_path_to_domain(row1['uid'])
print(row2)
domain_path2 = get_path_to_domain(row2['uid'])
