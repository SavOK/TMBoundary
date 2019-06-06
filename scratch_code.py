from optparse import OptionParser, OptionValueError
import xml.etree.ElementTree as ET
from pathlib import Path
import re
import subprocess
from subprocess import Popen, PIPE
import sys

from TMBoundrary import XMLParser
from TMBoundrary import RowSQL
from TMBoundrary import PDBParser
from TMBoundrary import Domain
from TMBoundrary import ProteinChain

def _process_range(reg: str):
    ptr = re.compile(r"\w+[:](?:(?P<start>\d+)[-](?P<end>\d+))")
    match = ptr.findall(reg)
    if len(match) == 1:
        start = int(match[0][0])
        end = int(match[0][0])
        return [(start-5, end+5)]
    start_prev = int(match[0][0])
    end_prev = int(match[0][1])
    region = []
    for reg in match[1:]:
        start_curr = int(reg[0])
        end_curr = int(reg[1])
        if (start_curr - end_prev) < 20:
            end_prev = end_curr
        else:
            region.append((start_prev-5, end_prev+5))
            start_prev = start_curr
            end_prev = end_curr
    region.append((start_prev-5, end_prev+5))
    return region

def _set_strucutre_dir(CWD: Path = None,
                       dir_name: str = None, exist_ok: bool = False):
    if dir_name is None:
        dir_name = 'TMstructures'
    if CWD is None:
        CWD = Path().absolute()
    str_dir = CWD / dir_name
    if not str_dir.is_dir():
        str_dir.mkdir(parents=True, exist_ok=exist_ok)
    return str_dir


def _check_inputFile(option, opt_str, value, parser):
    f_path = Path(value)
    if not f_path.is_file():
        raise OptionValueError(f"Cannot get {str(f_path)} file")
    setattr(parser.values, option.dest, Path(f_path))
    parser.values.saved_infile = True


def _check_inputDir(option, opt_str, value, parser):
    f_path = Path(value)
    if not f_path.is_dir():
        f_path = f_path.expanduser()
        if not f_path.is_dir():
            raise OptionValueError(f"Cannot get {str(f_path)} directory")
    setattr(parser.values, option.dest, Path(f_path))
    parser.values.saved_infile = True


# if __name__ == "__main__":
args = ['-i', './test_data/5y6p_bL.develop201.blast_summ.xml',
        '-w', '~/Projects/TMBoundrary']
options_parser = OptionParser()
options_parser.add_option("-i", "--input", dest="input_xml_filepath", type='str',
                          help="input blast_summ.xml FILE", metavar="FILE",
                          action='callback', callback=_check_inputFile)
options_parser.add_option("-w", "--work_dir", dest="work_dir", type='str',
                          help="DIR where structure files will be stored $DIR/TMfiles", metavar="DIR",
                          action='callback', callback=_check_inputDir, default=None)
(options, args) = options_parser.parse_args(args)
if options.input_xml_filepath is None:
    options_parser.print_help()
    sys.exit()

XML_Info = XMLParser(options.input_xml_filepath)
str_dir = _set_strucutre_dir(options.work_dir)

domain_id = XML_Info.hh_run['hits'][5]['domain_id']
test2 = XML_Info.hh_run['hits'][3]

sql = RowSQL()
row_data = sql.get_domain_row(domain_id)


test_pdb1 = Path('C:\\Users\\Saveliy Belkin\\Projects\\TMBoundrary\\test_data\\001720163.pdb') 
test_pdb2 = Path('C:\\Users\\Saveliy Belkin\\Projects\\TMBoundrary\\test_data\\out_file.pdb') 

with Popen(['/home/saveliy/prog/TMalign/TMalign', str(test_pdb1), str(test_pdb2)]) as proc:
    print(proc)
# _process_range(range_test_1)

# domain = Domain()
# path_to_domain = domain.get_structure_path(row_data['uid'])
# protein = ProteinChain(str_dir)


# print(test_domain)
# pdb = PDBParser(test_domain)
# out_file = Path('./test_data/test1.pdb')
# pdb.get_region(out_file, 1, 40)
