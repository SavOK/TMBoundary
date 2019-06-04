from optparse import OptionParser, OptionValueError
import xml.etree.ElementTree as ET
from pathlib import Path
import subprocess
from subprocess import Popen, PIPE
import sys

from TMBoundrary import XMLParser
from TMBoundrary import RowSQL
from TMBoundrary import PDBParser
from TMBoundrary import Domain


def _set_strucutre_dir(CWD: Path = None,
                       dir_name: str = 'TMStructures', exist_ok: bool = False):
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
test1 = XML_Info.hh_run['hits'][5]['domain_id']
test2 = XML_Info.hh_run['hits'][3]['domain_id']

sql = RowSQL()
domain = Domain()

row1 = sql.get_domain_row(test1)
test_domain = domain.get_structure_path(row1['uid'])
print(test_domain)
pdb = PDBParser(test_domain)
out_file = Path('./test_data/test1.pdb')
pdb.get_region(out_file, 1, 40)

test_dir = Path('/home/saveliy/Projects/TMBoundrary/TM')
if not test_dir.is_dir():
    test_dir.mkdir()

proc1 = subprocess.check_output(
    f"cd {str(test_dir)}; ~rschaeff/bin/generate_pc_pdb.pl 5y6p bL", shell=True)
print(proc1)
