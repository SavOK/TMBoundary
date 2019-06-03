from optparse import OptionParser, OptionValueError
import xml.etree.ElementTree as ET
from pathlib import Path
import sys

import psycopg2
import psycopg2.extras

conn_param_dict = {
        "dbname":"ecod", 
        "user":"ecodweb", 
        "host":"129.112.32.63", 
        "port":"45000", 
        "password":"serveruse421",
        }

from TMBoundrary import XMLParser

def _check_inputFile(option, opt_str, value, parser):
    f_path = Path(value)
    if not f_path.is_file():
        raise OptionValueError(f"Cannot read input file {f_path.absolute()}")
    setattr(parser.values, option.dest, Path(value))
    parser.values.saved_infile = True


def get_domain_row(domain_id:str, c):
    c.execute('SELECT * FROM domain where id=%s', (domain_id,))
    results = c.fetchone()
    if results is None:
        print(f"cannot find {domain_id}")
    res_dict = {k:v for k,v in results.items()}
    return res_dict


def get_path_to_domain(domain_uid:str, domain_root:Path=None):
    def _fill_uid(uid:str):
        return '0'*(9-len(uid))+uid
        
    if domain_root is None:
        domain_root=Path("/data/ecod/domain_data")
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
test = XML_Info.hh_run['hits'][5]['domain_id']


conn = psycopg2.connect(**conn_param_dict)
cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
row = get_domain_row(test, cur)
cur.close()
conn.close()

print(get_path_to_domain(row['uid']))
