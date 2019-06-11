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
from TMBoundrary import TMalign


def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def _process_range(reg: str):
    # ptr = re.compile(r"\w+[:](?:(?P<start>\d+)[-](?P<end>\d+))")
    ptr = re.compile(r"(?P<start>\d+)[-](?P<end>\d+)")
    match = ptr.findall(reg)
    # if len(match) < 1:
    #     match = ptr_simp.findall(reg)
    if len(match) == 1:
        start = int(match[0][0])-5
        end = int(match[0][1])+5
        if start < 1:
            start = 1
        return [(start, end)]
    start_prev = int(match[0][0])
    end_prev = int(match[0][1])
    region = []
    for reg in match[1:]:
        start_curr = int(reg[0])
        end_curr = int(reg[1])
        if (start_curr - end_prev) < 20:
            end_prev = end_curr
        else:
            start = start_prev - 5
            if start < 1:
                start = 1
            region.append((start, end_prev+5))
            start_prev = start_curr
            end_prev = end_curr
    start = start_prev - 5
    if start < 1:
        start = 1
    region.append((start, end_prev+5))
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


def _set_output_files(infile: Path):
    DIR = infile.parent
    p = infile.name.split('.')
    new_name = f"{p[0]}.{p[1]}.blast_summ_tm.xml"
    return DIR / new_name


def _process_chain_blast(hit: dict, WD: Path, query_structure: Path):
    query_parser = PDBParser(query_structure)
    query_region = _process_range(hit['query_reg'])
    query_reg_filename = f"{query_structure.stem}_{hit['query_reg']}.pdb"
    query_region_file = WD / query_reg_filename
    query_parser.get_region(out_file=query_region_file, regions=query_region)

    hit_gen = ProteinChain(wd=str(WD))
    hit_structure = hit_gen.get_chain_file(pdb=hit['pdb_id'],
                                           chain=hit['chain_id'])
    hit_region = _process_range(hit['hit_reg'])
    hit_reg_filename = f"{hit['pdb_id']}_{hit['chain_id']}_{hit['hit_reg']}.pdb"
    hit_region_file = WD/hit_reg_filename
    hit_parser = PDBParser(hit_structure)
    hit_parser.get_region(out_file=hit_region_file, regions=hit_region)

    TM = TMalign()
    TM_data = TM.run_align(query_region_file, hit_region_file)
    ali_reg = TM_data['ali_reg']
    query_dif = ali_reg[0][0] - query_region[0][0]
    query_reg = ''
    for r in ali_reg:
        query_reg += f'{r[0]-query_dif}-{r[1]-query_dif},'
    query_reg = query_reg.strip(',')
    hit_dif = ali_reg[0][0]-hit_region[0][0]
    hit_reg = ''
    for r in ali_reg:
        hit_reg += f'{r[0]-hit_dif}-{r[1]-hit_dif},'
    hit_reg = hit_reg.strip(',')
    clean_data = {}
    clean_data['pdb_id'] = hit['pdb_id']
    clean_data['chain_id'] = hit['chain_id']
    clean_data['tm_score_norm'] = TM_data['tmN']
    clean_data['tm_score_query'] = TM_data['tm1']
    clean_data['tm_score_hit'] = TM_data['tm2']
    clean_data['query_reg'] = query_reg
    clean_data['hit_reg'] = hit_reg
    clean_data['query_seq'] = TM_data['query_seq']
    clean_data['hit_seq'] = TM_data['hit_seq']
    clean_data['tm_align'] = TM_data['tm_align']
    return(clean_data)


def _process_domain_blast(hit: dict, WD: Path, query_structure: Path):
    query_parser = PDBParser(query_structure)
    query_region = _process_range(hit['query_reg'])
    query_reg_filename = f"{query_structure.stem}_{hit['query_reg']}.pdb"
    query_region_file = WD / query_reg_filename
    query_parser.get_region(out_file=query_region_file, regions=query_region)

    sql = RowSQL()
    domain_info = sql.get_domain_row(hit['domain_id'])
    domain = Domain()
    domain_path = domain.get_structure_path(domain_info['uid'])

    TM = TMalign()
    TM_data = TM.run_align(query_region_file, domain_path)
    ali_reg = TM_data['ali_reg']
    query_dif = ali_reg[0][0] - query_region[0][0]
    query_reg = ''
    for r in ali_reg:
        query_reg += f'{r[0]-query_dif}-{r[1]-query_dif},'
    query_reg = query_reg.strip(',')
    hit_dif = ali_reg[0][0]-1
    hit_reg = ''
    for r in ali_reg:
        hit_reg += f'{r[0]-hit_dif}-{r[1]-hit_dif},'
    hit_reg = hit_reg.strip(',')
    clean_data = {}
    clean_data['domain_id'] = hit['domain_id']
    clean_data['domain_uid'] = domain_info['uid']
    clean_data['tm_score_norm'] = TM_data['tmN']
    clean_data['tm_score_query'] = TM_data['tm1']
    clean_data['tm_score_hit'] = TM_data['tm2']
    clean_data['query_reg'] = query_reg
    clean_data['hit_reg'] = hit_reg
    clean_data['query_seq'] = TM_data['query_seq']
    clean_data['hit_seq'] = TM_data['hit_seq']
    clean_data['tm_align'] = TM_data['tm_align']
    return(clean_data)


def _process_domain_hh(hit: dict, WD: Path, query_structure: Path):
    query_parser = PDBParser(query_structure)
    query_region = _process_range(hit['query_reg'])
    query_reg_filename = f"{query_structure.stem}_{hit['query_reg']}.pdb"
    query_region_file = WD / query_reg_filename
    query_parser.get_region(out_file=query_region_file, regions=query_region)

    sql = RowSQL()
    domain_info = sql.get_domain_row(hit['domain_id'])
    domain = Domain()
    domain_path = domain.get_structure_path(domain_info['uid'])

    TM = TMalign()
    TM_data = TM.run_align(query_region_file, domain_path)
    ali_reg = TM_data['ali_reg']
    query_dif = ali_reg[0][0] - query_region[0][0]
    query_reg = ''
    for r in ali_reg:
        query_reg += f'{r[0]-query_dif}-{r[1]-query_dif},'
    query_reg = query_reg.strip(',')
    hit_dif = ali_reg[0][0]-1
    hit_reg = ''
    for r in ali_reg:
        hit_reg += f'{r[0]-hit_dif}-{r[1]-hit_dif},'
    hit_reg = hit_reg.strip(',')
    clean_data = {}
    clean_data['domain_id'] = hit['domain_id']
    clean_data['domain_uid'] = domain_info['uid']
    clean_data['tm_score_norm'] = TM_data['tmN']
    clean_data['tm_score_query'] = TM_data['tm1']
    clean_data['tm_score_hit'] = TM_data['tm2']
    clean_data['query_reg'] = query_reg
    clean_data['hit_reg'] = hit_reg
    clean_data['query_seq'] = TM_data['query_seq']
    clean_data['hit_seq'] = TM_data['hit_seq']
    clean_data['tm_align'] = TM_data['tm_align']
    return(clean_data)


def _create_domain_xml(hit, num):
    hit_xml = ET.Element('hit')
    hit_xml.attrib['num'] = num
    hit_xml.attrib['domain_id'] = hit['domain_id']
    hit_xml.attrib['domain_uid'] = hit['domain_uid']
    hit_xml.attrib['tm_score_query'] = hit['tm_score_query']
    hit_xml.attrib['tm_score_hit'] = hit['tm_score_hit']
    hit_xml.attrib['tm_score_norm'] = hit['tm_score_norm']
    query_reg = ET.SubElement(hit_xml, 'query_reg')
    query_reg.text = hit['query_reg']
    hit_reg = ET.SubElement(hit_xml, 'hit_reg')
    hit_reg.text = hit['hit_reg']
    query_seq = ET.SubElement(hit_xml, 'query_seq')
    query_seq.text = hit['query_seq']
    ali_seq = ET.SubElement(hit_xml, 'tm_align')
    ali_seq.text = hit['tm_align']
    hit_seq = ET.SubElement(hit_xml, 'hit_seq')
    hit_seq.text = hit['hit_seq']
    return hit_xml


def create_XML(old, new, Info):
    xml = ET.parse(str(old))
    root = xml.getroot()
    ET.SubElement(root, 'tm_summary')
    tm = root.find('tm_summary')
    ET.SubElement(tm, 'tm_blast_domain')
    tm_blast_domain = tm.find('tm_blast_domain')
    for ix, I in enumerate(Info['blast_domain']):
        ET.SubElement(tm_blast_domain, 'hit')
        hit = tm_blast_domain[ix]
        hit.attrib['num'] = str(ix+1)
        hit.attrib['domain_id'] = str(I['domain_id'])
        hit.attrib['domain_uid'] = str(I['domain_uid'])
        hit.attrib['tm_score_query'] = str(I['tm_score_query'])
        hit.attrib['tm_score_hit'] = str(I['tm_score_hit'])
        hit.attrib['tm_score_norm'] = str(I['tm_score_norm'])
        ET.SubElement(hit, 'query_reg')
        query_reg = hit.find('query_reg')
        query_reg.text = str(I['query_reg'])
        ET.SubElement(hit, 'hit_reg')
        hit_reg = hit.find('hit_reg')
        hit_reg.text = str(I['hit_reg'])
        ET.SubElement(hit, 'query_seq')
        query_seq = hit.find('query_seq')
        query_seq.text = I['query_seq']
        ET.SubElement(hit, 'align_seq')
        ali_seq = hit.find('align_seq')
        ali_seq.text = str(I['tm_align'])
        ET.SubElement(hit, 'hit_seq')
        hit_seq = hit.find('hit_seq')
        hit_seq.text = str(I['hit_seq'])
    ET.SubElement(tm, 'tm_hh_domain')
    tm_hh_domain = tm.find('tm_hh_domain')
    for ix, I in enumerate(Info['hh_domain']):
        ET.SubElement(tm_hh_domain, 'hit')
        hit = tm_hh_domain[ix]
        hit.attrib['num'] = str(ix+1)
        hit.attrib['domain_id'] = str(I['domain_id'])
        hit.attrib['domain_uid'] = str(I['domain_uid'])
        hit.attrib['tm_score_query'] = str(I['tm_score_query'])
        hit.attrib['tm_score_hit'] = str(I['tm_score_hit'])
        hit.attrib['tm_score_norm'] = str(I['tm_score_norm'])
        ET.SubElement(hit, 'query_reg')
        query_reg = hit.find('query_reg')
        query_reg.text = str(I['query_reg'])
        ET.SubElement(hit, 'hit_reg')
        hit_reg = hit.find('hit_reg')
        hit_reg.text = str(I['hit_reg'])
        ET.SubElement(hit, 'query_seq')
        query_seq = hit.find('query_seq')
        query_seq.text = I['query_seq']
        ET.SubElement(hit, 'align_seq')
        ali_seq = hit.find('align_seq')
        ali_seq.text = str(I['tm_align'])
        ET.SubElement(hit, 'hit_seq')
        hit_seq = hit.find('hit_seq')
        hit_seq.text = str(I['hit_seq'])
    indent(root)
    xml.write(str(new))
    return root


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


# Setup structure directory
str_dir = _set_strucutre_dir(options.work_dir)
# Setup classes
protein_gen = ProteinChain(wd=str(str_dir))  # protein generator
# Parse XML summary
XML_Info = XMLParser(options.input_xml_filepath)
# Set query structure
query_structure = protein_gen.get_chain_file(XML_Info.pdb, XML_Info.chain)
quary_parser = PDBParser(query_structure)

Info = {}
# Process blast_chain
Info['blast_chain'] = []
for hit in XML_Info.chain_blast['hits'][-10:]:
    Info['blast_chain'].append(
        _process_chain_blast(hit, str_dir, query_structure))

Info['blast_domain'] = []
for hit in XML_Info.domain_blast['hits'][-10:]:
    Info['blast_domain'].append(
        _process_domain_blast(hit, str_dir, query_structure))

Info['hh_domain'] = []
for hit in XML_Info.hh_run['hits'][-10:]:
    Info['hh_domain'].append(_process_domain_hh(hit, str_dir, query_structure))

out_file = _set_output_files(options.input_xml_filepath)
print(out_file)
t = create_XML(options.input_xml_filepath, out_file, Info)
# _process_range(range_test_1)

# domain = Domain()
# path_to_domain = domain.get_structure_path(row_data['uid'])
# protein = ProteinChain(str_dir)


# print(test_domain)
# pdb = PDBParser(test_domain)
# out_file = Path('./test_data/test1.pdb')
# pdb.get_region(out_file, 1, 40)
