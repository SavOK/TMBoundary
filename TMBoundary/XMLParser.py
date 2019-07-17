import xml.etree.ElementTree as ET
from pathlib import Path


class ErrorXMLParser(Exception):
    pass


class NumberBlastRuns(ErrorXMLParser):
    def __inti__(self, msg):
        self.msg = msg

    def __repr__(self):
        return(f'Error: {self.msg}')


class XMLParser:
    def __init__(self, filepath: Path):
        self.filepath = filepath
        self._parse_xml()

    def _parse_hit_info(self, xml: ET.Element, recType: str = 'c'):
        if recType not in {'c', 'd', 'h'}:
            raise ValueError
        hit_info = dict()
        hit_info['num'] = int(xml.attrib['num'])
        hit_info['query_reg'] = xml.find('./query_reg').text
        hit_info['hit_reg'] = xml.find('./hit_reg').text
        if recType in {'c', 'd'}:
            hit_info['pdb_id'] = xml.attrib['pdb_id']
            hit_info['chain_id'] = xml.attrib['chain_id']
            hit_info['hsp_count'] = int(xml.attrib['hsp_count'])
            try:
                hit_info['evalues'] = float(xml.attrib['evalues'])
            except ValueError:
                values = [float(x) for x in xml.attrib['evalues'].split(',')]
                hit_info['evalue'] = min(values)
            if recType == 'd':
                hit_info['domain_id'] = xml.attrib['domain_id']
            elif recType == 'c':
                hit_info['query_seq'] = xml.find('./query_seq').text
                hit_info['hit_seq'] = xml.find('./hit_seq').text
        elif recType == 'h':
            hit_info['domain_id'] = xml.attrib['domain_id']
            hit_info['hh_prob'] = float(xml.attrib['hh_prob'])
            hit_info['hh_score'] = float(xml.attrib['hh_score'])
            hit_info['hit_cover'] = float(xml.attrib['hit_cover'])
        return hit_info

    def _parse_chain_blast_run_xml(self, xml: ET.Element):
        self.chain_blast = {}
        self.chain_blast['blast_query'] = xml.find('./blast_query').text
        self.chain_blast['query_len'] = xml.find('./query_len').text
        self.chain_blast['hits'] = [self._parse_hit_info(hit, 'c')
                                    for hit in xml.findall('./*hit')]

    def _parse_domain_blast_run_xml(self, xml: ET.Element):
        self.domain_blast = {}
        self.domain_blast['blast_query'] = xml.find('./blast_query').text
        self.domain_blast['query_len'] = xml.find('./query_len').text
        self.domain_blast['hits'] = [self._parse_hit_info(hit, 'd')
                                     for hit in xml.findall('./*hit')]

    def _parse_hhsearch_run_xml(self, xml: ET.Element):
        self.hh_run = {}
        self.hh_run['hits'] = [self._parse_hit_info(hit, 'h')
                               for hit in xml.findall('./*hit')]

    def _parse_xml(self):
        tree = ET.parse(self.filepath)
        root = tree.getroot()

        blast_summ = root.find('./blast_summ')
        self.pdb = blast_summ.attrib['pdb']
        self.chain = blast_summ.attrib['chain']

        # chain blast parsing
        chain_blast_run_xml = root.findall('./*chain_blast_run')
        if chain_blast_run_xml is None or len(chain_blast_run_xml) == 0:
            raise NumberBlastRuns(
                f'No chain blast found in file {self.filepath.stem}')
        if len(chain_blast_run_xml) != 1:
            errStr = (f'Too many chain blast found ' +
                      f'{len(chain_blast_run_xml)} in ' +
                      f'file {self.filepath.stem}')
            raise NumberBlastRuns(errStr)
        self._parse_chain_blast_run_xml(chain_blast_run_xml[0])

        # domain blast parsing
        domain_blast_run_xml = root.findall('./*blast_run')
        if domain_blast_run_xml is None or len(domain_blast_run_xml) == 0:
            raise NumberBlastRuns(
                f'No domain blast found in file {self.filepath.stem}')
        if len(domain_blast_run_xml) != 1:
            errStr = (f'Too many domain blast found ' +
                      f'{len(domain_blast_run_xml)} in ' +
                      f'file {self.filepath.stem}')
            raise NumberBlastRuns(errStr)
        self._parse_domain_blast_run_xml(domain_blast_run_xml[0])

        # hhsearch run parsing
        hh_run_xml = root.findall('./*hh_run')
        if hh_run_xml is None or len(hh_run_xml) == 0:
            raise NumberBlastRuns(
                f'No hh_run found in file {self.filepath.stem}')
        if len(hh_run_xml) != 1:
            errStr = (f'Too many hhsearch found ' +
                      f'{len(hh_run_xml)} in ' +
                      f'file {self.filepath.stem}')
            raise NumberBlastRuns(errStr)
        self._parse_hhsearch_run_xml(hh_run_xml[0])

    def __repr__(self):
        Str = f"Info from XML {self.filepath.stem},"
        Str += f" pdb: {self.pdb} chain: {self.chain}\n"
        Str += f"number of hits"
        Str += f" chain blast: {len(self.chain_blast['hits'])},"
        Str += f" domain blast: {len(self.domain_blast['hits'])},"
        Str += f" hhsearch: {len(self.hh_run['hits'])}"
        return Str
