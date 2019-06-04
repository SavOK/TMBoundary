from pathlib import Path


class DomainFileError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        return self.msg


class Domain:
    def __init__(self, domain_root_path: Path = None):
        if domain_root_path is None:
            domain_root_path = Path("/data/ecod/domain_data")
        self.root_path = domain_root_path

    def _fill_uid(self, uid: str):
        uid = str(uid)
        return '0'*(9-len(uid))+uid

    def _get_path_to_domain_location(self, domain_uid: str):
        str_id = self._fill_uid(domain_uid)
        dirpath = self.root_path/str_id[2:7]/str_id
        return dirpath

    def get_structure_path(self, uid: str):
        pathtostructure = self._get_path_to_domain_location(uid)
        filepath = pathtostructure / f'{pathtostructure.name}.pdb'
        if filepath.is_file():
            return filepath
        filepath = pathtostructure / f'{pathtostructure.name}.pdbnum.pdb'
        if filepath.is_file():
            return filepath
        filepath = pathtostructure / f'{pathtostructure.name}.csv.pdb'
        if filepath.is_file():
            return filepath
        filepath = pathtostructure / f'{pathtostructure.name}.seqres.pdb'
        if filepath.is_file():
            return filepath
        raise DomainFileError(f'cannot find pdb file for domain {uid}')
