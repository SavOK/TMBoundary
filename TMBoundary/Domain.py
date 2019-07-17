from pathlib import Path
import subprocess


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

    def _grep_file(self, uid: str):
        path_to_strc = self._get_path_to_domain_location(uid)
        grep_str = ('grep' + ' -El ' + f'\'^.{{13}}(CA )\' ' +
                    f"{str(path_to_strc)}/*pdb")
        proc = subprocess.check_output(grep_str,
                                       stderr=subprocess.STDOUT, shell=True)
        file_list = [p.strip() for p in proc.decode().strip().split('\n')]
        if len(file_list) == 0:
            raise DomainFileError(f'Cannot find pdb file {uid}')
        return(file_list[0])

    def get_structure_path(self, uid: str):
        pathtostructure = self._get_path_to_domain_location(uid)
        fileP = Path(self._grep_file(uid))
        filepath = pathtostructure / f'{fileP.name}'
        if filepath.is_file():
            return filepath
        raise DomainFileError(f'cannot find pdb file for domain {uid}')
