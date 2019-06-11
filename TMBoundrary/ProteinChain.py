from pathlib import Path
import subprocess


class ProteinChainError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        return self.msg


class ProteinChain:
    def __init__(self, wd: Path, pdb: str = None, chain: str = None,
                 prog: str = None):
        if prog is None:
            prog = '~rschaeff/bin/generate_pc_pdb.pl'
        if isinstance(wd, str):
            wd = Path(wd)
        self.wd = wd
        self.pdb = pdb
        self.chain = chain
        self.prog = prog

    def _chain_loc(self):
        return self.wd/f'{self.pdb}_{self.chain}.pdb'

    def get_chain_file(self, pdb: str = None, chain: str = None):
        if pdb is None:
            pdb = self.pdb
        else:
            self.pdb = pdb
        if chain is None:
            chain = self.chain
        else:
            self.chain=chain
        if chain is None or pdb is None:
            raise ProteinChainError('PDB or Chain is not set')
        if self.prog is None:
            raise ProteinChainError('Program path is not set')
        out_path = self._chain_loc()
        if out_path.is_file():
            return out_path
        if not self.wd.is_dir():
            self.wd.mkdir()
        proc_str = f'cd {str(self.wd)} ; {self.prog} {pdb} {chain}'
        print(proc_str)
        print(out_path)
        proc = subprocess.check_output(proc_str, shell=True)
        if not out_path.is_file():
            errStr = f"Cannot get chain {pdb} {chain}"
            raise ProteinChainError(errStr)
        return out_path
