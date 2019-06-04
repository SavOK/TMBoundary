import subprocess

class ProteinChainError(Exeption):
    def __init__(self, msg):
        self.msg=msg

    def __repr__(self):
        return self.msg

class ProteinChain:
    def __init__(self, wd:str, pdb, chain, prog:str=None ):
        if prog is None:
            prog = '~rschaeff/bin/generate_pc_pdb.pl'
        self.wd = wd    
        self.pdb = pdb
        self.chain = chain
        self.prog = prog
        

    def _chain_loc (self):
        return  Path(f'{self.wd}/{self.pdb}_{self.chain}.pdb')

    def get_chain_file(self):
        out_path = self._chain_loc()
        if out_path.is_file():
            return out_path
        if not self.wd.is_dir():
            self.wd.mkdir()
        proc_str=f'cd {self.wd} ;{self.prog} {self.pdb} {self.chain}' 
        proc = subprocess.check_output(proc_str,shell=True)
        if not out_path.is_file():
            errStr = f"Cannot get chain {self.pdb} {self.chain}"
            raise ProteinChainError(errStr)
        return out_path



    
