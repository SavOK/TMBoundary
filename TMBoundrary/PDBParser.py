from pathlib import Path


class PDBError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        return self.msg


class PDBParser:
    def __init__(self, input_file: Path):
        self.input_file = input_file
        self.pdb_list = self._read_pdb()

    def _read_pdb(self):
        with open(self.input_file) as oFile:
            pdb_list = []
            res_list = []
            for ix, line in enumerate(oFile):
                if not line.startswith('ATOM'):
                    continue
                if ix == 0:
                    chain_prev = line[21]
                    res_prev = line[22:28].strip()
                    res_list.append(line)
                    continue
                chain_curr = line[21]
                res_curr = line[22:28].strip()
                if chain_curr != chain_prev:
                    raise PDBError(
                        f'{self.input_file.stem} Multi chain protein/domain')
                if res_curr != res_prev:
                    pdb_list.append(res_list)
                    res_list = []
                res_list.append(line)
                res_prev = res_curr
                chain_prev = chain_curr
            pdb_list.append(res_list)
        return pdb_list

    def _cut_region(self, start: int, end: int):
        out_list = []
        if end > len(self.pdb_list):
            end = len(self.pdb_list)
        if start < 1:
            start = 1
        for res in self.pdb_list[start-1:end]:
            out_list += [line for line in res if line[13:16].strip() == 'CA']
        return out_list

    def get_region(self, out_file: Path, regions: list):
        pdb_out = []
        for r in regions:
            pdb_out += self._cut_region(r[0], r[1])
        oFile = open(out_file, 'w')
        for l in pdb_out:
            oFile.write(l)
        oFile.close()
