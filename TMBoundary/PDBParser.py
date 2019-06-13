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

    def _residues_to_keep(self, regions: list):
        res_to_keep = set()
        for r in regions:
            R = set(range(r[0], r[1]+1))
            res_to_keep = res_to_keep.union(R)
        return res_to_keep

    def get_region(self, out_file: Path, regions: list):
        pdb_out = []
        out_map = {}
        res_to_keep = self._residues_to_keep(regions)
        for ix, r in enumerate(res_to_keep):
            if r < 1 or r > len(self.pdb_list):
                continue
            pdb_out += [line for line in self.pdb_list[r-1]
                        if line[13:16].strip() == 'CA']
            out_map.update({ix+1:r})
        oFile = open(out_file, 'w')
        for l in pdb_out:
            oFile.write(l)
        oFile.close()
        return out_map
