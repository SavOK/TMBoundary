import re
from subprocess import PIPE, Popen
from pathlib import Path


class TMalign:
    def __init__(self, prog: str = None):
        if prog is None:
            prog = '/usr7/TMalign/TMalign'
        self.prog = Path(prog)

    def _get_align_region(self, seq: str, ali: str):
        region = []
        curr_start = 0
        curr_end = 0
        # filter to remove gaps
        L = [x for x in zip(seq, ali)if x[0] != '-']
        # filter to remove misalign
        en = [(ix+1, x) for ix, x in enumerate(L) if x[1] != ' ']
        curr_start = en[0][0]
        curr_end = en[0][0]
        for ix, x in en[1:]:
            if ix - curr_end == 1:  # continue
                curr_end = ix
            else:  # break
                region.append([curr_start, curr_end])
                curr_start = ix
                curr_end = ix
        region.append([curr_start, curr_end])
        return region

    def _parse_TMalign(self, output: list):
        out_dict = {}
        TMscore1_ptr = re.compile(
            r'TM-score.+?(?P<score>\d+[.]\d+).+?(Chain_1)')
        TMscore2_ptr = re.compile(
            r'TM-score.+?(?P<score>\d+[.]\d+).+?(Chain_2)')
        TMscoreN_ptr = re.compile(
            r'TM-score.+?(?P<score>\d+[.]\d+).+?(scaled by user)')
        for line in output:
            match = TMscore1_ptr.match(line)
            if match:
                out_dict['tm1'] = float(match.group('score'))
            match = TMscore2_ptr.match(line)
            if match:
                out_dict['tm2'] = float(match.group('score'))
            match = TMscoreN_ptr.match(line)
            if match:
                out_dict['tmN'] = float(match.group('score'))
        out_dict['query_seq'] = output[-5].strip('\n')
        out_dict['hit_seq'] = output[-3].strip('\n')
        out_dict['query_reg'] = self._get_align_region(
            out_dict['query_seq'], output[-4])
        out_dict['hit_reg'] = self._get_align_region(
            out_dict['hit_seq'], output[-4])
        out_dict['tm_align'] = output[-4].strip('\n')
        return out_dict

    def run_align(self, query: Path, hit: Path, cut: int = None):
        if cut is None:
            cut = 5
        args = [str(self.prog), str(query), str(hit), '-d', str(cut)]
#        print(args)
        with Popen(args=args, stdout=PIPE, stderr=PIPE) as proc:
            outputlines = [l.strip('\n')
                           for l in proc.stdout.read().decode().split('\n')]
        if proc.returncode == 0:
            out_dict = self._parse_TMalign(outputlines)
        else:
            print(f"NOTE: CANNOT align {query} {hit}")
            return None
        return out_dict
