import re
from subprocess import PIPE, Popen
from pathlib import Path


class TMalign:
    def __init__(self, prog: str = None):
        if prog is None:
            prog = '/usr7/TMalign/TMalign'
        self.prog = Path(prog)

    def _get_align_reg(self, ali):
        reg = []
        start = end = 0
        regF = need_add = False
        for ix, c in enumerate(ali):
            if c == ':' or c == '.':
                if not regF:
                    need_add = True
                    regF = True
                    start = ix+1
                    end = start
                if regF:
                    end += 1
            else:
                if regF:
                    reg.append((start, end))
                    regF = False
                    need_add = False
                else:
                    continue
        if need_add:
            reg.append((start, end-1))
        return reg

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
        out_dict['ali_reg'] = self._get_align_reg(output[-4])
        out_dict['tm_align'] = output[-4].strip('\n')
        return out_dict

    def run_align(self,query: Path, hit: Path, cut: int = None):
        if cut is None:
            cut = 5
        args = [str(self.prog), str(query), str(hit), '-d', str(cut)]
        print(args)
        with Popen(args=args, stdout=PIPE) as proc:
            outputlines = [l.strip('\n')
                           for l in proc.stdout.read().decode().split('\n')]
            out_dict = self._parse_TMalign(outputlines)
        return out_dict
