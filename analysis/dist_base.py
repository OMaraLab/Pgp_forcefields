import glob
import os
from utils import analysis_path

class Sequence:
    @classmethod
    def from_alignment(cls, filename=analysis_path('hPgp-mPgp.txt')):
        with open(filename, 'r') as f:
            contents = f.read().split('\n\n')[3:]
        contents = [x for x in contents if x.startswith('Query')]
        seq = cls()
        for section in contents:
            human, mid, mouse = section.split('\n')
            q, h0, hseq, h1 = human.split()
            s, m0, mseq, m1 = mouse.split()
            seq.add(h0, h1, hseq, m0, m1, mseq)
        return seq

    def __init__(self):
        self.h2m = {}
        self.m2h = {}
        self.human = {}
        self.mouse = {}
    
    def add(self, h0, h1, hseq, m0, m1, mseq):
        hidx = list(range(int(h0), int(h1)+1))
        midx = list(range(int(m0), int(m1)+1))
        
        for hx, mx in zip(hseq, mseq):
            if hx == '-':
                mi = midx.pop(0)
                continue
            if mx == '-':
                hi = hidx.pop(0)
                continue
            hi = hidx.pop(0)
            mi = midx.pop(0)
            self.h2m[hi] = mi
            self.m2h[mi] = hi
            self.human[hi] = hx
            self.mouse[mi] = mx
            
    def get_mouse_from_human(self, hidx):
        mi = self.h2m[hidx]
        return (mi, self.mouse[mi])
    
    def get_human_from_mouse(self, midx):
        hi = self.m2h[midx]
        return (hi, self.human[hi])

    def human2mouse(self, indices):
        mi = [self.h2m[i] for i in indices]
        ma = [self.mouse[i] for i in mi]
        return mi, ma

    def mouse2human(self, indices):
        hi = [self.m2h[i] for i in indices]
        ha = [self.human[i] for i in hi]
        return hi, ha
        

class PgpSheet:
    def __init__(self, filename='$HOME/OneDrive*/pgp-forcefields/P-gp.xlsx'):
        filename = os.path.expandvars(filename)
        filename = glob.glob(filename)[0]
        self.sca = 2.8
        self.epr = pd.read_excel(filename, sheet_name='EPR_clean')
        self.domains = pd.read_excel(filename, sheet_name='Mouse Domains')
        self.xlinks = pd.read_excel(filename, sheet_name='Cross_links_clean')
        self.linkers = pd.read_excel(filename, sheet_name='Linker_lengths')