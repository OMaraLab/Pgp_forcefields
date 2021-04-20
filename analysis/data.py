import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns

from collections import OrderedDict

from matplotlib import rcParams

from base import *

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica Neue']
sns.set(font_scale=2, font='sans-serif')
sns.set_style("ticks")

DOMAIN_NDX = get_path('data/ca_domains.ndx')

CRYSTAL = get_path('data/4m1m.pdb.gz')

class Selection:
    def __init__(self, name, indices):
        name = name.replace('_', '-')
        self.name = name
        self.indices = np.asarray(indices)

class Selections:
    def __init__(self):
        self.selections = {}
        self.keys = set()
        
    def load_file(self, file):
        with open(file, 'r') as f:
            contents = f.read()
        _split = contents.split('[')[1:]
        _split = [x.split(']') for x in _split]
        for section, numbers in _split:
            name = section.strip().replace('-', '_')
            numbers = numbers.replace('\n', ' ').split()
            self.selections[name] = Selection(name, list(map(lambda x: int(x)-1, numbers)))
            setattr(self, name, self.selections[name])

class Domains(Selections):
    def __init__(self, file):
        super().__init__()
        self.load_file(file)
        self.tmds = [self.TMD, self.TMD1, self.TMD2]
        self.nbds = [self.NBD, self.NBD1, self.NBD2]
        self.helix_names = ['TM{}'.format(x) for x in range(1, 13)]
        self.helices = [self.selections[n] for n in self.helix_names]

    def indices(self, k):
        return self.selections[k].indices

