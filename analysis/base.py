import MDAnalysis as mda
import os
import numpy as np
import seaborn as sns
from utils import get_path, data_path, cached

# def get_ndx(file):
#     selections = {}
#     with open(file, 'r') as f:
#         contents = f.read()
#     _split = contents.split('[')[1:]
#     _split = [x.split(']') for x in _split]
#     for section, numbers in _split:
#         name = section.strip()
#         numbers = numbers.replace('\n', ' ')
#         selections[name] = np.array(list(map(lambda x: int(x)-1, numbers.split())))
#     return selections

DOMAIN_NDX = data_path('ca_domains.ndx')

NAMES = {
    'amber': 'AMBER 99SB-ILDN',
    'charmm': 'CHARMM 36',
    'opls': 'OPLS-AA/L',
    'gromos': 'GROMOS 54A7',
    'martini': 'MARTINI'
}

COLORS = {
    'amber': ['#e5d973', '#ccb24c', '#b28c26', '#996600'],  # '#ffff99' too light
    'gromos': ['#ff6666', '#d94c4c', '#b23333', '#8c1919', '#660000'],
    'charmm': ['#86dbdb', '#73d9f2', '#4cb2e5', '#268cd9', '#1844a1'],
    'opls': ['#cc99ff', '#a673d9', '#7f4cb2', '#59268c', '#330066'],
    'martini': ['#ccff99', '#99d973', '#66b24c', '#328c26', '#006600']
}

# patterns = {
#     # 'PDB': 'data/{ff}/rep{rep:d}/ca20ns.pdb',
#     'PDB': 'data/ca20ns_amber_1.pdb', # MDA doesn't load this anyway
#     'TRAJ200': 'data/{ff}/rep{rep:d}/ca200ns_200ps.xtc',
#     'TRAJ500': 'data/{ff}/rep{rep:d}/ca500ns_200ps.xtc',
#     'TRAJ200ALL': '',
#     'TRAJ500ALL': '',
#     'NOH2O_500': 'data/{ff}/rep{rep:d}/noh2o_500ns_200ps.xtc',
#     'NOH2O_PDB': 'data/{ff}/rep{rep:d}/noh2o_20ns.pdb'
# }


DATA_FILES = {
    'all': {
        'pdb': 'ca20ns_amber_1.pdb',
        'rep': '{ff}/rep{rep}/ca{time}ns_200ps.xtc',
        'all': 'all_ca{time}ns_200ps.xtc',
        'all_skip': 'all_ca{time}ns_1ns.xtc',
    },
    'tmh': {
        'pdb': 'tmh_ca20ns_amber_1.pdb',
        'rep': '{ff}/rep{rep}/tmh_{time}ns_200ps.xtc',
        'all': 'tmh_{time}ns_200ps.xtc',
        'all_skip': 'tmh_{time}ns_1ns.xtc',
    },
    'tmd': {
        'pdb': 'tmd_ca20ns_amber_1.pdb',
        'rep': '{ff}/rep{rep}/tmd_{time}ns_200ps.xtc',
        'all': 'tmd_ca{time}ns_200ps.xtc',
        'all_skip': 'tmd_ca{time}ns_1ns.xtc',
    },
    'noh2o': {
        'pdb': '{ff}/rep{rep}/noh2o_20ns.pdb',
        'rep': '{ff}/rep{rep}/noh2o_{time}ns_200ps.xtc',
    },
    'no_o3_all': {
        'pdb': 'ca20ns_amber_1.pdb',
        'rep': '{ff}/rep{rep}/ca{time}ns_200ps.xtc',
        'all': 'no_o3_ca{time}ns_200ps.xtc',
    },
    'no_o3_tmh': {
        'pdb': 'tmh_ca20ns_amber_1.pdb',
        'rep': '{ff}/rep{rep}/tmh_{time}ns_200ps.xtc',
        'all': 'no_o3_tmh_{time}ns_200ps.xtc',
    },
    'no_o3_tmd': {
        'pdb': 'tmd_ca20ns_amber_1.pdb',
        'rep': '{ff}/rep{rep}/tmd_{time}ns_200ps.xtc',
        'all': 'no_o3_tmd_{time}ns_200ps.xtc',
    }
}

DATA_CONFIG = {
    'amber': [1, 2, 3],
    'charmm': [1, 2, 3],
    'opls': [1, 2, 3],
    'gromos': [1, 2, 3],
    'martini': [1, 2, 3, 4, 5],
}

COLORS = {
    'amber': {1: '#e5d973', 2: '#bf9f39', 3: '#996600'},
    'charmm': {1: '#86dbdb', 2: '#4cb2e5', 3: '#1844a1'},
    'gromos': {1: '#ff6666', 2: '#b23333', 3: '#660000'},
    'opls': {1: '#cc99ff', 2: '#7f4cb2', 3: '#330066'},
    'martini': {1: '#ccff99', 2: '#99d973', 3: '#66b24c', 
                4: '#328c26', 5: '#006600'}
}

class ForceField:
    def __repr__(self):
        return f'ForceField(name={self.name}, n_reps={self.n_reps})'

    @property
    def _universes(self):
        return self.reps

    @cached
    def reps(self):
        return [self.first] + [mda.Universe(*fn) for fn in self.filenames[1:]]

    def __len__(self):
        return self.n_reps

    def __init__(self, ff, domain='all', time=500, reps=[1, 2, 3]):
        self._cache = {}
        colors = COLORS[ff]
        self.name = NAMES[ff]
        self._ff = ff
        if not reps:
            self.colors = []
        else:
            self.colors = [COLORS[ff][i] for i in reps]
        
        if len(COLORS[ff]) == 3:
            self.color = COLORS[ff][2]
        else:
            self.color = COLORS[ff][3]
        self.n_reps = len(reps)
        self.filenames = []

        traj = DATA_FILES[domain]['rep']
        top = DATA_FILES[domain]['pdb']
        
        for i in reps:
            path = data_path(traj.format(ff=ff, rep=i, time=time))
            top = data_path(top.format(ff=ff, rep=i, time=time))
            # self.reps.append(mda.Universe(top, path))
            self.filenames.append((top, path))

        self.color_labels = self.colors
        self.rep_labels = np.arange(1, self.n_reps+1)
        self.sim_labels = self.sim_names = [f'{self.name} #{i}' for i in self.rep_labels]
        self.ff_labels = [self.name]*self.n_reps

        if self.n_reps:
            self.first = mda.Universe(*self.filenames[0])
            self._n_frames = len(self.first.trajectory)
        else:
            self.first = None
            self._n_frames = 0
        self.color_labels_rep = np.repeat(self.colors, self._n_frames)
        self.rep_labels_rep = np.repeat(self.rep_labels, self._n_frames)
        self.sim_labels_rep = np.repeat(self.sim_labels, self._n_frames)
        self.ff_labels_rep = np.repeat(self.ff_labels, self._n_frames)



class Dataset:    
    _forcefields = ('amber', 'charmm', 'opls', 'gromos', 'martini')

    def __iter__(self):
        for u in self._universes:
            yield u
    
    def __getitem__(self, key):
        return self._universes[key]

    @cached
    def _universes(self):
        return [y for x in self.forcefields for y in x.reps]
    
    def __init__(self, domain='all', time=500, mod={}):
        self._cache = {}
        self.forcefields = []

        CONFIG = DATA_CONFIG.copy()
        CONFIG.update(mod)

        for ff in self._forcefields:
            self.forcefields.append(ForceField(ff, domain=domain, time=time,
                                               reps=CONFIG[ff]))
        self.first = self.forcefields[0].first
        self.n_reps = [len(x) for x in self.forcefields]
        self.n_traj = sum(self.n_reps)

        self.n_frames = len(self.first.trajectory)
        self.n_atoms = len(self.first.atoms)
        self.colors = [ff.color for ff in self.forcefields]

        # LABELS
        self.color_labels = [y for x in self.forcefields for y in x.color_labels]
        self.sim_labels = [y for x in self.forcefields for y in x.sim_labels]
        self.rep_labels = [y for x in self.forcefields for y in x.rep_labels]
        self.ff_labels = [y for x in self.forcefields for y in x.ff_labels]
        self._ff_labels = [y for x in self.forcefields for y in [x._ff]*x.n_reps]

        self.color_labels_rep = [y for x in self.forcefields for y in x.color_labels_rep]
        self.sim_labels_rep = [y for x in self.forcefields for y in x.sim_labels_rep]
        self.rep_labels_rep = [y for x in self.forcefields for y in x.rep_labels_rep]
        self.ff_labels_rep = [y for x in self.forcefields for y in x.ff_labels_rep]




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

