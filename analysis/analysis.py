import os
import functools
import subprocess

import numpy as np
import pandas as pd
import seaborn as sns
import MDAnalysis as mda
from matplotlib import rcParams
import matplotlib.pyplot as plt

from base import Domains, Dataset, DOMAIN_NDX, DATA_FILES
from plots import violinplot
from utils import (try_mkdir, read_csv, png,
                   force_verbose, datafile, cached,
                   data_path)


rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica Neue']
sns.set(context='paper', style='ticks', font_scale=2, font='sans-serif')



class AnalysisMeta(type):
    """
    I change my mind a lot. This just sets decorators.
    """
    plot_decorators = []
    get_decorators = [force_verbose, datafile]

    def __init__(cls, name, bases, classdict):
        for name, method in classdict.items():
            if name.startswith('plot_'):
                for dec in cls.plot_decorators[::-1]:
                    method = dec(method)
            if name.startswith('get_'):
                for dec in cls.get_decorators[::-1]:
                    method = dec(method)
            setattr(cls, name, method)

        type.__init__(type, name, bases, classdict)

class AnalysisBase(metaclass=AnalysisMeta):
    def __init__(self, data_dir='data', image_dir='images', 
                 prefix='', verbose=True, force=False,
                 png=False):
        self._cache = {}
        self._verbose = verbose
        self._force = force
        self.png = png
        self.prefix = prefix

        # set up directories
        self.dir = os.getcwd()
        self.image_dir = try_mkdir(self.dir, image_dir)
        self.data_dir = try_mkdir(self.dir, data_dir)
        

    def save_fig(self, filename, prefix=True, g=None, plotly=None, dpi=300):
        filename = filename.replace(' ', '-')
        if self.png or plotly:
            filename = filename[:-3] + 'png'
            dpi = 120
        if prefix:
            filename = self.prefix + filename
        path = os.path.join(self.image_dir, filename)
        if g is not None:
            g.savefig(path, dpi=dpi)
        elif plotly is not None:
            plotly.write_image(path)
        else:
            plt.savefig(path, dpi=dpi)
        print('Saved to', path)
    
    def try_get_data(self, filename, prefix=True,
                     force=False, verbose=False):
        if prefix:
            filename = self.prefix + filename
        path = os.path.join(self.data_dir, filename)
        suffix = path.split('.')[-1]
        
        if not force:
            # what kind of data
            if suffix == 'csv':
                loader = read_csv
            elif suffix in ('dat', 'txt'):
                loader = np.loadtxt
            elif suffix in ('npy', 'npz'):
                loader = np.load
            else:
                raise ValueError(f"Can't find loader for {suffix} file")
            
            try:
                data = loader(path)
            except:
                if verbose:
                    print(f'Could not load data from {path}: rerun.')
            else:
                try:
                    _data = dict(**data)
                    data.close()
                    data = _data
                except:
                    pass
                if verbose:
                    print(f'Loaded from {path}.')
                    return data, path
        
        return None, path

    def save_data(self, data, path, comments=None):
        suffix = path.split('.')[-1]
        if suffix == 'csv':
            data.to_csv(path)
        elif suffix in ('dat', 'txt'):
            np.savetxt(path, data, comments=comments)
        elif suffix == 'npy':
            np.save(path, data)
        elif suffix == 'npz':
            np.savez_compressed(path, **data)
        else:
            raise ValueError(f"Can't find saver for {suffix} file")
        print('Saved to', path)

    def ff_labels_from_sims(self, sim_labels):
        return [x.split('#')[0].strip() for x in sim_labels]

    # def plot_average(self, x, y, data=None, window=5, **kwargs):
    #     arr = data[y].values
    #     means = [arr[i-window:i].mean() for i in range(window, len(arr))]
    #     df = pd.DataFrame({y: means, x: data[x][window:]})
    #     sns.lineplot(data=df, x=x, y=y, **kwargs)

    def lineplot(self, x, y, data=None, window=None, **kwargs):
        sns.lineplot(data=data, x=x, y=y, **kwargs)

    def violinplot(self, g: sns.FacetGrid, x: str=None, y: str=None,
                   orient: str='h', xticks: list=[], **kwargs):
        _y = y
        if orient == 'h':
            y = x
            x = _y
            set_ylabels = g.set_xlabels
            set_xticks = g.set_yticklabels
        else:
            set_ylabels = g.set_ylabels
            set_xticks = g.set_xticklabels

        g.map_dataframe(violinplot, x=x, y=y, **kwargs)
        g = set_ylabels(_y)
        g = set_xticks(xticks)
        return g

    def reset_margin_titles(self, g, col_template='{col_name}',
                            row_template=''):
        for ax in g.axes.flat:
            plt.setp(ax.texts, text="")
        g.set_titles(col_template=col_template, row_template=row_template)
        return g

class DatasetAnalysis(AnalysisBase):
    def __init__(self, time=500, ndx=DOMAIN_NDX, prefix='',
                 domain='all', mod={}, **kwargs):
        super(DatasetAnalysis, self).__init__(prefix=prefix, **kwargs)
        
        # data details
        self._time = time
        self._data = Dataset(domain=domain, time=time, mod=mod)
        self.n_frames = self._data.n_frames
        self.time_ns = np.arange(self.n_frames) * 0.2
        self.n_atoms = self._data.n_atoms
        self.n_traj = self._data.n_traj
        self.first = self._data.first
        self.nfsim = self.n_frames * self.n_traj

        # colours and labels
        self.ff_colors = [ff.color for ff in self._data.forcefields]
        self.ff_names = [ff.name for ff in self._data.forcefields]
        self.ff_names_multiline = [n.replace(' ', '\n') for n in self.ff_names]
        self.sim_colors = self._data.color_labels
        self.sim_names = self.sim_labels = self._data.sim_labels
        self.sim_names_multiline = [n.replace(' ', '\n') for n in self.sim_names]

        self.short2long = {'ff': 'Force field', 'sim': 'Simulation'}
        self.long2short = {'Simulation': 'sim', 'Force field': 'ff'}

        # easy indices
        self.domains = Domains(ndx)
        hel = self.domains.helices
        tm1_start, tm1_end = hel[0].indices[0], hel[5].indices[-1]
        tm2_start, tm2_end = hel[6].indices[0], hel[-1].indices[-1]

        self.tmd1 = np.arange(tm1_start, tm1_end+1)
        self.tmh1 = self.domains.TMD1.indices
        self.tmd2 = np.arange(tm2_start, tm2_end+1)
        self.tmh2 = self.domains.TMD2.indices
        self.tmh = self.domains.TMD.indices
        self.tmd = np.concatenate([self.tmd1, self.tmd2], axis=0)
        self.nbd1 = self.domains.NBD1.indices
        self.nbd2 = self.domains.NBD2.indices
        self.nbd = self.domains.NBD.indices

    @cached
    def crystal(self):
        return mda.Universe(data_path('4m1m.pdb.gz'))

class DomainAnalysis(DatasetAnalysis):
    def __init__(self, domain='all', time=500, prefix='', **kwargs):

        super(DomainAnalysis, self).__init__(prefix=prefix+f'{domain}_{time:d}_',
                                             domain=domain,
                                             time=time, **kwargs)
        self._domain = domain

class ConcatData(DomainAnalysis):
    template = ""

    def __init__(self, domain='all', time=500, cpptraj_dir='cpptraj', 
                 traj_sel='all', **kwargs):
        super(ConcatData, self).__init__(time=time, domain=domain, **kwargs)
        self.cpp_dir = try_mkdir(self.data_dir, cpptraj_dir)
        self.cpp_in = self.cpp_path('cpptraj.in', run_on_missing=False)
        self.cpp_out = self.cpp_path('cpptraj.out', run_on_missing=False)
        self._top = data_path(DATA_FILES[domain]['pdb'])
        self._traj = data_path(DATA_FILES[domain][traj_sel].format(time=time))

    @cached
    def traj(self):
        return mda.Universe(self._top, self._traj)

    def cpp_path(self, *paths, run_on_missing=True, prefix=True):
        p = self.prefix if prefix else ''
        paths = [p+x for x in paths]
        path = os.path.join(self.cpp_dir, *paths)
        if not os.path.exists(path) and run_on_missing:
            self.run()
        return path

    def run(self, **kwargs):
        with open(self.cpp_in, 'w') as f:
            f.write(self.template.format(top=self._top, traj=self._traj,
                                         prefix=self.prefix, 
                                         data_dir=self.data_dir+'/',
                                         cpp_dir=self.cpp_dir+'/',
                                         **kwargs))
        if self._verbose:
            print('Running cpptraj')
        f = open(self.cpp_out, 'w')
        proc = subprocess.run(['cpptraj', '-i', self.cpp_in], stdout=f,
                              cwd=self.cpp_dir)
        f.close()