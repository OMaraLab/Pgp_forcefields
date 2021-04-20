import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FormatStrFormatter

from analysis import DatasetAnalysis, cached
from info import DOMAINS


palette = ListedColormap(sns.xkcd_palette(['aqua', 'cerulean', 'light navy',
                                                   'green/yellow', 'blue/green', 'grey', #'pumpkin',
                                                   'raspberry', 'amber']))

class DSSP(DatasetAnalysis):

    def __init__(self, *args, **kwargs):
        super(DSSP, self).__init__(*args, domain='noh2o',
                                   **kwargs)

        prot = self.first.select_atoms('protein').residues
        self.n_residues = len(prot)
        self.resids = prot.resids
        self.ss_colors = {
            'H': 'aqua',
            'G': 'cerulean',
            'I': 'light navy',
            'E': 'green/yellow',
            'B': 'blue/green',
            'C': 'grey',
            'S': 'raspberry',
            'T': 'amber'
        }
        self.index_to_code = {
            0: 'H',
            1: 'G',
            2: 'I',
            3: 'E',
            4: 'B',
            5: 'C',
            6: 'S',
            7: 'T'
        }
        self.code_to_index = dict((v, k) for k, v in self.index_to_code.items())
        self.codes = [v for k, v in sorted(list(self.index_to_code.items()))]
        
        self.sim_labels = self._data.sim_labels[:-5]
        self.ff_labels = self._data.ff_labels[:-5]
        self.tmh_indices = []
        self.resindices = []
        
        self.tmh_labels = []
        self.tmd_labels = []
        for k, (v1, v2) in DOMAINS.items():
            if k.startswith("TM"):
                n = v2 - v1 + 1
                if int(k[2:]) <= 6:
                    tm = "TMD1"
                else:
                    tm = "TMD2"
                self.tmh_labels.extend([k]*n)
                self.tmd_labels.extend([tm]*n)
                self.tmh_indices.extend(np.arange(v1, v2+1))
                rix = self._data.first.select_atoms(f"resnum {v1}:{v2}").residues.resindices
                self.resindices.append(rix)


        # self.tmh_labels = [n for h in self.domains.helices for n in [h.name]*len(h.indices)]
        # self.tmd_labels = ['TMD1']*len(self.domains.TMH1.indices) + ['TMD2']*len(self.domains.TMH2.indices)
        
    @cached
    def helices(self):
        helices = []
        for h in self.resindices:
            helices.append(self.ss_codes[:, :, h])
        return helices

    @cached
    def simple(self):
        simple = []
        for h in self.resindices:
            simple.append(self.ss_simple[:, :, h])
        return simple

    @cached
    def results(self):
        # ss_codes, ss_names, ss_simple, ss_mode, simple_mode
        return self.get_dssp(filename=("ss_codes.npy", "ss_names.npy",
                                       "ss_simple.npy", "ss_mode.npy",
                                       "simple_mode.npy"))

    @cached
    def ss_codes(self):
        a, b = self.try_get_data("ss_codes.npy", force=self._force, verbose=self._verbose)
        if a is None:
            return self.results[0]
        return a

    @cached
    def ss_names(self):
        return self.results[1]

    @cached
    def ss_simple(self):
        return self.results[2]

    @cached
    def ss_mode(self):
        return self.results[3]

    @cached
    def simple_mode(self):
        return self.results[4]

    @cached
    def crystal_ss(self):
        return self.get_crystal(filename=('crystal_ss_codes.npy',
                                          'crystal_ss_simple.npy'))

    @cached
    def crystal_ss_codes(self):
        return self.crystal_ss[0]

    @cached
    def crystal_ss_simple(self):
        return self.crystal_ss[1]

    @cached
    def helical_fraction(self):
        return self.get_helical_fraction(filename="helical_fraction.csv")

    @cached
    def helix_codes(self):
        return self.get_helix_codes(filename='helix_dssp.csv')

    @cached
    def helix_ints(self):
        return self.get_helix_ints(filename='helix_ints.csv')

    def get_crystal(self, verbose=False, **kwargs):
        # this installation of MDAnalysis is a custom branch
        # https://github.com/lilyminium/mdanalysis/tree/dssp
        from MDAnalysis.analysis import secondary_structure as ss
        d = ss.DSSP(self.crystal, select='backbone', verbose=verbose).run()
        return d.ss_codes, d.ss_simple

    def get_dssp(self, **kwargs):
        # this installation of MDAnalysis is a custom branch
        # https://github.com/lilyminium/mdanalysis/tree/dssp
        from MDAnalysis.analysis import secondary_structure as ss

        results = {
            'ss_codes': np.empty((n_u, self.n_frames, self.n_residues), 
                                 dtype=object),
            'ss_names': np.empty((n_u, self.n_frames, self.n_residues),
                            dtype=object),
            'ss_simple': np.empty((n_u, self.n_frames, self.n_residues), dtype=object),
            'ss_mode': np.empty((n_u, self.n_residues), dtype=object),
            'simple_mode': np.empty((n_u, self.n_residues), dtype=object),
        }

        for i, u in enumerate(self._data._universes[:-5]):
            d = ss.DSSP(u, select='backbone').run(verbose=verbose)
            for prop in results:
                attr = getattr(d, prop)
                results[prop] = results[prop].astype(attr.dtype)
                results[prop][i] = attr
        
        return (results["ss_codes"], results["ss_names"], results["ss_simple"],
                results["ss_mode"], results["simple_mode"])

    # def compute_dssp(self, force=False, verbose=False, **kwargs):
    #     from MDAnalysis.analysis import secondary_structure as ss

    #     verbose = verbose or self._verbose
    #     force = force or self._force
    #     n_u = self.n_traj-5


    #     results = {
    #         'ss_codes': np.empty((n_u, self.n_frames, self.n_residues), 
    #                              dtype=object),
    #         'ss_names': np.empty((n_u, self.n_frames, self.n_residues),
    #                         dtype=object),
    #         'ss_simple': np.empty((n_u, self.n_frames, self.n_residues), dtype=object),
    #         'ss_mode': np.empty((n_u, self.n_residues), dtype=object),
    #         'simple_mode': np.empty((n_u, self.n_residues), dtype=object),
    #     }
    #     dtypes = dict.fromkeys(results)
    #     paths = {}
    #     for prop in results:
    #         filename = f'{prop}.npy'
    #         data, path = self.try_get_data(filename, force=force, verbose=verbose)
    #         if data is not None:
    #             results[prop] = self._cache[prop] = data
    #         paths[prop] = path

    #     if all(x in self._cache for x in results):
    #         return results
        
    #     if verbose:
    #         print('Re-running')

    #     for i, u in enumerate(self._data._universes[:-5]):
    #         if verbose:
    #             print(os.path.dirname(u.trajectory.filename))
    #         d = ss.DSSP(u, select='backbone').run(verbose=verbose)
    #         for prop in results:
    #             attr = getattr(d, prop)
    #             results[prop] = results[prop].astype(attr.dtype)
    #             results[prop][i] = attr
        
    #     for prop, arrays in results.items():
    #         self._cache[prop] = arrays
    #         self.save_data(arrays, paths[prop])
        
    #     return results

    def _get_helical_fraction(self, data):
        t1i, t1j = self.domains.TMD1.indices[[0, -1]]
        t2i, t2j = self.domains.TMD2.indices[[0, -1]]

        tm1 = data.T[t1i:t1j].T  #  (n_traj, n_frame, n_res)
        tm1_frac = (tm1=='Helix').sum(axis=-1)/tm1.shape[-1]  #  (n_traj, n_frame)
        
        tm2 = data.T[t2i:t2j].T  # (n_traj, n_frame, n_res)
        tm2_frac = (tm2=='Helix').sum(axis=-1)/tm2.shape[-1]  #  (n_traj, n_frame)

        arr = np.concatenate([tm1_frac, tm2_frac], axis=-1)   # (n_traj, 2*n_frame)
        return arr

    def get_helical_fraction(self, verbose=False, force=False, window=9, **kwargs):

        arr = self._get_helical_fraction(self.ss_simple)
        
        df = pd.DataFrame(arr.T, columns=self.sim_labels)
        df['TMD'] = np.repeat(['TMD1', 'TMD2'], len(self.time_ns))
        df['Time (ns)'] = np.concatenate([self.time_ns, self.time_ns])

        tidy = pd.melt(df, id_vars=['Time (ns)', 'TMD'], 
                       var_name='Simulation',
                       value_name='Helix fraction')

        arrs = []
        front = window // 2
        back = window - front
        for sim in tidy["Simulation"].unique():
            sdf = tidy[tidy["Simulation"] == sim]
            for tmd in sdf["TMD"].unique():
                subdf = sdf[sdf["TMD"] == tmd]
                col = subdf["Helix fraction"]
                n = len(col)
                new = np.zeros(n)
                for i in range(n):
                    new[i] = col.values[max(i-front, 0):min(i+back, n)].mean()
                arrs.append(new)
        
        tidy["Helix fraction (avg)"] = np.concatenate(arrs)


        tidy['Force field'] = self.ff_labels_from_sims(tidy['Simulation'])
        return tidy

        

    def plot_helical_fraction(self, height=3, aspect=1.3, margin_titles=True,
                              crystal=True, linewidth=3, linestyle='--',
                              linecolor='grey', alpha=0.3, average=False):

        def plot_crystal(data=None, tmd1=1, tmd2=1, **kwargs):
            tmd = data['TMD'].unique()
            if tmd == 'TMD1':
                val = tmd1
            else:
                val = tmd2

            plt.axhline(val, linewidth=linewidth, linestyle=linestyle,
                        color=linecolor, alpha=alpha, zorder=-2)
        
        if average:
            self.helical_fraction["Helical fraction"] = self.helical_fraction["Helix fraction (avg)"]

        g = sns.FacetGrid(self.helical_fraction, col='Force field', 
                          row='TMD', hue='Simulation',
                          palette=self._data.color_labels, sharey=True, sharex=True,
                          aspect=aspect, height=height, margin_titles=margin_titles,
                          ylim=(0, 1))

        if crystal:
            tm1, tm2 = self._get_helical_fraction(self.crystal_ss_simple)
            g.map_dataframe(plot_crystal, tmd1=tm1, tmd2=tm2)

        if average:
            y = "Helical fraction"
        else:
            y = "Helix fraction"
        g.map(sns.lineplot, 'Time (ns)', y, ci=None)
        g = g.set_titles(col_template='', row_template="")

        filename = 'helical_fraction.tif'
        if average:
            filename = 'helical_fraction_avg.tif'
        self.save_fig(filename, g=g)
        return g

    def get_helix_codes(self, **kwargs):
        trajs = np.concatenate(self.helices, axis=-1)  #  (n_traj, n_frame, n_res)
        together = np.concatenate(trajs.T, axis=0)  #  (n_res*n_frame, n_traj)

        df = pd.DataFrame(together, columns=self.sim_labels)
        df['Time (ns)'] = np.tile(self.time_ns, trajs.shape[-1])
        df['Residue'] = np.repeat(self.tmh_indices, trajs.shape[1])
        df['Helix'] = np.repeat(self.tmh_labels, trajs.shape[1])
        df['TMD'] = np.repeat(self.tmd_labels, trajs.shape[1])

        tidy = pd.melt(df, id_vars=['Time (ns)', 'TMD', 'Helix', 'Residue'], 
                       var_name='Simulation',
                       value_name='Structure')
        tidy['Force field'] = self.ff_labels_from_sims(tidy['Simulation'])

        return tidy

    def get_helix_ints(self, **kwargs):
        df = self.helix_codes.replace(self.code_to_index)
        return df

    def _plot_helix_colormap_single(self, df, filename, height=3, aspect=1.5, sharey='row', 
                                    margin_titles=True, xticklabels=500, yticklabels=5,
                                    square=False, rotation=0):
        def heatmap(x, y, z, **kwargs):
            data = kwargs.pop('data')
            d = data.pivot(index=y, columns=x, values=z)
            sns.heatmap(d, **kwargs)

        df["Sim"] = [f"{x.split()[0]} {x.split()[-1]}" for x in df["Simulation"]]


        g = sns.FacetGrid(df, row='Helix', col='Sim', sharey=sharey, 
                        sharex=True, margin_titles=margin_titles, height=height,
                        aspect=aspect)
                        
        g.map_dataframe(heatmap, 'Time (ns)', 'Residue', 'Structure', 
                        cbar=False, xticklabels=xticklabels, yticklabels=yticklabels,
                        cmap=self.palette ,square=square)

        for row in g.axes:
            sm = plt.cm.ScalarMappable(cmap=self.palette, norm=self.norm)
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=row, ticks=self.cbar_ticks,
                                boundaries=self.cbar_bounds)
            cbar.ax.set_yticklabels(self.codes)

        g = g.set_titles(col_template='{col_name}')

        n_x = np.arange(0, len(df['Time (ns)']), xticklabels)
        xlabels = ['{:.0f}'.format(x) for x in df['Time (ns)'].values[n_x]]
        g.set_xticklabels(xlabels, xticklabels)
        for ax in g.axes.flat:
            for lab in ax.get_yticklabels():
                lab.set_rotation(rotation)
        self.save_fig(filename, g=g)
        return g
    
    def plot_helix_colormaps(self, height=3, aspect=1.5, sharey='row', margin_titles=True,
                             verbose=False, plot_all=False, plot_separate=True,
                             xticklabels=500, plot_one=False, rotation=0):

        verbose = verbose or self._verbose
        images = []
        ffs = self.helix_ints['Force field'].unique()

        self.code_colors = [self.ss_colors[x] for x in self.codes]
        self.colors = self._data.color_labels[:-5]
        self.palette = ListedColormap(sns.xkcd_palette(self.code_colors))
        self.norm = mpl.colors.Normalize(vmin=0, vmax=len(self.code_colors))
        self.cbar_ticks = np.arange(0, len(self.code_colors))
        self.cbar_bounds = np.arange(-0.5, 0.5+len(self.code_colors), 1)

        for tmd in ['TMD1', 'TMD2']:
            dfi = self.helix_ints[self.helix_ints['TMD']==tmd]

            if plot_separate:
                for name in ffs:
                    if verbose:
                        print(name, tmd)

                    df = dfi[dfi['Force field']==name]
                    ff_name = name.replace(' ', '-').replace('/', '-').lower()
                    filename = 'hel_dssp_{}_{}.tif'.format(tmd, ff_name)
                    g = self._plot_helix_colormap_single(df, filename, height=height,
                                                        aspect=aspect, sharey=sharey,
                                                        margin_titles=margin_titles,
                                                        xticklabels=xticklabels,
                                                        rotation=rotation)
                    images.append(g)
                    del df
                    if plot_one:
                        return images
        
            if plot_all:
                filename = 'hel_dssp_all_{}.tif'.format(tmd)
                g = self._plot_helix_colormap_single(dfi, filename, height=height,
                                                    aspect=aspect, sharey=sharey,
                                                    margin_titles=margin_titles,
                                                    xticklabels=xticklabels,
                                                    rotation=rotation)
                images.append(g)
            del dfi
        
        return images

