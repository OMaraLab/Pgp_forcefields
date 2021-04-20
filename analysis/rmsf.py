import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from analysis import DatasetAnalysis, cached
from info import DOMAINS


class RMSF(DatasetAnalysis):

    def get_rmsfs(self, force=False, verbose=False, **kwargs):

        if verbose:
            print('Re-running')
        rmsfs = np.zeros((self.n_traj, self.n_atoms))

        for i, u in enumerate(self._data):
            if verbose:
                print(f'Running for {u.trajectory.filename}')
            a = align.AlignTraj(u, u, in_memory=True).run()
            xyz = u.trajectory.timeseries().mean(axis=1)
            ref = mda.Merge(u.atoms).load_new(xyz)

            a = align.AlignTraj(u, ref).run()

            r = rms.RMSF(u.atoms).run()
            rmsfs[i] = r.rmsf

        rmsfs /= 10  # nanometers

        df = pd.DataFrame(rmsfs.T, columns=self._data.sim_labels)
        resids = self.first.residues.resids
        df['Residue'] = resids
        tidy = df.melt(id_vars=['Residue'],
                       value_vars=self._data.sim_labels,
                       var_name='Simulation',
                       value_name='RMSF (nm)')
        tidy['Force field'] = [x.split('#')[0].strip() for x in tidy['Simulation']]
        return tidy

    @cached
    def rmsfs(self):
        return self.get_rmsfs(filename="rmsfs.csv")

    @staticmethod
    def plot_domains(alpha=0.1, nmax=1, color='grey', text_color='black',
                     no_tm=False, ymax=None, rel_height=None,
                     fontsize='small', fontdict=None, rotation=90, label=None, **kwargs):
        if "1" in label:
            for name, (xmin, xmax) in DOMAINS.items():
                if name.startswith("TM") or name.startswith("NBD"):
                    plt.axvspan(xmin, xmax, alpha=alpha, color=color,
                                zorder=-2, linewidth=0)
                    x = xmin+((xmax-xmin)/2)
                    if no_tm:
                        name = name.replace('TM', '')
                    plt.text(x, ymax*rel_height, name, zorder=0, ha='center', va='bottom',
                                fontsize=fontsize, rotation=rotation, fontdict=fontdict,
                                color=text_color)


    def plot(self, height=3, aspect=2, sharey='row',
             color='grey', alpha=0.1, nmax=1, rel_height=0.8, rotation=90,
             no_tm=False, fontsize='small', fontdict=None,
             text_color='black'):

        ymax = self.rmsfs['RMSF (nm)'].values.max()*nmax

        g = sns.FacetGrid(self.rmsfs, row='Force field',
                          hue='Simulation', palette=self._data.color_labels,
                          height=height, aspect=aspect, margin_titles=True,
                          sharey=sharey, ylim=(0, ymax), dropna=False)
        
        g.map(self.plot_domains, color=color, alpha=alpha, nmax=nmax,
              rel_height=rel_height, rotation=rotation, no_tm=no_tm,
              fontsize=fontsize, fontdict=fontdict, text_color=text_color,
              ymax=ymax)
        g.map(plt.plot, 'Residue', 'RMSF (nm)')
        g.set_titles(row_template='', col_template="")
        plt.tight_layout()
        self.save_fig('rmsfs.tif', g=g)
        return g

    def plot_tmd(self, height=3, aspect=1.5, sharey=True, hue='Simulation',
                 legend=True):

        if hue == 'Simulation':
            colors = self._data.color_labels
        else:
            colors = self._data.colors

        g = sns.FacetGrid(self.rmsfs, col='Helix', col_wrap=6,
                          hue=hue, palette=colors,
                          height=height, aspect=aspect, margin_titles=True,
                          sharey=sharey, dropna=True, legend_out=True,
                          sharex=False)
        
        g.map(sns.lineplot, 'Residue', 'RMSF (nm)', ci='sd')
        if legend:
            g.add_legend()
        g.set_titles('{col_name}')
        self.save_fig(f'rmsf_tmds_{hue.replace(" ", "")}.tif', g=g)
        return g
        

    def plot_nbd(self, height=3, aspect=1.5, sharey=True, hue='Simulation',
                 legend=True):

        if hue == 'Simulation':
            colors = self._data.color_labels
        else:
            colors = self._data.colors

        df = self.rmsfs[(self.rmsfs["Domain"] == "NBD1") | (self.rmsfs["Domain"] == "NBD2")]
        g = sns.FacetGrid(df, col='Domain', col_wrap=6,
                          hue=hue, palette=colors,
                          height=height, aspect=aspect, margin_titles=True,
                          sharey=sharey, dropna=True, legend_out=True,
                          sharex=False)
        
        g.map(sns.lineplot, 'Residue', 'RMSF (nm)', ci='sd')
        if legend:
            g.add_legend()
        g.set_titles('{col_name}')
        self.save_fig(f'rmsf_nbds_{hue.replace(" ", "")}.tif', g=g)
        return g
        