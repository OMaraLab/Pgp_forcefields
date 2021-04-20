import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms, align
from analysis import DatasetAnalysis, cached

class RMSD(DatasetAnalysis):

    def get_rmsds(self, **kwargs):
        rmsds = np.zeros((2, self.n_traj, self.n_frames))

        for i, u in enumerate(self._data):
            a = align.AlignTraj(u, u, in_memory=True).run()
            r = rms.RMSD(u, select='all').run()
            rmsds[0][i] = r.rmsd[:, -1]
            sel = 'index ' + ' '.join(list(map(str, self.tmd)))
            a = align.AlignTraj(u, u, select=sel).run()
            r = rms.RMSD(u.atoms[self.tmd]).run()
            rmsds[1][i] = r.rmsd[:, -1]

        rmsds /= 10  # nanometers

        arr = np.concatenate(rmsds, axis=1)  # (n_traj, 2*n_frames)
        df = pd.DataFrame(arr.T, columns=self._data.sim_labels)
        domains = [r'C$\alpha$', r'TMD C$\alpha$']
        domain_labels = np.repeat(domains, self.n_frames)
        df['Domain'] = domain_labels
        df['Time (ns)'] = np.concatenate([self.time_ns, self.time_ns])
        tidy = df.melt(id_vars=['Time (ns)', 'Domain'],
                       value_vars=self._data.sim_labels,
                       var_name='Simulation',
                       value_name='RMSD (nm)')
        tidy['Force field'] = [x.split('#')[0].strip() for x in tidy['Simulation']]

        return tidy

    @cached
    def rmsds(self):
        return self.get_rmsds(filename="rmsds.csv")

    def plot(self, height=3, aspect=1.5, sharey='row'):
        g = sns.FacetGrid(self.rmsds, row='Domain', col='Force field',
                          hue='Simulation', palette=self._data.color_labels,
                          height=height, aspect=aspect, margin_titles=True,
                          sharey=sharey)
        g.map(sns.lineplot, 'Time (ns)', 'RMSD (nm)', ci=None)
        g.set_titles('{col_name}')
        plt.tight_layout()
        self.save_fig('rmsds.tif', g=g)
        return g

    def plot_moving_average(self, height=3, aspect=1.5, sharey='row', window=5):

        g = sns.FacetGrid(self.rmsds, row='Domain', col='Force field',
                          hue='Simulation', palette=self._data.color_labels,
                          height=height, aspect=aspect, margin_titles=True,
                          sharey=sharey)
        g.map_dataframe(self.plot_average, x='Time (ns)', y='RMSD (nm)', ci=None, window=window)
        g.set_axis_labels("Time (ns)", "RMSD (nm)")
        g.set_titles(col_template="{col_name}", row_template="{row_name}")
        plt.tight_layout()
        self.save_fig(f'rmsds_moving_average_{window:d}.tif', g=g)
        return g

    def plot_average(self, data=None, x=None, y=None, ci=None, window=5, color=None, label=None, **kwargs):
        df = data
        values = df["RMSD (nm)"].values
        avg = []
        lower = window // 2
        for i in range(1, lower+1):
            x = values[:i].mean()
            avg.append(x)
        for i in range(len(values) - window + 1):
            avg.append(values[i:i+window].mean())
        for i in list(range(1, lower+1))[::-1]:
            x = values[-i:].mean()
            avg.append(x)
        sns.lineplot(x=df["Time (ns)"].values, y=avg, color=color, ci=None, label=label, **kwargs)
