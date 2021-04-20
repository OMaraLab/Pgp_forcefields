
import pandas as pd
import numpy as np
import colorcet as cc
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import (MultipleLocator, 
                               FormatStrFormatter, 
                               AutoMinorLocator) 
from MDAnalysis.analysis import rms

from analysis import ConcatData, cached

colors = ['#272994', '#FFC131', '#F9337B', '#2BC1E2', '#580E33', 
          '#095400', '#EC7F40', '#38B9A3', '#f0bcc7', '#0286D2', 
          '#aaffc3', '#BD38A7', '#6F6F20', '#800000', '#e6beff', 
          '#896444', '#bfef45', '#cc3311', '#f4c99f', '#225555',
          '#CECECE', '#000000']


CLUSTER_TEMPLATE = """
parm {top}
trajin {traj}
rms first
average crdset avg
run
cluster C1 dbscan minpoints 4 epsilon {eps} sievetoframe \\
	rms \\
	sieve 4 \\
	loadpairdist pairdist {kdist_dir}CpptrajPairDist \\
	out {prefix}cl_vs_time.dat \\
	info {prefix}info.dat \\
	repout rep repfmt pdb \\
	avgout avg avgfmt pdb \\
    summary {prefix}summary.dat \\
	summarysplit {prefix}split.dat \\
	splitframe {frame_split}
run
"""

class ClusterAnalysis(ConcatData):
    
    template = CLUSTER_TEMPLATE
    cmap = sns.color_palette(colors)

    @classmethod
    def from_kdist(cls, kdist, eps=2, **kwargs):
        dct = {'eps': eps, 
               'domain': kdist._domain, 
               'time': kdist._time,
               'png': kdist.png,
               'data_dir': kdist.data_dir,
               'image_dir': kdist.image_dir,
               'verbose': kdist._verbose,
               'force': kdist._force,
               'kdist_dir': kdist.cpp_dir}
        dct.update(kwargs)
        return cls(**dct)

    def __init__(self, eps=2, prefix='', data_dir='data', image_dir='images',
                 cpptraj_dir='cpptraj', domain='all', time=500, kdist_dir=None,
                 **kwargs):
        prefix += f'eps{eps}_'
        self._eps = eps
        cpptraj_dir = f'{domain}_{time:d}_cpptraj'
        if kdist_dir is None:
            kdist_dir = cpptraj_dir
        if kdist_dir[-1] != '/':
            kdist_dir += '/'
        self.kdist_dir = kdist_dir

        cpptraj_dir += f'/eps{eps}'

        super(ClusterAnalysis, self).__init__(prefix=prefix, data_dir=data_dir,
                                              image_dir=image_dir, 
                                              cpptraj_dir=cpptraj_dir, 
                                              domain=domain, time=time, **kwargs)
    
    def run(self):
        frames = np.arange(1, self.n_traj) * self.n_frames
        frame_split = ','.join(map(str, frames))
        super().run(eps=self._eps, frame_split=frame_split, 
                    kdist_dir=self.kdist_dir)

    @cached
    def summary(self):
        return self.get_summary(filename='split_df.csv')

    @cached
    def cl_vs_time(self):
        return self.get_clvtime(filename='clvtime_df.csv')

    @cached
    def cluster_rmsds(self):
        return self.get_cluster_rmsds(filename='cluster_rmsds.csv')
    
    def get_summary(self, **kwargs):
        info, split = self.cpp_path('info.dat'), self.cpp_path('split.dat')
        summary = self.cpp_path('summary.dat')

        df = pd.read_csv(summary, header=0, sep='\s+')
        df.rename(columns={'#Cluster': 'Cluster'}, inplace=True)
        df['Cluster'] = df['Cluster']+1
        df['Average distance between points (nm)'] = df['AvgDist']/10
        df['Standard deviation (nm)'] = df['Stdev']/10
        df["Average distance to other clusters (nm)"] = df['AvgCDist']/10
        return df

    def get_cluster_rmsds(self, **kwargs):
        sdf = self.summary
        rmsds = np.zeros((len(sdf.Cluster), len(sdf.Cluster)))
        frames = sdf.Centroid.values
        for i, f1idx in enumerate(frames):
            self.traj.trajectory[f1idx]
            xyz1 = self.traj.atoms.positions.copy()
            for j, f2idx in enumerate(frames[i+1:], i+1):
                self.traj.trajectory[f2idx]
                xyz2 = self.traj.atoms.positions.copy()
                r = rms.rmsd(xyz1, xyz2, superposition=True)
                rmsds[i, j] = rmsds[j, i] = r/10
        df = pd.DataFrame(data=rmsds, index=sdf.Cluster, columns=sdf.Cluster)
        return df


    def get_clvtime(self, **kwargs):
        fn = self.cpp_path('cl_vs_time.dat')
        df = pd.read_csv(fn, sep='\s+', header=0)
        df.rename(columns={'#Frame': 'Frame', 'C1': 'Cluster'}, inplace=True)
        df['Force field'] = self._data.ff_labels_rep
        df['Simulation'] = self._data.sim_labels_rep
        df['Replicate'] = self._data.rep_labels_rep
        df['Time (ns)'] = np.tile(self.time_ns, self.n_traj)
        df['Cluster'] = df['Cluster']+1
        return df

    def plot_cluster_frames(self, cmap=None, filename='cl_frames.tif',
                            figscale=2):
        fs = plt.rcParams.get('figure.figsize')
        figsize = (figscale*fs[0], figscale*fs[1])
        fig = plt.figure(figsize=figsize)
        if cmap is None:
            cmap = sns.color_palette(colors)
        else:
            cmap = getattr(cc, cmap)
        ax = sns.barplot(data=self.summary, x='Cluster', y='Frac',
                         palette=cmap, alpha=0.9)
        self.save_fig(filename)
        return ax

    def plot_cluster_avg_distance(self, cmap=None, filename="cl_distance_to_cl.tif",
                                  figscale=2):
        fs = plt.rcParams.get('figure.figsize')
        figsize = (figscale*fs[0], figscale*fs[1])
        fig = plt.figure(figsize=figsize)                 
        if cmap is None:
            cmap = sns.color_palette(colors)
        else:
            cmap = getattr(cc, cmap)
        

        ax = sns.barplot(data=self.summary, x='Cluster',
                         y='Average distance to other clusters (nm)',
                         alpha=0.9,
                         palette=cmap)
        self.save_fig(filename)
        return ax

    def plot_cluster_distances(self, cmap=None, filename='cl_distances.tif',
                               figscale=2):
        fs = plt.rcParams.get('figure.figsize')
        figsize = (figscale*fs[0], figscale*fs[1])
        fig = plt.figure(figsize=figsize)                 
        if cmap is None:
            cmap = sns.color_palette(colors)
        else:
            cmap = getattr(cc, cmap)
        

        ax = sns.barplot(data=self.summary, x='Cluster',
                         y='Average distance between points (nm)',
                         yerr=self.summary['Standard deviation (nm)'].values,
                         alpha=0.9,
                         palette=cmap, error_kw={'capsize': 5})
        ax.yaxis.set_major_formatter(FormatStrFormatter('% .2f'))
        self.save_fig(filename)
        return ax

    def plot_cluster_rmsds(self, filename='cluster_rmsds.tif', figscale=3,
                           annot=True, linewidths=0.5, fmt='.2f', vmax=None):
        fs = plt.rcParams.get('figure.figsize')
        figsize = (figscale*fs[0], figscale*fs[1])
        fig = plt.figure(figsize=figsize)
        n_ticks = int(self.cluster_rmsds.columns.max())
        g = sns.heatmap(self.cluster_rmsds, cmap='magma',
                        square=True, annot=annot, linewidths=linewidths,
                        fmt=fmt, vmax=vmax,
                        cbar_kws={'label': 'RMSD (nm)'})
        self.save_fig(filename)
        return g

    
    def plot_cluster_map(self, cmap=None, 
                         filename='clustermap.tif', figsize=(12, 12),
                         xticklabels=500,
                         cbar_pos=(0.02, 0.069, 0.05, 0.74)):

        cols = ['Simulation', 'Time (ns)', 'Cluster']
        df = self.cl_vs_time[cols].pivot(*cols).reindex(self.sim_names)
        ticks = sorted(self.cl_vs_time['Cluster'].unique())[1:]
        n_ticks = max(ticks)

        mask = (df == 0).values.astype(bool)
        if cmap is None:
            cmap = ListedColormap(sns.color_palette(colors, n_colors=n_ticks))
        else:
            cmap = getattr(cc, cmap)

        g = sns.clustermap(df, figsize=figsize, cmap=cmap, mask=mask,
                           xticklabels=xticklabels, vmin=1, vmax=n_ticks,
                           col_cluster=False, row_cluster=False,
                           cbar_pos=cbar_pos,
                           cbar_kws={'label': 'Cluster',
                                     'drawedges': True,
                                    #  'shrink': cbar_shrink,
                                     'ticks': ticks})
        self.save_fig(filename)
        return g

    def plot_heatmap(self, cmap=None, filename='heatmap.tif', figscale=3,
                     xticklabels=500):
        fs = plt.rcParams.get('figure.figsize')
        figsize = (figscale*fs[0], figscale*fs[1])
        fig = plt.figure(figsize=figsize)  

        cols = ['Simulation', 'Time (ns)', 'Cluster']
        df = self.cl_vs_time[cols].pivot(*cols).reindex(self.sim_names)
        ticks = sorted(self.cl_vs_time['Cluster'].unique())[1:]
        n_ticks = max(ticks)

        mask = (df == 0).values.astype(bool)
        if cmap is None:
            cmap = sns.color_palette(colors, n_colors=n_ticks)
        else:
            cmap = getattr(cc, cmap)

        df[mask] = np.nan

        g = sns.heatmap(df, cmap=cmap,
                        xticklabels=xticklabels, vmin=1, vmax=n_ticks,
                        cbar_kws={'label': 'Cluster',
                                  'drawedges': True,
                                  'ticks': ticks})
        self.save_fig(filename)
        df[mask] = 0
        return g
