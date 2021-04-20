import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from analysis import ConcatData, cached

RMS2D = """
parm {top}
trajin {traj}
rms first
run
rms2d out {cpp_dir}{prefix}pairwise.dat
run
"""

class Pairwise(ConcatData):
    template = RMS2D
    
    @cached
    def pairwise(self):
        fn = self.cpp_path('pairwise.dat')
        if self._verbose:
            print(f'Loading from {fn}')
        arr = np.loadtxt(fn, skiprows=1)[:, 1:]/10  # skip frames
        return arr

    @cached
    def spaced(self):
        return self.get_spaced(filename='spaced.npy')

    def __init__(self, *args, traj_sel='all', **kwargs):
        # waaaaay too big to analyse without skipping as well
        super(Pairwise, self).__init__(*args, traj_sel='all_skip', **kwargs)

    def get_spaced(self, spacing=0.05, **kwargs):
        nf = int(self.pairwise.shape[0]/self.n_traj)
        n_space = int(spacing*nf)
        nfs = n_space + nf
        n_spaced = self.pairwise.shape[0] + (n_space*(self.n_traj-1))
        arr = np.zeros((n_spaced, n_spaced))
        arr[:] = np.nan

        for i in range(self.n_traj):
            for j in range(self.n_traj):
                pi, ai = i*nf, i*nfs
                pj, aj = j*nf, j*nfs
                arr[ai:ai+nf, aj:aj+nf] = self.pairwise[pi:pi+nf, pj:pj+nf]

        return arr

    def plot_pairwise(self, scale_figsize=3, cmap='magma', vmax=None, 
                      filename='pairwise.tif', **kwargs):
        fs = plt.rcParams.get('figure.figsize')
        figsize = (scale_figsize*fs[0], scale_figsize*fs[1])
        fig = plt.figure(figsize=figsize)

        if vmax is None:
            vmax = np.floor(self.pairwise.max()*10)/10 + 0.1

        pair_x = self.pairwise.shape[0]
        space_x = self.spaced.shape[0]
        n_frames = int(pair_x/self.n_traj)
        halfway = int(np.floor(n_frames/2))
        n_space = int((space_x-pair_x)/(self.n_traj-1))
        ticks = np.arange(self.n_traj)*(n_frames+n_space) + halfway

        ax = sns.heatmap(self.spaced, vmin=0, vmax=vmax, cmap=cmap, square=True,
                         cbar_kws={'label': 'RMSD (nm)'})
        ax.set_xticks(ticks)
        ax.set_xticklabels(self.sim_labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(self.sim_labels)
        plt.tight_layout()

        self.save_fig(filename)
        return ax
