import os

import pandas as pd
import seaborn as sns
import plotly.express as px

from analysis import ConcatData
from utils import cached, data_path

KDIST = """
parm {top}
trajin {traj}
rms first
average crdset avg
run
cluster C1 dbscan kdist {k} rms sieve 4 savepairdist
run
"""


class Kdist(ConcatData):
    template = KDIST
    
    def __init__(self, k=4, cpptraj_dir='cpptraj', kdist_dir='',
                 domain='all', time=500,
                 **kwargs):

        cpptraj_dir = f'{domain}_{time:d}_' + cpptraj_dir
        if kdist_dir:
            cpptraj_dir = cpptraj_dir + '/' + kdist_dir
        super(Kdist, self).__init__(cpptraj_dir=cpptraj_dir, domain=domain,
                                    time=time, **kwargs)
        self.k = k

    @cached
    def df(self):
        filename = self.cpp_path(f'Kdist.{self.k:d}.dat', prefix=False)
        try:
            return pd.read_csv(filename, header=0, sep='\s+')
        except (OSError, FileNotFoundError):
            self.run()
            return pd.read_csv(filename, header=0, sep='\s+')

    def run(self):
        return super().run(k=self.k)

    def plot_kdist(self, filename='dist.tif'):
        # px.scatter(self.df, x='#Point' y='4-dist', title='4-dist plot')
        ax = sns.lineplot(data=self.df, x='#Point', y=f'{self.k}-dist', ci=None)
        self.save_fig(filename)
        return ax

    def plotly_kdist(self, filename='pdist.png'):
        fig = px.scatter(self.df, x='#Point', y=f'{self.k}-dist')
        self.save_fig(filename, plotly=fig)
        return fig



