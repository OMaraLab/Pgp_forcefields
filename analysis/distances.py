import os
import glob
import warnings

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib import distances as dist
from analysis import DatasetAnalysis, cached, force_verbose, datafile
from dist_base import Sequence, PgpSheet
from plots import violinplot, kdeplot, axspan, axline
from info import DOMAINS_REV


class Distances(DatasetAnalysis):
    def __init__(self, excel_filename='$HOME/OneDrive*/pgp-forcefields/P-gp.xlsx',
                 **kwargs):
        super(Distances, self).__init__(**kwargs)
        self.excel_filename = glob.glob(os.path.expandvars(excel_filename))[0]
        self.sca = 2.8  # from pymol
        self.ra = 'Residue nA'
        self.rb = 'Residue nB'
        self.da = 'Domain A'
        self.db = 'Domain B'
        self.halves = [' <250 ns', '>=250 ns']

    def read_excel(self, sheet):
        return pd.read_excel(self.excel_filename, sheet_name=sheet)

    @cached
    def xlinks_data(self):
        return self.get_data('Cross_links_clean',
                              self._get_ca_distances_xlinks_nm,
                              filename='xlinks_data.csv')

    @cached
    def epr_data(self):
        return self.get_data('EPR_clean',
                              self._get_ca_distances_epr_nm,
                              filename='epr_data.csv')

    def _get_ca_distances_xlinks_nm(self, xlinks, linkers):
        lnk = {k:v for k, v in zip(linkers.Reagent, 
                                   linkers['Spacer arm / S-S distance (A)'])}
        ss_ranges = np.zeros((len(xlinks.Linker), 2))
        for i, entry in enumerate(xlinks.Linker):
            if not isinstance(entry, str):
                if xlinks.Min_ca.values[i]:
                    ss_ranges[i, 0] = xlinks.Min_ca.values[i] - self.sca*2
                    ss_ranges[i, 1] = xlinks.Max_ca.values[i] - self.sca * 2
                else:
                    ss_ranges[i] = [0, 2.02]  # average S-S bond
            elif entry.strip():
                links = [l.strip() for l in entry.split(',')]
                low, high = links[0], links[-1]
                if low.startswith('M') and not low.endswith('M'):
                    low += 'M'
                ss_ranges[i] = [lnk[low], lnk[high]]

        # get ca-ca distances
        ss_ranges += self.sca*2
        return ss_ranges.T/10

    def _get_ca_distances_epr_nm(self, epr, *args):
        return epr.Min_l/10, epr.Max_l/10

    def get_data(self, sheet_name, ca_func, **kwargs):
        seq = Sequence.from_alignment()
        df = self.read_excel(sheet_name)
        linkers = self.linkers

        for x in 'AB':
            ha, hn = f'H{x}_aa', f'H{x}_n'
            ma, mn = f'M{x}_aa', f'M{x}_n'

            hn0, ha0 = df[hn][0], df[ha][0]
            ha0 = ha0 if (isinstance(ha0, str) and ha0) else None
            mn0, ma0 = df[mn][0], df[ma][0]
            ma0 = ma0 if (isinstance(ma0, str) and ma0) else None
            if np.isnan(hn0) or not ha0:
                df[hn], df[ha] = seq.mouse2human(df[mn].values)
            if np.isnan(mn0) or not ma0:
                df[mn], df[ma] = seq.human2mouse(df[hn].values)

        
        # better names
        df['Entry'] = df.index
        df['Residue nA'] = df['MA_n']
        df['Residue A'] = df['MA_aa'] + df['MA_n'].map(str)
        df['Residue nB'] = df['MB_n']
        df['Residue B'] = df['MB_aa'] + df['MB_n'].map(str)
        df['Residues'] = df['Residue A'] + '-' + df['Residue B']

        if 'Environment' in df.columns:
            df['Res'] = df['Residues'] + ' (' + df['Environment'] + ')'
        else:
            inttemp = df['Min temp (° C)'].astype(int).map(str)
            df['Res'] = df['Residues'] + ' (' + inttemp + '°C)'

        # get domains
        df['Domain A'] = [DOMAINS_REV[x] for x in df.MA_n]
        df['Domain B'] = [DOMAINS_REV[x] for x in df.MB_n]
        df['Domains'] = df['Domain A'] + '-' + df['Domain B']
    
        # get ca-ca distance
        df['low'], df['high'] = ca_func(df, linkers)
        df['Min CA-CA distance (nm)'] = df['low']
        df['Max CA-CA distance (nm)'] = df['high']
        return df
        

    @cached
    def linkers(self):
        return self.read_excel('Linker_lengths')

    @cached
    def mouse_domains(self):
        return self.read_excel('Mouse Domains')
    
    @cached
    def dist_xlinks(self):
        return self.get_distances(data_type='xlinks',
                                  filename='dist_xlinks.csv')
    
    @cached
    def dist_epr(self):
        return self.get_distances(data_type='epr',
                                  filename='dist_epr.csv')

    @cached
    def nbd_nbd(self):
        return self.get_nbd_nbd_distances(filename="nbd_nbd.csv")

    
    def get_nbd_nbd_distances(self, verbose=False, **kwargs):
        results = np.zeros((self.n_traj, self.n_frames))
        for i, u in enumerate(self._data._universes):
            nbd1 = u.atoms[359:597]
            nbd2 = u.atoms[946:1188]
            for j, ts in enumerate(u.trajectory):
                com1 = nbd1.center_of_geometry()
                com2 = nbd2.center_of_geometry()
                results[i, j] = dist.calc_bonds(com1, com2, box=ts.dimensions)
        
        xnbd1 = self.crystal.atoms[359:597].center_of_geometry()
        xnbd2 = self.crystal.atoms[946:1188].center_of_geometry()
        xdist = dist.calc_bonds(xnbd1, xnbd2, box=self.crystal.dimensions)

        df = pd.DataFrame(results.T/10, columns=self._data.sim_labels)
        df["Crystal"] = xdist/10
        df["Time (ns)"] = self.time_ns
        tidy = pd.melt(df, id_vars=["Time (ns)"], var_name="Simulation",
                       value_name="Distance (nm)")
        tidy["Force field"] = [x.split("#")[0].strip() for x in tidy["Simulation"]]
        return tidy
        
    
    def _calculate_distances(self, rid_a, rid_b, verbose=False):
        verbose = verbose or self._verbose
        ca = self.first.select_atoms('name CA')
        aga = ca.select_atoms(*[f'resid {a}' for a in rid_a])
        agb = ca.select_atoms(*[f'resid {b}' for b in rid_b])
        ida, idb = aga.indices, agb.indices

        results = np.zeros((self.n_traj, self.n_frames, len(aga)))
        for i, u in enumerate(self._data._universes):
            if verbose:
                print(f'Analysing {self._data.sim_labels[i]}')
            aga = u.atoms[ida]
            agb = u.atoms[idb]
            for j, ts in enumerate(u.trajectory):
                ca = dist.calc_bonds(aga.positions, agb.positions,
                                     box=ts.dimensions)
                results[i, j] = ca
        
        xca = self.crystal.select_atoms('segid A and name CA')
        xa = xca.select_atoms(*[f'resid {a}' for a in rid_a])
        xb = xca.select_atoms(*[f'resid {b}' for b in rid_b])
        crystal = dist.calc_bonds(xa.positions, xb.positions,
                                  box=self.crystal.dimensions)
        return results, crystal

    def get_distances(self, data_type='xlinks', **kwargs):
        data_df = getattr(self, data_type+'_data')
        rid_a = data_df['MA_n'].values.astype(int)
        rid_b = data_df['MB_n'].values.astype(int)
        results, crystal = self._calculate_distances(rid_a, rid_b)
        
        # n_sim, n_frame, n_res
        arr = np.concatenate(results)  # n_sim * n_frame, n_res

        df = pd.DataFrame(arr/10, columns=data_df.Residues)
        df['Simulation'] = self._data.sim_labels_rep
        df['Force field'] = self._data.ff_labels_rep
        df['Time (ns)'] = np.tile(self.time_ns, self.n_traj)
        df['Half'] = [self.halves[int(x >= df['Time (ns)'].max()/2)]
                      for x in df['Time (ns)']]
        
        # => n_sim * n_frame * n_res
        df = pd.melt(df, id_vars=['Simulation', 'Force field', 'Time (ns)', 'Half'],
                     var_name='Residues', 
                     value_name='Distance (nm)')
        
        # reshape crystal from (n_res) => (n_sim * n_frame * n_res, 1)
        df['Crystal (nm)'] = np.repeat(crystal/10, self.nfsim)

        # add other stuff
        for label in ( 'Residue A', 'Residue B', 'Domain A', 'Domain B', 
                       'Domains', 'Residue pairs', 
                       'Min CA-CA distance (nm)', 'Max CA-CA distance (nm)',
                       'low', 'high',  # convenience
                       'Entry', 'Environment', 'n_peaks', 'Linker', 
                       'Peak range (Å)', 'Min temp (° C)', 'ATP', 
                       'Notes', 'Source', 'DOI', 'Res'):
            if label in data_df.columns:
                df[label] = np.repeat(data_df[label].values, self.nfsim)
                if label == 'Min temp (° C)':
                    df['temp'] = df[label]
                elif label == 'Peak range (Å)':
                    df['range'] = df[label]

        return df
    
    def count_domains(self, df='xlinks'):
        if isinstance(df, str):
            if not df.endswith('_data') or df.startswith('dist_'):
                df += '_data'
            df = getattr(self, df)
        ids, counts = np.unique(df['Domains'], return_counts=True)
        counters = {k:v for k, v in zip(ids, counts)}
        return counters
    
    def res_subdf(self, df='xlinks', a=None, b=None):
        if isinstance(df, str):
            df = getattr(self, 'dist_'+df)
        if a is None:
            a = df[self.ra][0]
        if b is None:
            b = df[self.rb][0]
        
        subdf = df[(df[self.ra]==a) & (df[self.rb]==b)]
        return subdf
    
    def dom_subdf(self, df='xlinks', a=None, b=None):
        if isinstance(df, str):
            df = getattr(self, 'dist_'+df)
        if a is None:
            a = df[self.da][0]
        if b is None:
            b = df[self.db][0]
        
        subdf = df[(df[self.da]==a) & (df[self.db]==b)]
        return subdf
        
    
    def timeplot(self, df='xlinks', a=None, b=None, dom='res'):            
        if dom == 'res':
            func = self.res_subdf
        else:
            func = self.dom_subdf
        subdf = func(df, a, b)
        ax = sns.lineplot(data=subdf,
                          x='Time (ns)', y='Distance (nm)', 
                          hue='Simulation', ci=None,
                          palette=self._data.color_labels_rep,
                          hue_order=self._data.sim_labels_rep)
        plt.title(subdf.Residues.values[0])
        return ax
    
    def kdeplot(self, df='xlinks', a=None, b=None, dom='res'):            
        if dom == 'res':
            func = self.res_subdf
        else:
            func = self.dom_subdf
        subdf = func(df, a, b)

        fig, ax = plt.subplots()
        for _name, _col in zip(self.ff_names, self.ff_colors):
            arr = subdf[subdf['Force field']==_name]['Distance (nm)'].values
            sns.kdeplot(arr, color=_col, label=_name, ax=ax)

        plt.title(subdf.Residues.values[0])
        plt.legend()
        return ax, arr

    def plot_nbd_nbd(self, aspect=2, height=2, filename="nbd_nbd.tif"):
        xtal = self.nbd_nbd[self.nbd_nbd['Force field'] == "Crystal"]["Distance (nm)"]
        df = self.nbd_nbd[self.nbd_nbd['Force field'] != "Crystal"]

        g = sns.FacetGrid(df, col="Force field", hue="Simulation",
                          margin_titles=True, aspect=aspect,
                          height=height, sharey=True, sharex=True,
                          palette=self._data.color_labels,
                          hue_order=self._data.sim_labels)
        for ax in g.axes[0]:
            ax.axhline(xtal.values[0], color="grey", alpha=0.55, ls="--",
                       lw=3)
        g.map_dataframe(sns.lineplot, y="Distance (nm)", x="Time (ns)", ci=None)
        g.axes[0][0].set_ylabel("Distance (nm)")
        for ax in g.axes[0]:
            ax.set_xlabel("Time (ns)")
        g.set_titles(col_template="")
        self.save_fig(filename, g=g)
        return g

    def plot_dist(self, df='xlinks', col='Domains', row='Force field',
                  value='NBD1-NBD2', df_col='Res',
                  col_wrap=4, margin_titles=True, aspect=2,
                  sharey=True, sharex='col', bg_color='grey',
                  bg_alpha=0.2, crystal_color='grey',
                  crystal_ls='--', crystal_lw=3, shade=False,
                  crystal_alpha=0.55, height=1, x='Distance (nm)',
                  suffix=''):

        if isinstance(df, str):
            df_str = df
            df = getattr(self, 'dist_'+df)
        else:
            df_str = ['epr', 'xlinks'][int('Linker' in df.columns)]
        
        
        subdf = df[df[col]==value]
        g = sns.FacetGrid(subdf, col=df_col, row=row, hue='Simulation',
                          margin_titles=margin_titles, aspect=aspect,
                          sharey=sharey, sharex=sharex, dropna=False,
                          palette=self._data.color_labels,
                          hue_order=self._data.sim_labels)

        g.map_dataframe(axspan, color=bg_color, orient='h',
                        alpha=bg_alpha)

        g.map_dataframe(axline, x='Crystal (nm)', color=crystal_color,
                        orient='h', alpha=crystal_alpha,
                        ls=crystal_ls, lw=crystal_lw)

        g.map_dataframe(kdeplot, x=x, shade=shade)

        g = self.reset_margin_titles(g)
        g.set_xlabels(x)
        plt.tight_layout()

        filename = f'kde_{df_str}{suffix}_{value}_{col}.tif'
        self.save_fig(filename, g=g)
        return g

    def plot_violins(self, df='xlinks', col='Domains', row=None,
                     value='NBD1-NBD2', df_col='Res',
                     col_wrap=4, margin_titles=True, aspect=1.5,
                     sharey=True, sharex=False, orient='h',
                     inner='box', hatch=None, bg_color='grey',
                     bg_alpha=0.2, crystal_color='grey',
                     crystal_ls='--', crystal_lw=3,
                     crystal_alpha=0.55, legend=True, scale_hue=True,
                     width=1.5, height=2, xticks=True,
                     gridsize=200, x='Force field', y='Distance (nm)',
                     suffix=''):

        if isinstance(df, str):
            df_str = df
            df = getattr(self, 'dist_'+df)
        else:
            df_str = ['epr', 'xlinks'][int('Linker' in df.columns)]
        
        
        subdf = df[df[col]==value]
        g = sns.FacetGrid(subdf, col=df_col, row=row, col_wrap=col_wrap,
                          margin_titles=margin_titles, aspect=aspect,
                          sharey=sharey, sharex=sharex, dropna=False)

        g.map_dataframe(axspan, x='low', y='high', color=bg_color, 
                        orient=orient, alpha=bg_alpha)

        g.map_dataframe(axline, x='Crystal (nm)', color=crystal_color,
                        orient=orient, alpha=crystal_alpha,
                        ls=crystal_ls, lw=crystal_lw)

        xticks = self.ff_names if xticks else []
        g = self.violinplot(g, x=x, y=y, xticks=xticks, hue='Half', 
                            hue_order=self.halves, palette=self.ff_colors,
                            order=self.ff_names, orient=orient, inner=inner,
                            hatch=hatch, split_palette=True, legend=legend,
                            scale_hue=scale_hue, width=width,
                            gridsize=gridsize)

        g = self.reset_margin_titles(g)

        filename = f'viol_{df_str}{suffix}_{value}_{col}.tif'
        self.save_fig(filename, g=g)
        return g
        
    def plot_violins_dom(self, df='xlinks', row=None,
                     domain='TM8', df_col='Res',
                     col_wrap=4, margin_titles=True, aspect=1.5,
                     sharey=True, sharex=False, orient='h',
                     inner='box', hatch=None, bg_color='grey',
                     bg_alpha=0.2, crystal_color='grey',
                     crystal_ls='--', crystal_lw=3,
                     crystal_alpha=0.55, legend=True, scale_hue=True,
                     width=1.5, height=2, xticks=True,
                     gridsize=200, x='Force field', y='Distance (nm)',
                     suffix=''):

        if isinstance(df, str):
            df_str = df
            df = getattr(self, 'dist_'+df)
        else:
            df_str = ['epr', 'xlinks'][int('Linker' in df.columns)]
        
        
        subdf = df[(df["Domain A"] == domain) | (df["Domain B"] == domain)]
        g = sns.FacetGrid(subdf, col=df_col, row=row, col_wrap=col_wrap,
                          margin_titles=margin_titles, aspect=aspect,
                          sharey=sharey, sharex=sharex, dropna=False)

        g.map_dataframe(axspan, x='low', y='high', color=bg_color, 
                        orient=orient, alpha=bg_alpha)

        g.map_dataframe(axline, x='Crystal (nm)', color=crystal_color,
                        orient=orient, alpha=crystal_alpha,
                        ls=crystal_ls, lw=crystal_lw)

        xticks = self.ff_names if xticks else []
        g = self.violinplot(g, x=x, y=y, xticks=xticks, hue='Half', 
                            hue_order=self.halves, palette=self.ff_colors,
                            order=self.ff_names, orient=orient, inner=inner,
                            hatch=hatch, split_palette=True, legend=legend,
                            scale_hue=scale_hue, width=width,
                            gridsize=gridsize)

        g = self.reset_margin_titles(g)

        filename = f'viol_{df_str}{suffix}_{domain}_{col}.tif'
        self.save_fig(filename, g=g)
        return g