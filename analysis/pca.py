import glob
import os

import scipy
import numpy as np
import pandas as pd
import seaborn as sns
from MDAnalysis.analysis import pca
import matplotlib.pyplot as plt

from utils import analysis_path
from analysis import ConcatData, cached

CPCA = """
parm {top}
trajin {traj}
rms first
average crdset avg
run
rms ref avg
matrix covar name all-covar
createcrd CRD1
run
runanalysis diagmatrix all-covar out {prefix}evecs.dat name myEvecs nmwiz nmwizvecs 3 nmwizfile vis.nmd
runanalysis modes eigenval name myEvecs out {prefix}evalfrac.dat
runanalysis modes name myEvecs trajout mode_1.xtc pcmin -100 pcmax 100 tmode 1 trajoutfmt xtc
runanalysis modes name myEvecs trajout mode_2.xtc pcmin -100 pcmax 100 tmode 2 trajoutfmt xtc
runanalysis modes name myEvecs trajout mode_3.xtc pcmin -100 pcmax 100 tmode 3 trajoutfmt xtc
runanalysis modes name myEvecs trajout mode_4.xtc pcmin -100 pcmax 100 tmode 4 trajoutfmt xtc
runanalysis modes name myEvecs trajout mode_5.xtc pcmin -100 pcmax 100 tmode 5 trajoutfmt xtc

crdaction CRD1 projection modes myEvecs out proj.dat beg 1 end {n}

run
"""

class PCA(ConcatData):

    template = CPCA

    def __init__(self, *args, n_components=100, pca_path='pca2', cpptraj_dir='cpptraj',
                 domain='all', time=500, **kwargs):
        super(PCA, self).__init__(*args, domain=domain, time=time,
                                  cpptraj_dir=f'{domain}_{time:d}_'+cpptraj_dir,
                                  **kwargs)
        self.n_components = n_components
        self._h0 = '0-500ns'
        self._h1 = '0-250ns'
        self._h2 = '250-500ns'
        self.halves = [' <250 ns', '>=250 ns']

    def run(self):
        return super().run(n=self.n_components)
        
    def __getattr__(self, key):
        try:
            return self._cache[key]
        except KeyError:
            pass

        for prefix, func in [('rmsip_pc', self.get_rmsip),
                             ('rmsip_norm_pc', self.get_rmsip_norm),
                             ]:
            if key.startswith(prefix):
                n_components = int(key.strip(prefix))
                self._cache[key] = func(n_components=n_components, 
                                        filename=key+'.csv')
                return self._cache[key]
        
        for prefix, func in [('_dres_pc', self.get_dres),
                            ]:
            if prefix in key:
                dom, npc = key.split(prefix)
                n_components = int(npc)

                self._cache[key] = func(n_components=n_components, 
                                        col=self.short2long[dom],
                                        filename=key+'.csv')
                return self._cache[key]

        raise AttributeError(f'No attribute {key}')

    @staticmethod
    def pc_cols(n_components):
        return [f'PC {i+1}' for i in range(n_components)]
    
    @cached
    def cpp_pca_wide(self):
        return self.get_cppdf_wide(filename='cpp_pca_wide.csv')
    
    @cached
    def cpp_pca(self):
        return self.get_cppdf(filename='cpp_pca.csv')

    @cached
    def variance(self):
        return self.get_variance(filename='variance.csv')

    @cached
    def _individual_pc(self):
        return self.get_pca(filename=('traj_components.csv',
                                      'traj_variance.csv'))
    
    @cached
    def traj_components(self):
        return self._individual_pc[0]
    
    @cached
    def traj_variance(self):
        return self._individual_pc[1]

    def get_cppdf_wide(self, *args, **kwargs):
        names = ['Frame'] + self.pc_cols(self.n_components)
        fn = self.cpp_path('proj.dat')
        df = pd.read_csv(fn, sep='\s+', header=None, index_col=0,
                         skiprows=1, names=names)
        df['Time (ns)'] = np.tile(self.time_ns, self.n_traj)
        df['Simulation'] = self._data.sim_labels_rep
        df['Replicate'] = self._data.rep_labels_rep
        df['Force field'] = self._data.ff_labels_rep
        df['Half'] = [self.halves[int(x >= df['Time (ns)'].max()/2)]
                      for x in df['Time (ns)']]
        return df

    def get_cppdf(self, *args, **kwargs):
        tidy = pd.melt(self.cpp_pca_wide, var_name='PC', value_name='Value',
                       id_vars=['Time (ns)', 'Simulation', 'Force field', 
                                'Replicate', 'Half'])
        return tidy
    
    def get_variance(self, **kwargs):
        path = self.cpp_path('evalfrac.dat')
        cols = ['PC', 'Fraction', 'Cumulative', 'Eigenvalue']
        df = pd.read_csv(path, sep='\s+', header=None,
                         skiprows=1, names=cols)
        df['% variance'] = df['Fraction']*100
        df['Variance'] = np.cumsum(df.Fraction.values)
        return df

    def get_traj_pca(self, ff, start=None, stop=None, verbose=False, **kwargs):
        na3 = self.n_atoms*3
        n_reps = len(ff.reps)
        # last row is eigenvals
        pca_arr = np.zeros((n_reps, na3+1, na3))
        for i, u in enumerate(ff.reps):
            p = pca.PCA(u, select='all', align=True, verbose=verbose)
            p.run(start=start, stop=stop)
            # var: (n_components, )
            # p_comp: (n_atoms*3, n_components)
            pca_arr[i] = np.r_[p.p_components, [p.variance]]
        index = np.tile(list(np.arange(na3)) + ['Fraction'], n_reps)
        arr = np.concatenate(pca_arr)
        labels = np.repeat(ff.sim_labels, na3+1)
        cols = ['PC {}'.format(i+1) for i in range(na3)]
        df = pd.DataFrame(arr, index=index, columns=cols)
        df['Simulation'] = labels
        df['Force field'] = ff.name
        return df

    def get_pca(self, verbose=False, **kwargs):
        import gc
        half = int(np.ceil(self.n_frames/2))
        components = None
        variance = None
        for ff in self._data.forcefields:
            evects_ = []
            evals_ = []
            for start, stop, time in ((None, half, self._h1),
                                      (half, None, self._h2),
                                      (None, None, self._h0),):
                fn = f'{ff._ff}_{time}_pca.csv'
                df = self.get_traj_pca(ff, start=start, stop=stop,
                                             verbose=verbose,
                                             filename=fn)
                df['Range'] = time
                mask = df.index == 'Fraction'
                # evects_.append(df[~mask])
                # evals_.append(df[mask])
                print(df.values)

                if components is None:
                    components = df[~mask]
                else:
                    components = pd.concat([components, df[~mask]])
                
                if variance is None:
                    variance = df[mask]
                else:
                    variance = pd.concat([variance, df[mask]])

                del df

            # evects.append(pd.concat(evects_))
            # del evects_
            # evals.append(pd.concat(evals_))
            # del evals_
            # gc.collect()
            
        # components = pd.concat(evects)
        # del evects
        # variance = pd.concat(evals)
        # del evals

        return components, variance
    
    def get_rmsip(self, n_components=10, verbose=False, **kwargs):
        df = self.traj_components
        pc_cols = self.pc_cols(n_components)
        
        arr = np.zeros((self.n_traj, self.n_traj))
        for i, name in enumerate(self.sim_names):
            if verbose:
                print(f'Computing for {name}')
            sdf = df[df.Simulation==name]

            # (n_features, n_components) => (n_components, n_features)
            h0 = sdf[sdf.Range==self._h0][pc_cols].values.T
            h1 = sdf[sdf.Range==self._h1][pc_cols].values.T
            h2 = sdf[sdf.Range==self._h2][pc_cols].values.T
            arr[i, i] = pca.rmsip(h1, h2)

            for j, name in enumerate(self.sim_names[i+1:], i+1):
                j0 = df[(df.Simulation==name) & 
                        (df.Range==self._h0)][pc_cols].values.T
                arr[i, j] = arr[j, i] = pca.rmsip(h0, j0)
            
        df = pd.DataFrame(arr, index=self.sim_names, columns=self.sim_names)
        return df

    def get_rmsip_norm(self, n_components=10, **kwargs):
        df = getattr(self, f'rmsip_pc{n_components}')
        arr = df.values
        diag = np.diag_indices_from(arr)
        norm_factor = np.outer(arr[diag], arr[diag]) ** 0.5
        normalised = arr/norm_factor
        normalised[diag] = 0
        df = pd.DataFrame(arr, index=self.sim_names, columns=self.sim_names)
        return df

    def get_dres(self, col='Simulation', n_components=10, n_resample=1000,
                 verbose=False, **kwargs):
        df = self.cpp_pca_wide
        pc_cols = self.pc_cols(n_components)
        labels = df[col].unique()

        # get gaussian kde
        kdes = []
        resamples = []
        pdfs = []
        logpdfs = []
        div = np.zeros((len(labels), len(labels)))

        for j, name in enumerate(labels):
            if verbose:
                print('  Computing for ', name)
            df_ = df[df[col]==name][pc_cols]
            # transpose to (n_components, n_samples)
            dist_a = scipy.stats.gaussian_kde(df_.values.T)
            kdes.append(dist_a)
            res_a = dist_a.resample(size=n_resample)
            resamples.append(res_a)

            p_pdf = dist_a.evaluate(res_a)
            pdfs.append(p_pdf)
            ln_p = np.log(p_pdf).mean()
            logpdfs.append(ln_p)

            for i in range(j):
                dist_b = kdes[i]
                res_b = resamples[i]
                q_pdf = pdfs[i]
                ln_q = logpdfs[i]

                ln_pq_exp_p = np.log(0.5*(p_pdf+dist_b.evaluate(res_a))).mean()
                ln_pq_exp_q = np.log(0.5*(q_pdf+dist_a.evaluate(res_b))).mean()

                js = 0.5*(ln_p - ln_pq_exp_p + ln_q - ln_pq_exp_q)

                div[i, j] = div[j, i] = js
            
        df = pd.DataFrame(div, index=labels, columns=labels)
        return df

    def plot_points(self, filename='points.tif', x=1, y=2,
                    figsize=(10, 10), hue='Simulation', legend=False,
                    alpha=1, linewidth=0, s=4, edgecolor='none',
                    **kwargs):
        fig = plt.figure(figsize=figsize)
        if hue == 'Simulation':
            palette = self._data.color_labels
        else:
            palette = self._data.colors

        xlabel = f'PC {x}'
        ylabel = f'PC {y}'
        
        ax = sns.scatterplot(x=xlabel, y=ylabel, data=self.cpp_pca_wide, hue=hue,
                             palette=palette, ci=None, legend=legend, alpha=alpha,
                             linewidth=linewidth, s=s, edgecolor=edgecolor,
                             **kwargs)
        varx = self.variance[self.variance.PC==x]['% variance'].values[0]
        vary = self.variance[self.variance.PC==y]['% variance'].values[0]
        ax.set_xlabel(f'{xlabel} ({varx:.2f}% variance)')
        ax.set_ylabel(f'{ylabel} ({vary:.2f}% variance)')

        filename = f'pc{x}_pc{y}_{filename}'
        self.save_fig(filename)
        return ax
        
    
    def plot_dist(self, filename='dist.tif', ymax=0.03,
                  adjust_left=0.2, legend_loc='upper left',
                  bbox_to_anchor=(1.3, 1), n_yticks=5,
                  shade=False, hue='Force field',
                  n_components=10, col_wrap=5, row=None,
                  height=3, aspect=1.3, sharex=False, sharey=False,
                  margin_titles=False, legend=True):
        
        def plot_kde(data=None, color=None, label=None):
            sns.kdeplot(data['Value'], color=color, label=label, shade=shade)
        
        pcs = self.pc_cols(n_components)
        df = self.cpp_pca[self.cpp_pca.PC.isin(pcs)]
        if hue == 'Force field':
            palette = self.ff_colors
            hue_order = self.ff_names
            filename = 'ff_' + filename
        else:
            palette = self.sim_colors
            hue_order = self.sim_names
            filename = 'sim_' + filename

        g = sns.FacetGrid(df, col='PC', row=row, hue=hue,
                          col_wrap=col_wrap, palette=palette,
                          hue_order=hue_order, 
                          height=height, aspect=aspect,
                          sharex=sharex, sharey=sharey,
                          margin_titles=margin_titles)
        g.map_dataframe(plot_kde)
        g = self.reset_margin_titles(g)

        if legend:
            plt.subplots_adjust(left=adjust_left)
            plt.legend(bbox_to_anchor=bbox_to_anchor, loc=legend_loc,
                       borderaxespad=0.)
        self.save_fig(filename, g=g)
        return g
    
    def plot_violins(self, filename='viol.tif', x='Force field', 
                     n_components=10, col_wrap=5, row=None, orient='h',
                     height=4, aspect=1.3, sharex=False, sharey=True,
                     margin_titles=False, inner=None, width=1, scale_hue=True,
                     xticks=True, split=True, split_palette=True,
                     hue='Half'):
                
        pcs = self.pc_cols(n_components)
        df = self.cpp_pca[self.cpp_pca.PC.isin(pcs)]

        if row:
            col_wrap=None

        if hue == 'Half':
            hue_order = self.halves
        else:
            hue_order = np.arange(1, 5)

        g = sns.FacetGrid(df, col='PC', row=row,
                          col_wrap=col_wrap,
                          height=height, aspect=aspect,
                          sharex=sharex, sharey=sharey,
                          margin_titles=margin_titles)
        xticks = self.ff_names if xticks else []
        self.violinplot(g=g, x=x, y='Value', hue=hue, split=split,
                        hue_order=hue_order,
                        split_palette=split_palette, palette=self.ff_colors, 
                        order=self.ff_names, xticks=xticks)
        g = self.reset_margin_titles(g)
        self.save_fig(filename, g=g)
        return g

    
    def plot_scree(self, filename='scree.tif', y='Variance',
                   n_components=10, ymax=1):
        ax = sns.lineplot(data=self.variance[:n_components],
                          x='PC', y=y)
        plt.ylim(bottom=0, top=ymax)
        filename = y.lower() + '_' + filename
        self.save_fig(filename)
        return ax

    def plot_dres_rmsip(self, n_components=10, figsize=(24, 18),
                        dres_cmap='gist_heat', rmsip_cmap='bone_r',
                        linecolor='white', linewidths=0.1,
                        filename='dres_rmsip.tif', 
                        round_vmax=True, vmin=None,
                        **kwargs):
        
        dres = getattr(self, f'sim_dres_pc{n_components}').values
        rmsip = getattr(self, f'rmsip_pc{n_components}').values
        norm_rmsip = getattr(self, f'rmsip_norm_pc{n_components}').values

        arr = np.zeros((self.n_traj+2, self.n_traj))
        iu = np.triu_indices(self.n_traj, k=1)
        il = np.tril_indices(self.n_traj, k=-1)

        arr[iu] = dres[iu]
        arr[il] = rmsip[il]
        arr[-1] = np.diagonal(norm_rmsip)

        fig, ax = plt.subplots(figsize=figsize)
        labels = self.sim_names + ['', 'Self']

        # plot dres
        mask = np.ones_like(arr, dtype=bool)
        mask[iu] = False
        mask[np.diag_indices(mask.shape[1])] = False
        variance = self.variance.Variance[n_components] * 100
        dres_label = 'Jensen-Shannon divergence'

        vmax = 0.7 if round_vmax else np.log(2)
        if vmin is None:
            vmin = 0
        elif vmin == True:
            vmin = max(np.round(dres[iu].min(), decimals=1)-0.1, 0)

        sns.heatmap(arr, ax=ax, vmin=vmin, vmax=vmax,mask=mask,
                    square=True, cmap=dres_cmap, cbar=True, 
                    linecolor=linecolor, linewidths=linewidths,
                    cbar_kws={'label': dres_label})
        
        # plot rmsip
        mask = ~mask
        mask[-2] = True
        rmsip_label = 'RMSIP'

        sns.heatmap(arr, ax=ax, cmap=rmsip_cmap, vmin=0, vmax=1,
                    square=True, cbar=True, mask=mask,
                    linecolor=linecolor, linewidths=linewidths, 
                    xticklabels=self.sim_names, yticklabels=labels,
                    cbar_kws={'label': rmsip_label}) 

        title = f'{n_components:d} components ({variance:.0f}% variance)'
        plt.title(title)
        plt.tight_layout()

        path = filename[:-4] + f'_pc{n_components}' + filename[-4:]
        self.save_fig(path)
        return fig

    def plot_variance(self, filename='individual_variance.tif', 
                      n_components=10, plot_line=0.85,
                      height=1, aspect=1.3, sharex=True, sharey='row',
                      margin_titles=True, ylim=(0, 1),
                      alpha=0.4, **kwargs):

        def lineplot(data=None, color=None, label=None, **kwargs):
            var = np.cumsum(data.Fraction)
            sns.lineplot(data.PC, var, color=color, label=label)

        pc_cols = self.pc_cols(n_components)
        others = ['Force field', 'Range', 'Simulation']
        df = pd.melt(self.traj_variance, id_vars=others,
                     value_vars=pc_cols, var_name='PC',
                     value_name='Fraction')
        
        g = sns.FacetGrid(df, row='Range', col='Force field',
                          hue='Simulation', hue_order=self.sim_names,
                          row_order=[self.h0, self.h1, self.h2],
                          sharey=sharey, sharex=sharex, height=height,
                          aspect=aspect, margin_titles=margin_titles,
                          ylim=ylim)
        
        if plot_line:
            g.map(plt.axhline, y=plot_line, color='k', alpha=alpha)
        
        g.map(lineplot)
        g.reset_margin_titles(g)
        plt.tight_layout()
        self.save_fig(filename, g=g)
        return g
    
    def plot_dres(self, n_components=10, col='Simulation',
                  square=True, fmt='.2g', cmap='gist_heat',
                  linewidths=0.1, linecolor='k',
                  cbar_label='Jensen-Shannon divergence',
                  filename='dres.tif', figscale=3,
                  round_vmax=True, vmin=None, spines=True,
                  **kwargs):
        _col = self.long2short[col]
        dres = getattr(self, f'{_col}_dres_pc{n_components}')

        fig = plt.figure()
        size = fig.get_size_inches()
        fig = plt.figure(figsize=(figscale*size[0], figscale*size[1]))

        vmax = 0.7 if round_vmax else np.log(2)
        if vmin is None:
            vmin = 0
        elif vmin == True:
            iu = np.triu_indices(len(dres), k=1)
            vmin = max(np.round(dres.values[iu].min(), decimals=1)-0.1, 0)

        ax = sns.heatmap(dres, vmin=vmin, vmax=vmax, square=square,
                         linewidths=linewidths, linecolor=linecolor,
                         cmap=cmap, fmt=fmt, cbar_kws={'label': cbar_label})
        
        ax.spines['top'].set_visible(spines)
        ax.spines['right'].set_visible(spines)
        ax.spines['left'].set_visible(spines)
        ax.spines['bottom'].set_visible(spines)

        variance = self.variance.Variance[n_components] * 100
        title = f'{n_components:d} components ({variance:.0f}% variance)'
        plt.title(title)

        path = filename[:-4] + f'_pc{n_components}' + filename[-4:]
        self.save_fig(path)
        return ax
        






        

    

