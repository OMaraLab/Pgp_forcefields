from typing import Callable, Union

import seaborn as sns
import numpy as np
import pandas as pd
import colorsys
import matplotlib as mpl
import matplotlib.pyplot as plt
from seaborn.categorical import _ViolinPlotter
from seaborn.utils import remove_na

class SplitViolinPlotter(_ViolinPlotter):
    
    def __init__(self, *args, hatch=None, split_palette=False,
                 **kwargs):
        self.hatch = hatch
        self.split_palette = split_palette
        super().__init__(*args, **kwargs)
        
    def establish_colors(self, color, palette, saturation):
        if not self.split_palette:
            super().establish_colors(color, palette, saturation)
            # self.colors = [self.colors]* len(self.group_names)
            return
        n_colors = len(self.hue_names)
        if self.hatch:
            rgb_colors = [[x]*n_colors for x in sns.color_palette(palette)]
        else:
            rgb_colors = []
            for h in palette:
                pal = sns.blend_palette(['#ffffff', h, '#000000'], n_colors=n_colors*2+1)
                rgb_colors.append(pal[1::2])
        
        light_vals = [colorsys.rgb_to_hls(*c[-1][:3])[1] for c in rgb_colors]
        lum = min(light_vals) * .6
        gray = mpl.colors.rgb2hex((lum, lum, lum))

        # Assign object attributes
        self.colors = rgb_colors
        self.gray = gray
        
    def add_legend_data(self, ax, color, label, hatch=None):
        """Add a dummy patch object so we can get legend data."""
        rect = plt.Rectangle([0, 0], 0, 0,
                             linewidth=self.linewidth / 2,
                             edgecolor=self.gray,
                             facecolor=color,
                             label=label,
                             hatch=hatch)
        ax.add_patch(rect)
        
    def draw_violins(self, ax):
        fill_func = ax.fill_betweenx if self.orient == "v" else ax.fill_between
        for i, group_data in enumerate(self.plot_data):

            kws = dict(edgecolor=self.gray, linewidth=self.linewidth)

            # Option 1: we have a single level of grouping
            # --------------------------------------------

            if self.plot_hues is None:

                support, density = self.support[i], self.density[i]

                # Handle special case of no observations in this bin
                if support.size == 0:
                    continue

                # Handle special case of a single observation
                elif support.size == 1:
                    val = support.item()
                    d = density.item()
                    self.draw_single_observation(ax, i, val, d)
                    continue

                # Draw the violin for this group
                grid = np.ones(self.gridsize) * i
                fill_func(support,
                          grid - density * self.dwidth,
                          grid + density * self.dwidth,
                          facecolor=self.colors[i],
                          **kws)

                # Draw the interior representation of the data
                if self.inner is None:
                    continue

                # Get a nan-free vector of datapoints
                violin_data = remove_na(group_data)

                # Draw box and whisker information
                if self.inner.startswith("box"):
                    self.draw_box_lines(ax, violin_data, support, density, i)

                # Draw quartile lines
                elif self.inner.startswith("quart"):
                    self.draw_quartiles(ax, violin_data, support, density, i)

                # Draw stick observations
                elif self.inner.startswith("stick"):
                    self.draw_stick_lines(ax, violin_data, support, density, i)

                # Draw point observations
                elif self.inner.startswith("point"):
                    self.draw_points(ax, violin_data, i)

            # Option 2: we have nested grouping by a hue variable
            # ---------------------------------------------------

            else:
                offsets = self.hue_offsets
                for j, hue_level in enumerate(self.hue_names):

                    support, density = self.support[i][j], self.density[i][j]
                    kws["facecolor"] = self.colors[i][j]
                    
                    if self.hatch:
                        hatch = self.hatch[j]
                    else:
                        hatch = None

                    # Add legend data, but just for one set of violins
                    if not i and (self.hatch or not self.split_palette):
                        if self.hatch:
                            lhatch = hatch
                        else:
                            lhatch = None
                        
                        if self.split_palette:
                            hcol = '#ffffff'
                        else:
                            hcol = self.colors[i][j]
                    
                        self.add_legend_data(ax, hcol, hue_level, hatch=lhatch)

                    # Handle the special case where we have no observations
                    if support.size == 0:
                        continue

                    # Handle the special case where we have one observation
                    elif support.size == 1:
                        val = support.item()
                        d = density.item()
                        if self.split:
                            d = d / 2
                        at_group = i + offsets[j]
                        self.draw_single_observation(ax, at_group, val, d)
                        continue

                    # Option 2a: we are drawing a single split violin
                    # -----------------------------------------------
                    
                    grid = np.ones(self.gridsize) * i
                    if j:
                        fill_func(support,
                                grid,
                                grid + density * self.dwidth,
                                hatch=hatch,
                                **kws)
                    else:
                        fill_func(support,
                                grid - density * self.dwidth,
                                grid,
                                hatch=hatch,
                                **kws)

                    # Draw the interior representation of the data
                    if self.inner is None:
                        continue

                    # Get a nan-free vector of datapoints
                    hue_mask = self.plot_hues[i] == hue_level
                    violin_data = remove_na(group_data[hue_mask])

                    # Draw quartile lines
                    if self.inner.startswith("quart"):
                        self.draw_quartiles(ax, violin_data,
                                            support, density, i,
                                            ["left", "right"][j])

                    # Draw stick observations
                    elif self.inner.startswith("stick"):
                        self.draw_stick_lines(ax, violin_data,
                                            support, density, i,
                                            ["left", "right"][j])

                    # The box and point interior plots are drawn for
                    # all data at the group level, so we just do that once
                    if not j:
                        continue

                    # Get the whole vector for this group level
                    violin_data = remove_na(group_data)

                    # Draw box and whisker information
                    if self.inner.startswith("box"):
                        self.draw_box_lines(ax, violin_data,
                                            support, density, i)

                    # Draw point observations
                    elif self.inner.startswith("point"):
                        self.draw_points(ax, violin_data, i)

                    # Option 2b: we are drawing full nested violins
                    # -----------------------------------------------

                    else:
                        grid = np.ones(self.gridsize) * (i + offsets[j])
                        fill_func(support,
                                  grid - density * self.dwidth,
                                  grid + density * self.dwidth,
                                  **kws)

                        # Draw the interior representation
                        if self.inner is None:
                            continue

                        # Get a nan-free vector of datapoints
                        hue_mask = self.plot_hues[i] == hue_level
                        violin_data = remove_na(group_data[hue_mask])

                        # Draw box and whisker information
                        if self.inner.startswith("box"):
                            self.draw_box_lines(ax, violin_data,
                                                support, density,
                                                i + offsets[j])

                        # Draw quartile lines
                        elif self.inner.startswith("quart"):
                            self.draw_quartiles(ax, violin_data,
                                                support, density,
                                                i + offsets[j])

                        # Draw stick observations
                        elif self.inner.startswith("stick"):
                            self.draw_stick_lines(ax, violin_data,
                                                  support, density,
                                                  i + offsets[j])

                        # Draw point observations
                        elif self.inner.startswith("point"):
                            self.draw_points(ax, violin_data, i + offsets[j])



def violinplot(x=None, y=None, hue=None, data=None, order=None, hue_order=None,
               bw="scott", cut=2, scale="area", scale_hue=True, gridsize=100,
               width=.8, inner="box", split=False, dodge=True, orient=None,
               linewidth=None, color=None, palette=None, saturation=.75,
               ax=None, hatch=None, split_palette=False, **kwargs):

    plotter = SplitViolinPlotter(x, y, hue, data, order, hue_order,
                             bw, cut, scale, scale_hue, gridsize,
                             width, inner, split, dodge, orient, linewidth,
                             color, palette, saturation, hatch=hatch,
                             split_palette=split_palette)

    if ax is None:
        ax = plt.gca()

    plotter.plot(ax)
    return ax

def kdeplot(data: pd.DataFrame=None, x: str=None, color: str=None, 
            label: str=None, **kwargs):
    sns.kdeplot(data[x], color=color, label=label, **kwargs)

def axspan(data: pd.DataFrame=None, x: str=None, y: str=None, orient: str='h',
           **kwargs):
    func = plt.axvspan if orient=='h' else plt.axhspan
    low = data[x].min()
    if y is None:
        y = x
    high = data[y].max()
    func(low, high, **kwargs)

def axline(data: pd.DataFrame=None, x: Union[Callable, str]=None,
           color: str=None, label: str=None, orient: str='h', **kwargs):
    if isinstance(x, str):
        values = data[x].unique()
    else:
        values = x(data)

    func = plt.axvline if orient=='h' else plt.axhline

    for v in values:
        func(v, color=color, **kwargs)
    

    
