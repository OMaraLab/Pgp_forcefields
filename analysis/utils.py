import os
import functools

import numpy as np
import pandas as pd
from MDAnalysis.lib.util import asiterable, iterable


def get_path(path):
    root = os.path.abspath(os.getcwd()).split('pgp')[0] + 'pgp/'
    return os.path.join(root, path)

def analysis_path(path):
    path = os.path.join('analysis', path)
    return get_path(path)

def data_path(path):
    path = os.path.join('data', path)
    return get_path(path)

def cached(func):
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        key = func.__name__
        try:
            return self._cache[key]
        except KeyError:
            self._cache[key] = ret = func(self, *args, **kwargs)
            return ret

    return property(wrapper)

def try_mkdir(cwd, dirname=None):
    if dirname:
        path = os.path.join(cwd, dirname)
        try:
            os.mkdir(path)
        except:
            pass
        return path
    return cwd

def read_csv(path):
    return pd.read_csv(path, index_col=0, header=0)

def png(func):
    def wrapper(self, filename, png=False, dpi=300, **kwargs):
        if png:
            filename = filename[:-3] + 'png'
            dpi = 120
        return wrapper(self, filename, dpi=dpi, **kwargs)
    return wrapper

def force_verbose(func):
    def wrapper(self, *args, verbose=None, force=None, **kwargs):
        verbose = verbose if verbose is not None else self._verbose
        force = force if force is not None else self._force
        return func(self, *args, verbose=verbose, force=force, **kwargs)
    return wrapper

def datafile(func):
    @functools.wraps(func)
    def wrapper(self, *args, verbose=False, force=False, **kwargs):
        filename = kwargs.get('filename')
        if filename is None:
            return func(self, *args, verbose=verbose, force=force, **kwargs)
        fns = asiterable(filename)
        paths = []
        data = []
        for fn in fns:
            datum, path = self.try_get_data(fn, force=force, verbose=verbose)
            paths.append(path)
            data.append(datum)

        if all(d is not None for d in data):
            if len(data) == 1:
                data = data[0]
            return data

        data = func(self, *args, verbose=verbose, force=force, **kwargs)
        if (not iterable(data) or 
            isinstance(data, (np.ndarray, pd.DataFrame))):
            data = [data]
        
        for datum, path in zip(data, paths):
            comments = None
            if path.endswith('npy'):
                try:
                    datum, comments = datum
                except ValueError:
                    pass
            self.save_data(datum, path, comments=comments)
        
        if len(data) == 1:
            data = data[0]
        return data
    return wrapper