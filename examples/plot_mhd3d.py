#!/usr/bin/env python3




import os
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter




def directory_contents(directory):
    """
    Return the contents of the given directory, sorted alphabetically, and with
    the directory name appended to each.
    """
    return filter(lambda f: f.endswith('.h5'), [os.path.join(directory, filename) for filename in sorted(os.listdir(directory))])




def plot_variables(ax1, ax2, ax3, filename, edges=False, **kwargs):
    h5f = h5py.File(filename, 'r')
    density    = h5f['conserved'][...][0]

    ax1.step(vertices[:-1], density, **kwargs)
    ax2.step(vertices[:-1], pressure, **kwargs)
    ax3.step(vertices[:-1], gamma_beta, **kwargs)

    if edges:
        for v in vertices[:-1]:
            ax1 .axvline(v, ls='-', c='k', lw=0.5)




def plot_slice(fig, filenames, focus=0.0, edges=False):
    ax1 = fig.subplots(nrows=1, ncols=1)

    for filename in filenames:
        h5f = h5py.File(filename, 'r')
        nk = h5f['conserved'].shape[2]
        density = h5f['conserved'][:, :, nk // 2][:,:,0]
        ax1.imshow(density)




def make_movie(fig, directories, output, plot_function, **kwargs):
    filename_lists = [directory_contents(d) for d in directories]
    writer = FFMpegWriter(fps=10)

    with writer.saving(fig, output, dpi=200):
        for filename_list in zip(*filename_lists):
            print(list(filename_list))
            plot_function(fig, filename_list, **kwargs)
            writer.grab_frame()
            fig.clf()




if __name__ == "__main__":
    """
    1. Program generates a figure with each filename argument plotted as a
       single curve, all on the same figure

    2. Program generates a movie by zipping the contents of each of the
       directory arguments, plotting corresponding files on the same frame
    """
    plots = dict(profile=plot_slice)

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    parser.add_argument("--output", default='output.mp4')
    parser.add_argument("--plot", choices=plots.keys(), default='profile')
    args = parser.parse_args()

    fig = plt.figure(figsize=[12, 9])

    if all([os.path.isdir(p) for p in args.filenames]):
        make_movie(fig, args.filenames, args.output, plots[args.plot])

    else:
        plots[args.plot](fig, args.filenames)
        plt.show()
