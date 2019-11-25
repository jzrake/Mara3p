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




def load_slice(h5f, axis, variable):

    variables = ['density', 'px', 'py', 'pz', 'energy']

    if variable in variables:
        ivar = variables.index(variable)
        if axis == 'x':
            n = h5f['conserved'].shape[0]
            return h5f['conserved'][n // 2, :, :][..., ivar]
        if axis == 'y':
            n = h5f['conserved'].shape[1]
            return h5f['conserved'][:, n // 2, :][..., ivar]
        if axis == 'z':
            n = h5f['conserved'].shape[2]
            return h5f['conserved'][:, :, n // 2][..., ivar]

    if variable == 'flux':
        if axis == 'x':
            n = h5f['magnetic_flux_1'].shape[0]
            return h5f['magnetic_flux_1'][n // 2, :, :][...]
        if axis == 'y':
            n = h5f['magnetic_flux_2'].shape[1]
            return h5f['magnetic_flux_2'][:, n // 2, :][...]
        if axis == 'z':
            n = h5f['magnetic_flux_3'].shape[2]
            return h5f['magnetic_flux_3'][:, :, n // 2][...]




def plot_slice(fig, filenames, focus=0.0, edges=False):
    ax1, ax2, ax3 = fig.subplots(nrows=1, ncols=3)

    for filename in filenames:
        h5f = h5py.File(filename, 'r')

        dx = load_slice(h5f, 'x', 'energy')
        dy = load_slice(h5f, 'y', 'energy')
        dz = load_slice(h5f, 'z', 'energy')
        ax1.imshow(dx)
        ax2.imshow(dy)
        ax3.imshow(dz)

    ax1.set_title('X slice')
    ax2.set_title('Y slice')
    ax3.set_title('Z slice')

    for ax in [ax1, ax2, ax3]:
        ax.set_xticks([])
        ax.set_yticks([])




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

    fig = plt.figure(figsize=[15, 5])
    fig.subplots_adjust(top=1.0, bottom=0.0, left=0.01, right=0.99, hspace=0.0)

    if all([os.path.isdir(p) for p in args.filenames]):
        make_movie(fig, args.filenames, args.output, plots[args.plot])

    else:
        plots[args.plot](fig, args.filenames)
        plt.show()
