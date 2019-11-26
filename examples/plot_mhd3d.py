#!/usr/bin/env python3




import os
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter




hydro_variables = ['density', 'px', 'py', 'pz', 'energy']
field_variables = dict(bx='magnetic_flux_1', by='magnetic_flux_2', bz='magnetic_flux_3')




def directory_contents(directory):
    """
    Return the contents of the given directory, sorted alphabetically, and with
    the directory name appended to each.
    """
    return filter(lambda f: f.endswith('.h5'), [os.path.join(directory, filename) for filename in sorted(os.listdir(directory))])




def load_slice(h5f, axis, variable):

    if variable in hydro_variables:
        ivar = hydro_variables.index(variable)
        if axis == 'x':
            n = h5f['conserved'].shape[0]
            return h5f['conserved'][n // 2, :, :][..., ivar]
        if axis == 'y':
            n = h5f['conserved'].shape[1]
            return h5f['conserved'][:, n // 2, :][..., ivar]
        if axis == 'z':
            n = h5f['conserved'].shape[2]
            return h5f['conserved'][:, :, n // 2][..., ivar]

    if variable in field_variables:
        kvar = field_variables[variable]
        if axis == 'x':
            n = h5f[kvar].shape[0]
            return h5f[kvar][n // 2, :, :][...]
        if axis == 'y':
            n = h5f[kvar].shape[1]
            return h5f[kvar][:, n // 2, :][...]
        if axis == 'z':
            n = h5f[kvar].shape[2]
            return h5f[kvar][:, :, n // 2][...]




def plot_slice(fig, filenames=[], field=None, **kwargs):
    (ax1, ax2, ax3), (cb1, cb2, cb3) = fig.subplots(nrows=2, ncols=3, gridspec_kw={'height_ratios': [19, 1]})

    for filename in filenames:
        h5f = h5py.File(filename, 'r')

        dx = load_slice(h5f, 'x', field)
        dy = load_slice(h5f, 'y', field)
        dz = load_slice(h5f, 'z', field)
        m1 = ax1.imshow(dx)
        m2 = ax2.imshow(dy)
        m3 = ax3.imshow(dz)

    for m, c in zip([m1, m2, m3], [cb1, cb2, cb3]):
        fig.colorbar(m, cax=c, orientation='horizontal')

    ax1.set_title('X slice')
    ax2.set_title('Y slice')
    ax3.set_title('Z slice')

    for ax in [ax1, ax2, ax3]:
        ax.set_xticks([])
        ax.set_yticks([])




def make_movie(fig, plot_function, output=None, directories=[], **kwargs):
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
    parser.add_argument("--field", choices=hydro_variables + list(field_variables.keys()), default='density')
    parser.add_argument("--plot", choices=plots.keys(), default='profile')
    args = parser.parse_args()

    fig = plt.figure(figsize=[16, 8])
    fig.subplots_adjust(left=0.02, right=0.98)

    if all([os.path.isdir(p) for p in args.filenames]):
        make_movie(fig, plots[args.plot], **vars(args))

    else:
        plots[args.plot](fig, **vars(args))
        plt.show()
