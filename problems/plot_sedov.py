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




def light_pulse_position(h5f):
    """
    Return the position of a light-pulse that left the wind-ambient medium
    interface at t=1 (the start of the simulation).
    """
    t = h5f['time'][()]
    engine_onset = h5f['run_config']['engine_onset'][()]
    return t - engine_onset + 1.0




def configure_axes(ax1, ax2, ax3, filename, focus=0.0):
    h5f = h5py.File(filename, 'r')
    lc = light_pulse_position(h5f)

    for ax in [ax1, ax2, ax3]:
        ax.axvline(lc, c='k', ls='--', lw=0.5)

        if focus:
            ax.set_xlim(lc / (1 + focus), lc * (1 + 0.01))
            ax.set_yscale('log')
        else:
            ax.set_xlim(1.0, 1000.0)
            ax.set_xscale('log')
            ax.set_yscale('log')

    # ax1.set_ylim(1e-4, 1e2)
    # ax2.set_ylim(1e-8, 1e3)
    ax3.set_ylim(1e-2, 2e2)
    ax1.set_ylabel(r'$\rho$')
    ax2.set_ylabel(r'$p$')
    ax3.set_ylabel(r'$\Gamma \beta$')
    ax3.set_xlabel(r'Radius $(r / r_{\rm inner})$')
    ax1.set_title(r'$t = {0:.1f}$'.format(h5f['time'][()]))




def plot_variables(ax1, ax2, ax3, filename, edges=False, **kwargs):
    h5f = h5py.File(filename, 'r')
    vertices   = h5f['vertices']
    density    = h5f['primitive'][...][:,0]
    gamma_beta = h5f['primitive'][...][:,1]
    pressure   = h5f['primitive'][...][:,4]

    ax1.step(vertices[:-1], density, **kwargs)
    ax2.step(vertices[:-1], pressure, **kwargs)
    ax3.step(vertices[:-1], gamma_beta, **kwargs)

    if edges:
        for v in vertices[:-1]:
            ax1 .axvline(v, ls='-', c='k', lw=0.5)




def plot_radial_profile(fig, filenames, focus=0.0, edges=False):
    ax1, ax2, ax3 = fig.subplots(nrows=3, ncols=1)

    for filename in filenames:
        plot_variables(ax1, ax2, ax3, filename, edges=edges, label=filename)

        if filename is filenames[0]:
            configure_axes(ax1, ax2, ax3, filename, focus=focus)
    ax1.legend()




def plot_cell_size_histogram(fig, filenames, **kwargs):
    for filename in filenames:

        h5f = h5py.File(filename, 'r')
        vertices = h5f['vertices']
        dr = np.diff(vertices)

        ax1 = fig.subplots(nrows=1, ncols=1)
        ax1.hist(np.log10(dr), bins=200, histtype='step')




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
    plots = dict(profile=plot_radial_profile, cells=plot_cell_size_histogram)

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    parser.add_argument("--output", default='output.mp4')
    parser.add_argument("--edges", action='store_true')
    parser.add_argument("--focus", type=float, default=0.0)
    parser.add_argument("--plot", choices=plots.keys(), default='profile')
    args = parser.parse_args()

    fig = plt.figure(figsize=[12, 9])

    if all([os.path.isdir(p) for p in args.filenames]):
        make_movie(fig, args.filenames, args.output, plots[args.plot], focus=args.focus)

    else:
        plots[args.plot](fig, args.filenames, edges=args.edges, focus=args.focus)
        plt.show()
