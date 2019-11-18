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
    rwind = h5f['run_config']['rwind'][()]
    return rwind + (t - 1.0)




def configure_axes(ax1, ax2, ax3, filename, focus=0.0):
    h5f = h5py.File(filename, 'r')
    lc = light_pulse_position(h5f)

    for ax in [ax1, ax2, ax3]:
        ax.axvline(lc, c='k', ls='--', lw=0.5)

        if focus:
            ax.set_xlim(lc / (1 + focus), lc * (1 + 0.01))
            ax.set_yscale('log')
        else:
            # ax.set_xlim(1.0, 1000.0)
            ax.set_xscale('log')
            ax.set_yscale('log')

    # ax1.set_ylim(1e-4, 1e2)
    # ax2.set_ylim(1e-8, 1e3)
    ax3.set_ylim(1e-2, 2e2)
    ax1.set_ylabel(r'$\rho$')
    ax2.set_ylabel(r'$p$')
    ax3.set_ylabel(r'$\Gamma \beta$')
    ax3.set_xlabel(r'Radius $(r / r_{\rm inner})$')




def plot_variables(ax1, ax2, ax3, filename, edges=False):
    h5f = h5py.File(filename, 'r')
    vertices   = h5f['vertices']
    density    = h5f['primitive'][...][:,0]
    gamma_beta = h5f['primitive'][...][:,1]
    pressure   = h5f['primitive'][...][:,4]

    ax1.step(vertices[:-1], density)
    ax2.step(vertices[:-1], pressure)
    ax3.step(vertices[:-1], gamma_beta)

    if edges:
        for v in vertices[:-1]:
            ax1 .axvline(v, ls='-', c='k', lw=0.5)




def plot_radial(fig, filenames, focus=0.0, edges=False):
    ax1, ax2, ax3 = fig.subplots(nrows=3, ncols=1)

    for filename in filenames:
        plot_variables(ax1, ax2, ax3, filename, edges=edges)

        if filename is filenames[0]:
            configure_axes(ax1, ax2, ax3, filename, focus=focus)




def make_movie(fig, directories, output, **kwargs):
    filename_lists = [directory_contents(d) for d in directories]
    writer = FFMpegWriter(fps=10)

    with writer.saving(fig, output, dpi=200):
        for filename_list in zip(*filename_lists):
            print(list(filename_list))
            plot_radial(fig, filename_list, **kwargs)
            writer.grab_frame()
            fig.clf()




if __name__ == "__main__":
    """
    1. Program generates a figure with each filename argument plotted as a
       single curve, all on the same figure

    2. Program generates a movie by zipping the contents of each of the
       directory arguments, plotting corresponding files on the same frame
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    parser.add_argument("--output", default='output.mp4')
    parser.add_argument("--edges", action='store_true')
    parser.add_argument("--focus", type=float, default=0.0)
    args = parser.parse_args()

    fig = plt.figure(figsize=[12, 9])

    if all([os.path.isdir(p) for p in args.filenames]):
        make_movie(fig, args.filenames, args.output, focus=args.focus)

    else:
        plot_radial(fig, args.filenames, edges=args.edges, focus=args.focus)
        plt.show()
