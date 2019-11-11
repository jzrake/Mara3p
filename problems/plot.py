#!/usr/bin/env python3




import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt




def plot_radial(ax1, ax2, filename):
    h5f = h5py.File(filename, 'r')
    vertices   = h5f['vertices']
    density    = h5f['primitive'][...][:,0]
    gamma_beta = h5f['primitive'][...][:,1]

    ax1.step(vertices[:-1], density)
    ax2.step(vertices[:-1], gamma_beta)

    # for x in vertices:
    #     ax2.axvline(x, c='k', lw=0.5)

    for ax in [ax1, ax2]:
        ax.set_xlim(1.0, 1000.0)

    ax1.set_ylim(1e-8, 1e2)
    ax2.set_ylim(1e-3, 2e2)
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.set_ylabel(r'$\rho$')
    ax2.set_ylabel(r'$\Gamma \beta$')
    ax2.set_xlabel(r'Radius $(r / r_{\rm inner})$')




def make_movie(fig, filenames, output):
    from matplotlib.animation import FFMpegWriter

    dpi = 200
    res = 768

    writer = FFMpegWriter(fps=10)

    with writer.saving(fig, output, dpi):
        for filename in filenames:
            print(filename)
            ax1, ax2 = fig.subplots(nrows=2, ncols=1)
            plot_radial(ax1, ax2, filename)
            writer.grab_frame()
            fig.clf()




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    parser.add_argument("--movie", action='store_true')
    parser.add_argument("--output", default='output.mp4')
    args = parser.parse_args()

    fig = plt.figure(figsize=[12, 9])

    if args.movie:
        make_movie(fig, args.filenames, args.output)
    else:
        ax1, ax2 = fig.subplots(nrows=2, ncols=1)
        for filename in args.filenames:
            plot_radial(ax1, ax2, filename)

        plt.show()
