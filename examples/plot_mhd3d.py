#!/usr/bin/env python3




import os
import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter




hydro_variables = ['density', 'vx', 'vy', 'vz', 'pressure', 'bx', 'by', 'bz']
field_variables = dict(Bx='magnetic_flux_1', By='magnetic_flux_2', Bz='magnetic_flux_3')




def load_slice(h5f, axis, variable):
    if variable in hydro_variables:
        ivar = hydro_variables.index(variable)
        if axis == 'x':
            n = h5f['primitive'].shape[0]
            return h5f['primitive'][n // 2, :, :][..., ivar]
        if axis == 'y':
            n = h5f['primitive'].shape[1]
            return h5f['primitive'][:, n // 2, :][..., ivar]
        if axis == 'z':
            n = h5f['primitive'].shape[2]
            return h5f['primitive'][:, :, n // 2][..., ivar]

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




def plot_slice(fig, filename, field=None, **kwargs):
    (ax1, ax2, ax3), (cb1, cb2, cb3) = fig.subplots(nrows=2, ncols=3, gridspec_kw={'height_ratios': [19, 1]})

    h5f = h5py.File(filename, 'r')

    print("total mass ..... ", h5f['conserved'][...][..., 0].mean())
    print("total px ....... ", h5f['conserved'][...][..., 1].mean())
    print("total py ....... ", h5f['conserved'][...][..., 2].mean())
    print("total pz ....... ", h5f['conserved'][...][..., 3].mean())
    print("total energy ... ", h5f['conserved'][...][..., 4].mean())

    dx = load_slice(h5f, 'x', field).T
    dy = load_slice(h5f, 'y', field).T
    dz = load_slice(h5f, 'z', field).T
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




def make_figure():
    fig = plt.figure(figsize=[16, 8])
    fig.subplots_adjust(left=0.02, right=0.98)
    return fig




if __name__ == "__main__":
    plots = dict(profile=plot_slice)

    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    parser.add_argument("--output", default='output.mp4')
    parser.add_argument("--field", choices=hydro_variables + list(field_variables.keys()), default='density')
    parser.add_argument("--plot", choices=plots.keys(), default='profile')
    args = parser.parse_args()


    for filename in args.filenames:
        plots[args.plot](make_figure(), filename, **vars(args))

    plt.show()
