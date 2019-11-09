#!/usr/bin/env python3




import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt




def plot_radial(filename):
    h5f = h5py.File(filename, 'r')
    vertices = h5f['vertices']
    density = h5f['primitive'][...][:,0]

    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=[12, 10])
    ax1.step(vertices, density)
    ax1.set_xscale('log')
    ax1.set_xlabel(r'Radius $(r / r_{\rm inner})$')




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')

    args = parser.parse_args()

    for filename in args.filenames:
        plot_radial(filename)

    plt.show()
