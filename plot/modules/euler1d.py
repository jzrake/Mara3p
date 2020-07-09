"""
==============================================================================
Copyright 2020, Jonathan Zrake

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

==============================================================================
"""




import os
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from ..mesh import *




def plot(fig, args):
    h5f = h5py.File(args.filenames[0], 'r')
    ax1 = fig.add_subplot(1, 1, 1)

    x_array = []
    y_array = []

    for block, data in h5f['solution']['conserved'].items():
        level, index = block_index_from_string(block)
        x0, x1 = block_extent_1d(level, index)
        xv = np.linspace(x0, x1, data.shape[0] + 1)
        xc = 0.5 * (xv[1:] + xv[:-1])
        x_array += list(xc)
        y_array += list(data[...][:,0])

    ax1.plot(x_array, y_array, '-o', mfc='none')

    for edge in read_block_edges_1d(args.filenames[0]):
        ax1.axvline(edge, 0.0, 0.1, c='k', lw=0.5)




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()

    fig = plt.figure(figsize=[10, 8])
    plot(fig, args)
    plt.show()
