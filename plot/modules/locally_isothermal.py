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
from ..mesh import *




def print_quantiles(data):
    print()
    print('Quantiles:')
    for q in [0.0, 0.1, 1.0, 10.0, 90.0, 99.0, 99.9, 100.0]:
        print(r'{:.01f}%'.format(q).rjust(10), '.......', '{:.02e}'.format(np.percentile(data, q)))
    print()




def print_totals(data):
    print('total = {:.8f}'.format(data.mean()))




def field_names():
    return ['sigma', 'vx', 'vy']




def conserved(filename, field):
    if type(filename) is list:
        return np.array([conserved(f, field) for f in filename]).mean(axis=0)

    print('loading {}'.format(filename))

    blocks = []
    h5f = h5py.File(filename, 'r')

    for block in h5f['solution']['conserved']:
        blocks.append((block_index_from_string(block), h5f['solution']['conserved'][block][()][:,:,field]))

    depth = max(i[0] for i, v in blocks)
    mx = 1 << depth
    my = 1 << depth
    nx = max(v.shape[0] for i, v in blocks)
    ny = max(v.shape[1] for i, v in blocks)
    result = np.zeros([nx * mx, ny * my])

    for (level, i, j), v in blocks:
        bx = nx * (1 << (depth - level))
        by = ny * (1 << (depth - level))
        uv = upsample(v, depth - level)
        result[i*bx:i*bx+bx, j*by:j*by+by] = uv[...]

    return result




def formatted_time(filename):
    if len(filename) == 1:
        h5f = h5py.File(filename[0], 'r')
        return r'$t = {:.01f}$'.format(h5f['solution']['time'][()])
    orbit0 = h5py.File(filename[ 0], 'r')['solution']['time'][()]
    orbit1 = h5py.File(filename[-1], 'r')['solution']['time'][()]
    return r'$t = {:.01f} - {:0.1f}$'.format(orbit0, orbit1)




def plot_streamlines(ax, filename, args):
    sigma = conserved(filename, 0)
    vx = conserved(filename, 1) / sigma
    vy = conserved(filename, 2) / sigma
    R = h5py.File(filename[0], 'r')['run_config']['domain_radius'][()]
    X = np.linspace(-R, R, sigma.shape[0])
    Y = np.linspace(-R, R, sigma.shape[1])
    U = vx.T
    V = vy.T
    ax.streamplot(X, Y, U, V, density=3.0, linewidth=0.3, color='grey', arrowsize=0.6)




def plot(fig, filename, args):
    value = conserved(filename, field_names().index(args.field))
    color = value**args.index
    vmin, vmax = [a**args.index if a else a for a in eval(args.range or '[None,None]')]

    R = h5py.File(filename[0], 'r')['run_config']['domain_radius'][()]
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.imshow(color.T, cmap=args.cmap, vmin=vmin, vmax=vmax, origin='bottom', extent=[-R, R, -R, R])
    ax1.text(0.02, 0.02, formatted_time(filename), transform=ax1.transAxes, color='white')

    for x0, y0, x1, y1 in read_block_edges_2d(filename[0], domain_size=2 * R):
        ax1.plot([x0, x1], [y0, y1], lw=0.5, c='k')

    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_xlim(-R, R)
    ax1.set_ylim(-R, R)

    if args.streamlines > 0.0:
        plot_streamlines(ax1, filename, args)

    if args.print_quantiles:
        print_quantiles(value)
        print_totals(value)




def make_movie(fig, plot_fn, args):
    from matplotlib.animation import FFMpegWriter

    writer = FFMpegWriter(fps=15)

    with writer.saving(fig, 'output.mp4', dpi=200):
        for filename in args.filenames:
            print(filename)
            plot_fn(fig, filename, args)
            writer.grab_frame()
            fig.clf()

    print('writing', 'output.mp4')




def save_frames(fig, plot_fn, args):

    if args.archive and os.path.isfile('images.tar.gz'):
        raise IOError('images.tar.gz file already exists, please remove it')

    frames = []

    for filename in args.filenames:
        pngname = filename.replace('.h5', '.png')
        print('{} -> {}'.format(filename, pngname))
        frames.append(pngname)

        if os.path.isfile(pngname):
            continue

        plot_fn(fig, filename, args)
        fig.savefig(pngname)
        fig.clf()

    if args.archive:
        import tarfile

        with tarfile.open('images.tar.gz', "w:gz") as tar:
            for frame in frames:
                tar.add(frame, arcname=os.path.join('images', os.path.basename(frame)))

        for frame in frames:
            os.remove(frame)




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('--cmap',  metavar='inferno',       default='inferno',    help='Name of the matplotlib color map to use')
    parser.add_argument('--index', metavar=0.5, type=float, default=0.5,          help='Power-law index used in scaling the image')
    parser.add_argument('--range', metavar='[0.0,1.0]',     default=None,         help='Manual vmin/vmax range (use None for auto)')
    parser.add_argument('--field', choices=field_names(),   default='sigma',      help='Which hydrodynamic field to show')
    parser.add_argument('--print-quantiles', '-q', action='store_true',           help='Print a set of quantiles for the image data')
    parser.add_argument('--streamlines',           action='store_true',           help='Overlay a streamplot of the velocity field')
    parser.add_argument('--movie',                 action='store_true',           help='Use ffmpeg to make a movie from the provided files')
    parser.add_argument('--frame',                 action='store_true',           help='Write a PNG image for each file rather than raising a window')
    parser.add_argument('--no-window', '-nw',      action='store_true',           help='Do not plot anything')
    parser.add_argument('--archive',               action='store_true',           help='Compress images to images.tar.gz (implies --frame)')

    args = parser.parse_args()

    if args.frame or args.movie or args.no_window:
        import matplotlib
        matplotlib.use('agg')

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=[8, 8])
    fig.subplots_adjust(top=1, left=0, bottom=0, right=1)

    if args.movie:
        make_movie(fig, plot, args)
    elif args.frame or args.archive:
        save_frames(fig, plot, args)
    else:
        plot(fig, args.filenames, args)
        if not args.no_window:
            plt.show()
