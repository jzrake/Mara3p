#!/usr/bin/env python3




def conserved(h5f, field):
    blocks = []

    for block in h5f['solution']['conserved']:
        level, rest = block.split(':')
        index = [int(i) for i in rest.split('-')]
        blocks.append((index, h5f['solution']['conserved'][block]))

    mx = max(i[0] for i, v in blocks) + 1
    my = max(i[1] for i, v in blocks) + 1
    nx = max(v.shape[0] for i, v in blocks)
    ny = max(v.shape[1] for i, v in blocks)
    result = np.zeros([nx * mx, ny * my])

    for (i, j), v in blocks:
        result[i*nx:i*nx+nx, j*ny:j*ny+ny] = v[...][:,:,field]

    return result




def plot(fig, filename, args):
    h5f = h5py.File(filename, 'r')
    sigma = conserved(h5f, 0)
    vmin, vmax = [a**args.index if a else a for a in eval(args.range)]
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.imshow(sigma.T**args.index, cmap=args.cmap, vmin=vmin, vmax=vmax)
    ax1.text(0.02, 0.02, r'orbit {:.01f}'.format(h5f['solution']['time'][()] / 2 / np.pi), transform=ax1.transAxes, color='white')
    ax1.set_xticks([])
    ax1.set_yticks([])




def make_movie(fig, plot_fn, args):
    from matplotlib.animation import FFMpegWriter

    writer = FFMpegWriter(fps=15)

    with writer.saving(fig, args.output, dpi=200):
        for filename in args.filenames:
            print(filename)
            plot_fn(fig, filename, args)
            writer.grab_frame()
            fig.clf()

    print('writing', args.output)




def save_frames(fig, plot_fn, args):

    for filename in args.filenames:
        pngname = filename.replace('.h5', '.png')
        print('{} -> {}'.format(filename, pngname))
        plot_fn(fig, filename, args)
        fig.savefig(pngname)
        fig.clf()




if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-o', '--output', default='output.mp4')
    parser.add_argument('--movie', action='store_true')
    parser.add_argument('--frame', action='store_true')
    parser.add_argument('--cmap', default='inferno')
    parser.add_argument('--index', default=0.5, type=float)
    parser.add_argument('--range', default='[None,None]')
    args = parser.parse_args()

    if args.frame:
        import matplotlib
        matplotlib.use('agg')

    import h5py
    import numpy as np
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=[8, 8])
    fig.subplots_adjust(top=1, left=0, bottom=0, right=1)

    if args.movie:
        make_movie(fig, plot, args)
    elif args.frame:
        save_frames(fig, plot, args)
    else:
        plot(fig, args.filenames[0], args)
        plt.show()
