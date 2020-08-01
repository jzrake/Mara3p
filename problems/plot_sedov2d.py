#!/usr/bin/env python3

import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon
from matplotlib import cm




def field_names():
    return ['rho', 'ur', 'ut', 'up', 'p']




def cell_triangulation(fname):
    h5f = h5py.File(fname, 'r')
    track_triangles_x = []
    track_triangles_y = []
    num_triangles = 0
    for track in h5f['tracks']:
        primitive  = h5f['tracks'][track]['primitive'][()]
        face_radii = h5f['tracks'][track]['face_radii'][()]
        theta0     = h5f['tracks'][track]['theta0'][()]
        theta1     = h5f['tracks'][track]['theta1'][()]
        x00 = np.sin(theta0) * face_radii[:-1]
        y00 = np.cos(theta0) * face_radii[:-1]
        x01 = np.sin(theta1) * face_radii[:-1]
        y01 = np.cos(theta1) * face_radii[:-1]
        x10 = np.sin(theta0) * face_radii[+1:]
        y10 = np.cos(theta0) * face_radii[+1:]
        x11 = np.sin(theta1) * face_radii[+1:]
        y11 = np.cos(theta1) * face_radii[+1:]
        xa = np.array([x00, x01, x11]).T
        xb = np.array([x11, x00, x10]).T
        ya = np.array([y00, y01, y11]).T
        yb = np.array([y11, y00, y10]).T
        track_triangles_x.append(np.vstack([xa, xb]))
        track_triangles_y.append(np.vstack([ya, yb]))
        num_triangles += (len(face_radii) - 1) * 2
    tx = np.vstack(track_triangles_x).flatten()
    ty = np.vstack(track_triangles_y).flatten()
    return tri.Triangulation(tx, ty, np.arange(num_triangles * 3).reshape(num_triangles, 3))




def get_cell_patches(fname):
    patches = []
    h5f = h5py.File(fname, 'r')
    for track in h5f['tracks']:
        face_radii = h5f['tracks'][track]['face_radii'][()]
        theta0     = h5f['tracks'][track]['theta0'][()]
        theta1     = h5f['tracks'][track]['theta1'][()]
        x00 = np.sin(theta0) * face_radii[:-1]
        y00 = np.cos(theta0) * face_radii[:-1]
        x01 = np.sin(theta1) * face_radii[:-1]
        y01 = np.cos(theta1) * face_radii[:-1]
        x10 = np.sin(theta0) * face_radii[+1:]
        y10 = np.cos(theta0) * face_radii[+1:]
        x11 = np.sin(theta1) * face_radii[+1:]
        y11 = np.cos(theta1) * face_radii[+1:]
        patches += [Polygon(coords.reshape(4, 2)) for coords in np.vstack([x00, y00, x01, y01, x11, y11, x10, y10]).T]
    return patches




def get_primitive_for_triangles(fname, field):
    h5f = h5py.File(fname, 'r')
    result = []
    num_triangles = 0
    i = field_names().index(field)
    for track in h5f['tracks']:
        primitive  = h5f['tracks'][track]['primitive'][()]
        s = primitive[:,i].tolist()
        result += s
        result += s
    return np.array(result)




def get_primitive_for_quads(fname, field):
    h5f = h5py.File(fname, 'r')
    result = []
    i = field_names().index(field)
    for track in h5f['tracks']:
        primitive = h5f['tracks'][track]['primitive'][()]
        result += primitive[:,i].tolist()
    return np.array(result)




def get_track_edges(fname):
    edges = []
    h5f = h5py.File(fname, 'r')
    tracks = list(h5f['tracks'])
    for track in tracks:
        face_radii = h5f['tracks'][track]['face_radii'][()]
        theta0     = h5f['tracks'][track]['theta0'][()]
        theta1     = h5f['tracks'][track]['theta1'][()]
        x0 = np.sin(theta0) * face_radii[0]
        y0 = np.cos(theta0) * face_radii[0]
        x1 = np.sin(theta0) * face_radii[-1]
        y1 = np.cos(theta0) * face_radii[-1]
        edges.append([(x0, y0), (x1, y1)])
        if track == tracks[-1]:
            x0 = np.sin(theta1) * face_radii[0]
            y0 = np.cos(theta1) * face_radii[0]
            x1 = np.sin(theta1) * face_radii[-1]
            y1 = np.cos(theta1) * face_radii[-1]
            edges.append([(x0, y0), (x1, y1)])
    return np.array(edges)




def get_radial_faces(fname):
    h5f = h5py.File(fname, 'r')
    faces = []
    for track in h5f['tracks']:
        face_radii = h5f['tracks'][track]['face_radii'][()]
        theta0     = h5f['tracks'][track]['theta0'][()]
        theta1     = h5f['tracks'][track]['theta1'][()]
        x0 = np.sin(theta0) * face_radii
        y0 = np.cos(theta0) * face_radii
        x1 = np.sin(theta1) * face_radii
        y1 = np.cos(theta1) * face_radii
        faces += np.array([[x0, x1], [y0, y1]]).T.tolist()
    return np.array(faces)




def get_max_radius(fname):
    h5f = h5py.File(fname, 'r')
    return max([h5f['tracks'][track]['face_radii'][-1] for track in h5f['tracks']])




def plot_field(ax, fname, field, triangulation=None, transform=None, draw_edges=False, **kwargs):
    ax.set_aspect('equal')

    if draw_edges:
        track_edges  = LineCollection(get_track_edges(fname), linewidth=0.5)
        radial_faces = LineCollection(get_radial_faces(fname), linewidth=0.5)
        ax.add_collection(track_edges)
        ax.add_collection(radial_faces)

    t = triangulation or cell_triangulation(fname)
    z = get_primitive_for_triangles(fname, field)
    return ax.tripcolor(t, transform(z) if transform else z, **kwargs)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('--width', type=float, default=10.0)
    parser.add_argument('-c', '--print-run-config', action='store_true')
    parser.add_argument('-e', '--draw-edges', action='store_true')

    args = parser.parse_args()

    fname = args.filename
    h5f = h5py.File(fname, 'r')

    if args.print_run_config:
        for key, val in h5f['run_config'].items():
            print(f'{key} = {val[()]}')

    width = 12
    fig = plt.figure(figsize=[width, width * 8 / 7])
    spec = fig.add_gridspec(2, 2, height_ratios=[29, 1])
    ax1  = fig.add_subplot(spec[0, 0])
    ax2  = fig.add_subplot(spec[0, 1])
    cax1 = fig.add_subplot(spec[1, 0])
    cax2 = fig.add_subplot(spec[1, 1])
    triangulation = cell_triangulation(fname)
    cbar1 = plot_field(ax1, fname, 'rho', triangulation, draw_edges=args.draw_edges, cmap='inferno', transform=np.log10)
    cbar2 = plot_field(ax2, fname, 'ur',  triangulation, draw_edges=args.draw_edges, cmap='viridis')#, vmin=0, vmax=10)
    plt.colorbar(cbar1, cax=cax1, orientation='horizontal')
    plt.colorbar(cbar2, cax=cax2, orientation='horizontal')
    ax1.set_title(r'$\log_{10} \rho$')
    ax2.set_title(r'Radial $\gamma \beta$')
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    fig.suptitle(r'$t = {:0.2f}$'.format(h5f['time'][()]))
    fig.tight_layout()
    plt.show()
