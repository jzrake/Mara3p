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




import h5py




def block_index_from_string(index_string):
    level, rest = index_string.split(':')
    return [int(level)] + [int(i) for i in rest.split('-')]




def block_extent_1d(level, i, domain_size=1.0):
    i0 = i
    i1 = i + 1
    dl = float(1 << level)
    x0 = domain_size * (-0.5 + 1.0 * i0 / dl)
    x1 = domain_size * (-0.5 + 1.0 * i1 / dl)
    return x0, x1




def block_extent_2d(level, i, j, domain_size=1.0):
    i0, j0 = i, j
    i1, j1 = i + 1, j + 1
    dl = float(1 << level)
    x0 = domain_size * (-0.5 + 1.0 * i0 / dl)
    x1 = domain_size * (-0.5 + 1.0 * i1 / dl)
    y0 = domain_size * (-0.5 + 1.0 * j0 / dl)
    y1 = domain_size * (-0.5 + 1.0 * j1 / dl)
    return x0, x1, y0, y1




def block_edges_1d(block_list):
    edges = set()
    for block in block_list:
        x0, x1 = block_extent_1d(*block)
        edges.add(x0)
        edges.add(x1)
    return edges




def block_edges_2d(block_list):
    edges = set()
    for block in block_list:
        x0, x1, y0, y1 = block_extent_2d(*block)
        edges.add((x0, y0, x0, y1))
        edges.add((x0, y0, x1, y0))
        edges.add((x1, y0, x1, y1))
        edges.add((x0, y1, x1, y1))
    return edges




def read_block_edges_1d(filename):
    return block_edges_1d([block_index_from_string(block) for block in h5py.File(filename, 'r')['solution']['conserved']])




def read_block_edges_2d(filename):
    return block_edges_2d([block_index_from_string(block) for block in h5py.File(filename, 'r')['solution']['conserved']])




def upsample(data, fold):
    if fold == 0:
        return data
    nx, ny = data.shape
    result = np.zeros([nx * 2, ny * 2])
    result[0:nx*2:2, 0:ny*2:2] = data
    result[1:nx*2:2, 0:ny*2:2] = data
    result[0:nx*2:2, 1:ny*2:2] = data
    result[1:nx*2:2, 1:ny*2:2] = data
    return upsample(result, fold - 1)

