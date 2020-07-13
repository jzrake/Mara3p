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

import sys
import h5py




def run_module(module_name):
	if module_name == b'euler1d':
		from .modules import euler1d
		euler1d.main()
	elif module_name == b'euler2d':
		from .modules import euler2d
		euler2d.main()
	elif module_name == b'locally_isothermal':
		from .modules import locally_isothermal
		locally_isothermal.main()
	else:
		print('ploting support is not available for module', module_name)




module_names = set([h5py.File(f, 'r')['module'][()] for f in sys.argv if f.endswith('.h5')])

if not module_names:
	print('[mara.plot] no input files given')

elif len(module_names) > 1:
	print('[mara.plot] input files are from different modules:', module_names)

else:
	run_module([*module_names][0])
