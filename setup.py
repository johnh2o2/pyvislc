#!/usr/bin/env python

from distutils.core import setup

setup(name='pyvislc',
      version='1.0',
      description='Lightcurve visualizer using Tkinter',
      author='John Hoffman',
      author_email='jah5@princeton.edu',
      url='https://github.com/johnh2o2/pyvislc',
      package_dir = {'pyvislc': '.'},
      packages=['pyvislc'],
      requires=['matplotlib', 'shlex', 'pyfits', 'subprocess', 'numpy', 'scipy', 'pandas'],
      provides=['pyvislc']
      #py_modules = [ 'vislc', 'defaults', 'vlcutils']
      #package_data={'mypkg': ['data/*.dat']},
      #packages=['pyvislc' ],
     )