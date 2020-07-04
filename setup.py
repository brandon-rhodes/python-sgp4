import os
import sys
import platform
from distutils.core import setup, Extension
from textwrap import dedent

import sgp4, sgp4.model

description, long_description = sgp4.__doc__.split('\n', 1)
satdoc = dedent(sgp4.model.Satellite.__doc__.split('\n', 1)[1])
long_description = long_description.replace('entry.', 'entry.' + satdoc)

# Force compilation on Travis CI + Python 3 to make sure it keeps working.
optional = True
if sys.version_info[0] != 2 and os.environ.get('TRAVIS') == 'true':
    optional = False

# It is hard to write C extensions that support both Python 2 and 3, so
# we opt here to support the acceleration only for Python 3.


def get_compile_args():

    comp_args = {
        'extra_compile_args': ['-fopenmp', '-O3'],  # , '-std=c++11'
        'extra_link_args': ['-fopenmp'],
        'language': 'c++',
        }

    if platform.system().lower() == 'windows':
        comp_args['extra_compile_args'] = ['/openmp']
        del comp_args['extra_link_args']

    elif 'darwin' in platform.system().lower():

        extra_compile_args = [
            '-fopenmp', '-O3',  # '-std=c++11',
            '-mmacosx-version-min=10.7',
            ]

        # if c++11 stdlib support is necessary and clang compiler is used
        # from subprocess import getoutput
        # if 'clang' in getoutput('gcc -v'):
        #     extra_compile_args += ['-stdlib=libc++', ]

        comp_args['extra_compile_args'] = extra_compile_args

    return comp_args


ext_modules = []
if sys.version_info[0] == 3:
    ext_modules.append(Extension(
        'sgp4.vallado_cpp',
        optional=optional,
        sources=[
            'extension/SGP4.cpp',
            'extension/wrapper.cpp',
        ],
        **get_compile_args(),
        # TODO: can we safely figure out how to use a pair of options
        # like these, adapted to as many platforms as possible, to use
        # multiple processors when available?
        # extra_compile_args=['-fopenmp'],
        # extra_link_args=['-fopenmp'],
    ))

setup(name = 'sgp4',
      version = '2.12',
      description = description,
      long_description = long_description,
      license = 'MIT',
      author = 'Brandon Rhodes',
      author_email = 'brandon@rhodesmill.org',
      url = 'https://github.com/brandon-rhodes/python-sgp4',
      classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
      packages = ['sgp4'],
      package_data = {'sgp4': ['SGP4-VER.TLE', 'tcppver.out']},
      ext_modules = ext_modules,
)
