from distutils.core import setup

# Place this `setup.py` at the root of the project to build and
# distribute a pure-Python version of SGP4.

# Read the package's docstring and "__version__" without importing it.
path = 'sgp4/__init__.py'
with open(path, 'rb') as f:
    text = f.read().decode('utf-8')
text = text.replace('-*- coding: utf-8 -*-', '')  # for Python 2.7

namespace = {}
eval(compile(text, path, 'exec'), namespace)

version = namespace['__version__']
description, long_description = namespace['__doc__'].split('\n', 1)

# Avoid "`long_description` has syntax errors in markup" from extra
# whitespace that somehow creeps in.
description = description.strip()
long_description = """\
This distribution provides the `sgp4 package <https://pypi.org/project/sgp4/>`_
for computing Earth satellite positions, but without making any attempt to
compile and install the accelerated C++ version of the algorithm.  Instead,
it is implemented in (rather slow) pure Python code.
"""

setup(
    name = 'sgp4_pure_python',
    version = version,
    description = description,
    long_description = long_description,
    long_description_content_type = 'text/x-rst',
    license = 'MIT',
    author = 'Brandon Rhodes',
    author_email = 'brandon@rhodesmill.org',
    url = 'https://github.com/brandon-rhodes/python-sgp4',
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    packages = ['sgp4'],
    package_data = {'sgp4': ['SGP4-VER.TLE', 'sample*', 'tcppver.out']},
    provides = ['sgp4'],
)
