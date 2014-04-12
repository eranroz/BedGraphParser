from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "BedGraphReader",
    ext_modules = cythonize('*.pyx'),
	author='Eranroz',
	url='https://github.com/eranroz/BedGraphParser',
)

