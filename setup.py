from setuptools import setup
from Cython.Build import cythonize

setup(
        ext_modules = cythonize("fillH.pyx",annotate=True,compiler_directives={'language_level' : "3",'boundscheck':False,'wraparound':False})
)
setup(
        ext_modules = cythonize("fillH_largeN.pyx",annotate=True,compiler_directives={'language_level' : "3",'boundscheck':False,'wraparound':False})
)
setup(
        ext_modules = cythonize("expand_sfac.pyx",annotate=True,compiler_directives={'language_level' : "3",'boundscheck':False,'wraparound':False,'cdivision':True})
)

