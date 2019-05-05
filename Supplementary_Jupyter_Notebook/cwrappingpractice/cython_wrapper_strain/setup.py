from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

examples_extension = Extension(
    name="TDStrain",
    sources=["TDStrain.pyx"],
    libraries=["examples"],
    library_dirs=["lib"],
    include_dirs=["lib"]
)
setup(
    name="TDStrain",
    ext_modules=cythonize([examples_extension])
)
