from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

examples_extension = Extension(
    name="pyDispInHalfSpace",
    sources=["pyDispInHalfSpace.pyx"],
    libraries=["DispInHalfSpace"],
    library_dirs=["lib"],
    include_dirs=["lib"]
)
setup(
    name="pyDispInHalfSpace",
    ext_modules=cythonize([examples_extension])
)
