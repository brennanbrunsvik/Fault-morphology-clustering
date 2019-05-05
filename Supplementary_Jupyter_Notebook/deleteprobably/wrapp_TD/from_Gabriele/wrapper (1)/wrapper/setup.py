from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

wrapper_extension = Extension(
    name="pywrapper",
    sources=["pywrapper.pyx"],
    libraries=["wrapper"],
    library_dirs=["lib"],
    include_dirs=["lib"]
)
setup(
    name="pywrapper",
    ext_modules=cythonize([wrapper_extension])
)
