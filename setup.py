from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

# ext_modules = [Extension('fiber_properties.filter_image',
#                          ['fiber_properties/filter_image.pyx'],
#                          include_dirs=[np.get_include()])]
ext_modules = []
cmdclass = {'build_ext': build_ext}

setup(name='fiber_properties',
      version='1.0',
      author='Ryan Petersburg',
      author_email='ryan.petersburg@yale.edu',
      description='Image Analysis for Fiber Characterization Station',
      license='MIT',
      url='https://github.com/rpetersburg/fiber_properties',
      packages=['fiber_properties'],
      cmdclass=cmdclass,
      ext_modules=ext_modules
     )