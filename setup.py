from distutils.core import setup
from distutils.extension import Extension

USE_CYTHON = True
try:
    from Cython.Distutils import build_ext
except ImportError:
    USE_CYTHON = False

cmdclass = {}
ext_modules = []

if USE_CYTHON:
    ext_modules += [Extension('fiber_properties.cfilter_image', ['fiber_properties/cfilter_image.pyx'])]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [Extension('fiber_properties.cfilter_image', ['fiber_properties/cfilter_image.c'])]

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