# Fiber Properties

Fiber Properties is a python package that contains functions to calculate the focal ratio degradation, modal noise, and scrambling gain using simultaneous images from the ends of an optical fiber and the far field of the fiber output. It also contains example scripts that utilize these functions and formulate the appropriate data representations. The code was developed for use with the Fiber Characterization Station within the Yale Exoplanet Lab in New Haven, CT, but it may be used with other fiber imaging systems with some tweaking. Otherwise, the image calibration and fiber property measurement methods and algorithms definitely have applications outside of the FCS.

## Installation

Simply place the python packages in a file location which can be imported from your own script.

## Usage

There are four primary modules that can be imported into scripts:

* ImageAnalysis.py: Contains the ImageAnalysis class which is instantiated with a collection of images and a calibration object. The object then contains methods to find the fiber face's centroid and center as well as header information from the images.
  
* Calibration.py: Contains the Calibration class which is instantiated with any combination of dark, flat field, and ambient images and executes corrections to other images based on the calibration images.
  
* FiberProperties.py: Contains functions FRD(), scramblingGain(), and modalNoise() which calculate the respective fiber properties when given a single ImageAnalysis object (for FRD and Modal Noise) or a collection of ImageAnalysis objects (for Scrambling Gain)
  
* NumpyArrayHandler.py: Contains functions that handle two dimensional numpy.ndarray objects that represent the fiber and calibration images. Functions include converting image files to numpy.ndarrays, image cropping, function fitting, and image plotting.
  
Basic functionality includes importing the ImageAnlysis and Calibration class, instantiating the objects with image file locations, and then executing the relevant FiberProperties function on the object. For example:

`from ImageAnalysis import ImageAnalysis`

`from Calibration import Calibration`

`from FiberProperties import FRD`

`calibration = Calibration(dark=['dark_file_1.fit', dark_file_2.fit'], flat='flat_file.fit', ambient=None)`

`image_obj = ImageAnalysis(['image_file_' + str(i) + '.fit' for i in xrange(10)], calibration)`

`frd_output = FRD(image_obj)`

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request

## Credits

The Fiber Properties package was written by Ryan Petersburg (ryan.petersburg@yale.edu) with contributions from Saki Kamon (University of Texas) and Tyler McCracken (Yale University). Please contact Ryan if you have any interest in the code or would like to make a contribution to the project.
