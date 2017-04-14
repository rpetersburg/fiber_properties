# Fiber Properties

fiber_properties is a python package that contains functions to calculate the focal ratio degradation, modal noise, and scrambling gain using simultaneous images from the ends of an optical fiber and the far field of the fiber output. It also contains example scripts that utilize these functions and formulate the appropriate data representations. The code was developed for use with the Fiber Characterization Station within the Yale Exoplanet Lab in New Haven, CT, but it may be used with other fiber imaging systems with some tweaking. Otherwise, the image calibration and fiber property measurement methods and algorithms definitely have applications outside of the FCS.

## Installation

Download the package into any available folder and then run in your terminal/command line: `python setup.py install`
If you don't have administrative rights on the computer, install only for a single user: `python setup.py install --user`

## Usage

There are many modules that can be imported into scripts:

* fiber_image.py: Contains the FiberImage class which is instantiated with a collection of images and a calibration object. FiberImage is a subclass of the CalibratedImage (calibrated_image.py) and BaseImage (base_image.py) classes that contain methods to instantiate the image and place the image header values into object attributes for easy access. FiberImage contains methods to find the fiber face's center, centroid, and radius as well as fiber properties including modal noise and focal ratio degradation.

* scrambling_gain.py: Contains function that calculates the scrambling gain given a list of fiber input and fiber output images

* plotting.py: Contains functions that plot fiber images and the results from the fiber property measurements (e.g. scrambling gain, frd, modal noise).
  
* numpy_array_handler.py: Contains functions that handle two dimensional numpy.ndarray objects that represent the fiber and calibration images. Functions include image cropping, function fitting (polynomial, gaussian), fft window application, and generic image array creation (tophat, rectangle, gaussian).

* input_output.py: Contains functions that save and load FiberImage objects and the related data as well as image_list() that produces quick file lists for easy FiberImage instantiation.
  
Basic functionality includes importing the FiberImage class, instantiating an objects with fiber image and calibration image file locations, and then calling the respective getter from the object. For example:

```python
from fiber_properties import FiberImage, image_list
image = image_list('image_file_')
dark = ['dark_file_1.fit', 'dark_file_2.fit']
ambient = ['ambient_file_1.fit', 'ambient_file_2.fit', 'ambient_file_3.fit']
image_obj = FiberImage(image, dark=dark, ambient=ambient)
y0, x0 = image_obj.get_fiber_center(method='radius')
diameter = image_obj.get_fiber_diameter(method='radius')
modal_noise = image_obj.get_modal_noise(method='filter')
frd_info = image_obj.get_frd_info()
```

## Contributing

1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request

## Credits

The fiber_properties package was written by Ryan Petersburg (ryan.petersburg@yale.edu) with contributions from Saki Kamon (University of Texas) and Tyler McCracken (Yale University). Please contact Ryan if you have any interest in the code or would like to make a contribution to the project.
