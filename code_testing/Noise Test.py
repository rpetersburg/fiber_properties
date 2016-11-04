from FiberProperties import ImageAnalysis, Calibration

folder = 'FRD Measurements/Core Extension/2016-08-04 Reference Octagonal/'
calibration = Calibration([folder + 'Dark/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)],
                          None,
                          [folder + 'Ambient/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)])

images = [folder + 'Output 3.5/im_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
im_obj = ImageAnalysis(images, calibration, magnification=1)
im_obj.plotCrossSections()