import ImageAnalysis

calibration_folder = "Calibration/TIF/"
image_folder = '../Alignment Images/2016-07-12/'

nf_images = []
for i in xrange(10):
    #nf_images.append(image_folder + 'nf_laser_noagitation_' + str(i) + '_1.8ms.tif')
    nf_images.append(image_folder + 'nf_laser_agitation_' + str(i) + '_1.8ms.tif') 
nf_dark_images = []
nf_ambient_images = []
for i in xrange(3):
    nf_dark_images.append(calibration_folder + 'nf_dark_' + str(i) + '_10ms.tif')
    nf_ambient_images.append(image_folder + 'nf_ambient_' + str(i) + '_1.8ms.tif')

nf_flat_images = []
#for i in xrange(8):
#    nf_flat_images.append(calibration_folder + 'nf_flat_' + str(i) + '_1ms.tif')

imAnalysis = ImageAnalysis.ImageAnalysis(nf_images, nf_dark_images,
                                         nf_flat_images, nf_ambient_images,
                                         pixel_size=3.45,
                                         magnification=10,
                                         bits_adc=16)
tol = .2
test_range = 10

print 'Height:', imAnalysis.getImageHeight(), 'Width:', imAnalysis.getImageWidth()
imAnalysis.showImageArray()
print
print 'Centroid'
centroid_row, centroid_column = imAnalysis.getFiberCentroid()
print 'Centroid Row:', centroid_row, 'Centroid Column:', centroid_column
print
print 'Edge:'
center_y, center_x = imAnalysis.getFiberCenter(method='edge')
print 'Diameter:', imAnalysis.getFiberDiameterMicrons(method='edge'), 'microns'
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Circle:'
center_y, center_x = imAnalysis.getFiberCenterCircleMethod(tol=tol, test_range=test_range)
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Radius:'
center_y, center_x = imAnalysis.getFiberCenterRadiusMethod(tol=tol, test_range=test_range)
print 'Diameter:', imAnalysis.getFiberDiameterMicrons(method='radius'), 'microns'
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Gaussian:'
center_y, center_x = imAnalysis.getFiberCenter(method='gaussian')
print 'Diameter:', imAnalysis.getFiberDiameterMicrons(method='gaussian'), 'microns'
print 'Center Row:', center_y, 'Center Column:', center_x
