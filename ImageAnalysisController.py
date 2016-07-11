import ImageAnalysis

"""
nf_image = "../Alignment Images/2016-06-30/nf_sn67_0.30ms.tif"
nf_dark_images = ["../Alignment Images/2016-06-30/dark_sn67_100.0ms.tif"]
nf_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn67_0.98ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_3.tif"]

imAnalysis = ImageAnalysis.ImageAnalysis(nf_image, nf_dark_images, nf_flatfield_images)

print 'Height:', imAnalysis.getImageHeight(), 'Width:', imAnalysis.getImageWidth()
imAnalysis.showImageArray()

centroid_row, centroid_column = imAnalysis.getFiberCentroid()
print 'Centroid Row:', centroid_row, 'Centroid Column:', centroid_column
print
print 'Edge:'
center_y, center_x = imAnalysis.getFiberCenterEdgeMethod()
print 'Diameter:', imAnalysis.getFiberDiameterEdgeMethod()
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Circle:'
center_y, center_x = imAnalysis.getFiberCenterCircleIteration(tol=0.1, test_range=20)
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Radius:'
center_y, center_x = imAnalysis.getFiberCenterRadiusIteration(tol=0.1, test_range=20)
print 'Diameter:', imAnalysis.getFiberDiameterRadiusMethod()
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Gaussian:'
center_y, center_x = imAnalysis.getFiberCenterGaussianMethod()
print 'Diameter:', imAnalysis.getFiberDiameterGaussianMethod()
print 'Center Row:', center_y, 'Center Column:', center_x
"""

ff_image =  "../Alignment Images/2016-06-22 FCS First Proper Alignment/LED far field.tif"
ff_dark_images = ["../Alignment Images/2016-06-30/dark_sn19_100.0ms.tif"]
ff_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn19_1.00ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_3.tif"]

imAnalysis = ImageAnalysis.ImageAnalysis(ff_image, ff_dark_images, ff_flatfield_images)
imAnalysis.showImageArray()

'Edge:'
center_y, center_x = imAnalysis.getFiberCenterEdgeMethod()
print 'Diameter:', imAnalysis.getFiberDiameterEdgeMethod()
print 'Center Row:', center_y, 'Center Column:', center_x

print 'Radius:'
center_y, center_x = imAnalysis.getFiberCenterRadiusIteration(tol=0.1, test_range=None)
print 'Diameter:', imAnalysis.getFiberDiameterRadiusMethod()
print 'Center Row:', center_y, 'Center Column:', center_x