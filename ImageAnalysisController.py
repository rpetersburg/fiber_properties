import ImageAnalysis

nf_image = "../Alignment Images/2016-06-22 FCS First Proper Alignment/LED near field.tif"
nf_dark_images = ["../Alignment Images/2016-06-30/dark_sn67_100.0ms.tif"]
nf_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn67_0.98ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_3.tif"]

imAnalysis = ImageAnalysis.ImageAnalysis(nf_image, nf_dark_images, nf_flatfield_images)

print 'Height:', imAnalysis.getImageHeight(), 'Width:', imAnalysis.getImageWidth()
#imAnalysis.showImageArray()

centroid_row, centroid_column = imAnalysis.getFiberCentroid()
print 'Centroid Row:', centroid_row, 'Centroid Column:', centroid_column
print
print 'Edge:'
center_y, center_x = imAnalysis.getFiberCenterEdgeMethod()
print 'Diameter:', imAnalysis.getFiberDiameterEdgeMethod()
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Circle:'
center_y, center_x = imAnalysis.getFiberCenterCircleIteration()
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Radius:'
center_y, center_x = imAnalysis.getFiberCenterRadiusIteration()
print 'Diameter:', imAnalysis.getFiberDiameterRadiusMethod()
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Top Hat:'
center_y, center_x = imAnalysis.getFiberCenterTophatMethod()
print 'Diameter:', imAnalysis.getFiberDiameterTophatMethod()
print 'Center Row:', center_y, 'Center Column:', center_x
print
print 'Gaussian:'
center_y, center_x = imAnalysis.getFiberCenterGaussianMethod()
print 'Diameter:', imAnalysis.getFiberDiameterGaussianMethod()
print 'Center Row:', center_y, 'Center Column:', center_x
