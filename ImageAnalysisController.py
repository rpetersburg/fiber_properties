import ImageAnalysis

ff_image = "../Alignment Images/2016-06-30/FF_sn19_1.00ms.tif"
ff_dark_images = ["../Alignment Images/2016-06-30/dark_sn19_100.0ms.tif"]
ff_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn19_1.00ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_3.tif"]

imAnalysis = ImageAnalysis.ImageAnalysis(ff_image, ff_dark_images, ff_flatfield_images)

print 'Height:', imAnalysis.getImageHeight(), 'Width:', imAnalysis.getImageWidth()
imAnalysis.showImageArray()

centroidRow, centroidColumn = imAnalysis.getFiberCentroid()
print 'Centroid Row:', centroidRow, 'Centroid Column:', centroidColumn

centerRow, centerColumn = imAnalysis.getFiberCenterEdgeMethod()
print 'Edge: Center Row =', centerRow, 'Center Column =', centerColumn
print 'Fiber Width:', imAnalysis.getFiberDiameter()
imAnalysis.plotHorizontalCrossSection(imAnalysis.getImageArray(), centerRow)

centerY, centerX = imAnalysis.getFiberCenterCircleIteration()
print 'Circle: Center Row =', centerY, 'Center Column =', centerX

centerY, centerX = imAnalysis.getFiberCenterRadiusIteration()
print 'Diameter =', imAnalysis.getFiberDiameter(), 'Center Row =', centerY, 'Center Column =', centerX
