import ImageAnalysis

ff_image = "../Alignment Images/2016-06-22 FCS First Proper Alignment/Laser far field.tif"
ff_dark_images = ["../Alignment Images/2016-03-15/ff_19_dark.tif"]
ff_flatfield_images = ["../Alignment Images/2016-03-16/sn67_flatfield1.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield2.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield3.tif"]

imAnalysis = ImageAnalysis.ImageAnalysis(ff_image, ff_dark_images, ff_flatfield_images)

print 'Height:', imAnalysis.getImageHeight(), 'Width:', imAnalysis.getImageWidth()
imAnalysis.showImageArray()

centroidRow, centroidColumn = imAnalysis.getFiberCentroid()
print 'Centroid Row:', centroidRow, 'Centroid Column:', centroidColumn

centerRow, centerColumn = imAnalysis.getFiberCenterEdgeMethod()
print 'Edge: Center Row =', centerRow, 'Center Column =', centerColumn
print 'Fiber Width:', imAnalysis.getFiberDiameter()
imAnalysis.plotHorizontalCrossSection(imAnalysis.getImageArray(), centerRow)

centerY, centerX = imAnalysis.getFiberCenterCircleIteration(381.5)
print 'Circle: Center Row =', centerY, 'Center Column =', centerX

centerY, centerX = imAnalysis.getFiberCenterRadiusIteration()
print 'Diameter =', imAnalysis.getFiberDiameter, 'Center Row =', centerY, 'Center Column =', centerX
