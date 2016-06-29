import ImageAnalysis

nf_dark_images = ["../Alignment Images/2016-03-15/NF_67_dark_99ms.tif"]
nf_flatfield_images = ["../Alignment Images/2016-03-16/sn67_flatfield1.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield2.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield3.tif"]

imAnalysis = ImageAnalysis.ImageAnalysis("../Alignment Images/2016-03-15/NF_67_test1_3.5ms.tif",
                                         nf_dark_images, nf_flatfield_images)

imAnalysis.executeErrorCorrections()

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
