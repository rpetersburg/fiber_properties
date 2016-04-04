import ImageAnalysis

imAnalysis = ImageAnalysis.ImageAnalysis("../Fiber Images/2016-03-15/NF_67_test1_3.5ms.tif")

imAnalysis.setDarkImage("../Fiber Images/2016-03-15/NF_67_dark_99ms.tif")
imAnalysis.setFlatFieldImage("../Fiber Images/2016-03-16/sn67_flatfield1.tif",
                                "../Fiber Images/2016-03-16/sn67_flatfield2.tif",
                                "../Fiber Images/2016-03-16/sn67_flatfield3.tif")
imAnalysis.executeErrorCorrections()

centroidRow, centroidColumn = imAnalysis.getFiberCentroid()
print 'Centroid Row:', centroidRow, 'Centroid Column:', centroidColumn

centerRow, centerColumn = imAnalysis.getFiberCenterEdgeMethod()
print 'Edge: Center Row =', centerRow, 'Center Column =', centerColumn
print 'Fiber Width:', imAnalysis.getFiberDiameter()

centerY, centerX = imAnalysis.getFiberCenterCircleIteration(imAnalysis.getFiberDiameter() / 2)
print 'Circle: Center Row =', centerY, 'Center Column =', centerX

print
#==============================================================================================

imAnalysis = ImageAnalysis.ImageAnalysis("../Fiber Images/2016-03-15/NF_67_test2_3.5ms.tif")

imAnalysis.setDarkImage("../Fiber Images/2016-03-15/NF_67_dark_99ms.tif")
imAnalysis.setFlatFieldImage("../Fiber Images/2016-03-16/sn67_flatfield1.tif",
                                "../Fiber Images/2016-03-16/sn67_flatfield2.tif",
                                "../Fiber Images/2016-03-16/sn67_flatfield3.tif")
imAnalysis.executeErrorCorrections()

centroidRow, centroidColumn = imAnalysis.getFiberCentroid()
print 'Centroid Row:', centroidRow, 'Centroid Column:', centroidColumn

centerRow, centerColumn = imAnalysis.getFiberCenterEdgeMethod()
print 'Edge: Center Row =', centerRow, 'Center Column =', centerColumn
print 'Fiber Width:', imAnalysis.getFiberDiameter()

centerY, centerX = imAnalysis.getFiberCenterCircleIteration(imAnalysis.getFiberDiameter() / 2)
print 'Circle: Center Row =', centerY, 'Center Column =', centerX

print
#=============================================================================================

imAnalysis = ImageAnalysis.ImageAnalysis("../Fiber Images/2016-03-15/NF_67_test3_3.5ms.tif")

imAnalysis.setDarkImage("../Fiber Images/2016-03-15/NF_67_dark_99ms.tif")
imAnalysis.setFlatFieldImage("../Fiber Images/2016-03-16/sn67_flatfield1.tif",
                                "../Fiber Images/2016-03-16/sn67_flatfield2.tif",
                                "../Fiber Images/2016-03-16/sn67_flatfield3.tif")
imAnalysis.executeErrorCorrections()

print 'Height:', imAnalysis.height, 'Width:', imAnalysis.width

centroidRow, centroidColumn = imAnalysis.getFiberCentroid()
print 'Centroid Row:', centroidRow, 'Centroid Column:', centroidColumn

centerRow, centerColumn = imAnalysis.getFiberCenterEdgeMethod()
print 'Edge: Center Row =', centerRow, 'Center Column =', centerColumn
print 'Fiber Width:', imAnalysis.getFiberDiameter()

centerY, centerX = imAnalysis.getFiberCenterCircleIteration(imAnalysis.getFiberDiameter() / 2)
print 'Circle: Center Row =', centerY, 'Center Column =', centerX