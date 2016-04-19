import ImageAnalysis

imAnalysis = ImageAnalysis.ImageAnalysis("../Fiber Images/2016-03-15/NF_67_test1_3.5ms.tif")

imAnalysis.setDarkImage("../Fiber Images/2016-03-15/NF_67_dark_99ms.tif")
imAnalysis.setFlatFieldImage("../Fiber Images/2016-03-16/sn67_flatfield1.tif",
                                "../Fiber Images/2016-03-16/sn67_flatfield2.tif",
                                "../Fiber Images/2016-03-16/sn67_flatfield3.tif")
imAnalysis.executeErrorCorrections()

print 'Height:', imAnalysis.height, 'Width:', imAnalysis.width
imAnalysis.showImageArray()

centroidRow, centroidColumn = imAnalysis.getFiberCentroid()
print 'Centroid Row:', centroidRow, 'Centroid Column:', centroidColumn

centerRow, centerColumn = imAnalysis.getFiberCenterEdgeMethod()
print 'Edge: Center Row =', centerRow, 'Center Column =', centerColumn
print 'Fiber Width:', imAnalysis.getFiberDiameter()

centerY, centerX, minSum = imAnalysis.getFiberCenterCircleIteration(381.5)
print 'Circle: Center Row =', centerY, 'Center Column =', centerX

radius, centerY, centerX, minSum = imAnalysis.getFiberCenterRadiusIteration()
print 'Radius =', radius, 'Center Row =', centerY, 'Center Column =', centerX
