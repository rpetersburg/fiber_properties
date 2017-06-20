from FiberProperties import ImageAnalysis

imAnalysis = ImageAnalysis.ImageAnalysis("../Alignment Images/2016-06-28/Near field normal camera distance zoomed out.tif",
                                         ["../Alignment Images/2016-03-15/NF_67_dark_99ms.tif"],
                                         ["../Alignment Images/2016-03-16/sn67_flatfield1.tif",
                                          "../Alignment Images/2016-03-16/sn67_flatfield2.tif",
                                          "../Alignment Images/2016-03-16/sn67_flatfield3.tif"])
imAnalysis.executeErrorCorrections()

print 'Height:', imAnalysis.getImageHeight(), 'Width:', imAnalysis.getImageWidth()
imAnalysis.showImageArray()

centerRow, centerColumn = imAnalysis.getFiberCenterEdgeMethod()
print 'Edge: Center Row =', centerRow, 'Center Column =', centerColumn
print 'Fiber Width:', imAnalysis.setFiberDiameter()

#centerY, centerX, minSum = imAnalysis.getFiberCenterCircleIteration(106)
#print 'Circle: Center Row =', centerY, 'Center Column =', centerX

radius, centerY, centerX, minSum = imAnalysis.getFiberCenterRadiusIteration()
print 'Radius =', radius, 'Center Row =', centerY, 'Center Column =', centerX