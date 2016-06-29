import ImageAnalysis
import numpy as np
import time

time.clock()

imAnalysis = ImageAnalysis.ImageAnalysis("../Fiber Images/2016-03-15/NF_67_test1_3.5ms.tif")

radius = 300
x = 500.5
y = 600.5
circle = imAnalysis.circleArray(radius, x, y, 16)
print 'Array area:', np.sum(circle)

area = np.pi * radius**2
print 'Actual area:', area

imAnalysis.showImageArray(circle)