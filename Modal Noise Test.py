import FiberProperties as FP
import ImageAnalysis as IA

LED_image = "../Alignment Images/2016-06-28/Near field normal camera distance zoomed out.tif"
LED_image2 = "../Alignment Images/2016-03-15/NF_67_fiber.tif"
#LED_image = "../Alignment Images/2016-06-22 FCS First Proper Alignment/LED near field.tif"
laser_image = "../Alignment Images/2016-06-22 FCS First Proper Alignment/Laser near field.tif"
#laser_image = "../Alignment Images/2016-03-15/NF_67_test1_3.5ms.tif"
nf_dark_images = ["../Alignment Images/2016-03-15/NF_67_dark_99ms.tif"]
nf_flatfield_images = ["../Alignment Images/2016-03-16/sn67_flatfield1.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield2.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield3.tif"]


LED_image_analysis = IA.ImageAnalysis(LED_image, nf_dark_images, nf_flatfield_images)
LED_image_analysis2 = IA.ImageAnalysis(LED_image2, nf_dark_images, nf_flatfield_images)
laser_image_analysis = IA.ImageAnalysis(laser_image, nf_dark_images, nf_flatfield_images)

fpObject = FP.FiberProperties(None, LED_image_analysis, None)
print "LED:"
for i in xrange(10):
    i = 0.1 + 0.1 * i
    test = fpObject.getModalNoise(i)
    print "Factor:", i, "MN:", test
print

fpObject = FP.FiberProperties(None, LED_image_analysis2, None)
print "LED:"
for i in xrange(10):
    i = 0.1 + 0.1 * i
    test = fpObject.getModalNoise(i)
    print "Factor:", i, "MN:", test

print

fpObject = FP.FiberProperties(None, laser_image_analysis, None)
print "Laser:"
for i in xrange(10):
    i = 0.1 + 0.1 * i
    test = fpObject.getModalNoise(i)
    print "Factor:", i, "MN:", test