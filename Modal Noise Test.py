import FiberProperties as FP
import ImageAnalysis as IA

LED_image = "../Alignment Images/2016-06-22 FCS First Proper Alignment/LED near field.tif"
laser_image = "../Alignment Images/2016-06-22 FCS First Proper Alignment/Laser near field.tif"
nf_dark_images = ["../Alignment Images/2016-03-15/NF_67_dark_99ms.tif"]
nf_flatfield_images = ["../Alignment Images/2016-03-16/sn67_flatfield1.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield2.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield3.tif"]


LED_image_analysis = IA.ImageAnalysis(LED_image, nf_dark_images, nf_flatfield_images)
laser_image_analysis = IA.ImageAnalysis(laser_image, nf_dark_images, nf_flatfield_images)

fpObject = FP.FiberProperties(None, LED_image_analysis, None)
print "LED:", fpObject.getModalNoise()

fpObject = FP.FiberProperties(None, laser_image_analysis, None)
print "Laser:", fpObject.getModalNoise()