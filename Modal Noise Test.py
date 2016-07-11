import FiberProperties as FP
import ImageAnalysis as IA

#LED_image = "../Alignment Images/2016-06-28/Near field normal camera distance zoomed out.tif"
#LED_image2 = "../Alignment Images/2016-03-15/NF_67_fiber.tif"
#laser_image = "../Alignment Images/2016-03-15/NF_67_test1_3.5ms.tif"
#LED_nf_file = "../Alignment Images/2016-06-28/Near field normal camera distance zoomed out.tif"
LED_nf_file = "../Alignment Images/2016-06-22 FCS First Proper Alignment/LED near field.tif"
LED_ff_file = "../Alignment Images/2016-06-22 FCS First Proper Alignment/LED far field.tif"

laser_nf_file = "../Alignment Images/2016-06-22 FCS First Proper Alignment/Laser near field.tif"
laser_ff_file =  "../Alignment Images/2016-06-30/FF_sn19_1.00ms.tif"

nf_dark_images = ["../Alignment Images/2016-06-30/dark_sn67_100.0ms.tif"]
nf_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn67_0.98ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_3.tif"]


ff_dark_images = ["../Alignment Images/2016-06-30/dark_sn19_100.0ms.tif"]
ff_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn19_1.00ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_3.tif"]

LED_nf = IA.ImageAnalysis(LED_nf_file, nf_dark_images, nf_flatfield_images)
LED_ff = IA.ImageAnalysis(LED_ff_file, ff_dark_images, ff_flatfield_images)
laser_nf = IA.ImageAnalysis(laser_nf_file, nf_dark_images, nf_flatfield_images)
laser_ff = IA.ImageAnalysis(laser_ff_file, ff_dark_images, ff_flatfield_images)

laser_fp = FP.FiberProperties(None, laser_nf, laser_ff)
LED_fp = FP.FiberProperties(None, LED_nf, LED_ff)

print "LED Tests:"
print "FF MN Poly:", LED_fp.getModalNoise('far', 'polynomial')
print "NF MN Grad:", LED_fp.getModalNoise('far', 'gradient')
print "NF MN Poly:", LED_fp.getModalNoise('near', 'polynomial')
print "NF MN TopH:", LED_fp.getModalNoise('near', 'tophat')
print "NF MN Grad:", LED_fp.getModalNoise('near', 'gradient')


print "Laser Tests:"
print "FF MN Poly:", laser_fp.getModalNoise('far', 'polynomial')
print "NF MN Grad:", laser_fp.getModalNoise('far', 'gradient')
print "NF MN Poly:", laser_fp.getModalNoise('near', 'polynomial')
print "NF MN TopH:", laser_fp.getModalNoise('near', 'tophat')
print "NF MN Grad:", laser_fp.getModalNoise('near', 'gradient')
