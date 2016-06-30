import FiberProperties as FP
import ImageAnalysis as IA

laser_image = "../Alignment Images/2016-06-22 FCS First Proper Alignment/Laser far field.tif"

nf_dark_images = ["../Alignment Images/2016-03-15/NF_67_dark_99ms.tif"]
nf_flatfield_images = ["../Alignment Images/2016-03-16/sn67_flatfield1.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield2.tif",
                       "../Alignment Images/2016-03-16/sn67_flatfield3.tif"]

laser_image_analysis = IA.ImageAnalysis(laser_image, nf_dark_images, nf_flatfield_images)

fpObject = FP.FiberProperties(None, laser_image_analysis, None)

fpObject.showImageArray(fpObject.gaussianArray(laser_image_analysis))