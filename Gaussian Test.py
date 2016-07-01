import FiberProperties as FP
import ImageAnalysis as IA

ffr_image = "../Alignment Images/2016-06-30/FF_sn19_1.00ms.tif"

ff_dark_images = ["../Alignment Images/2016-06-30/dark_sn19_100.0ms.tif"]
ff_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn19_1.00ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_3.tif"]

ff_image_analysis = IA.ImageAnalysis(ff_image, ff_dark_images, ff_flatfield_images)

center_y, center_x = ff_image_analysis.getFiberCenterGaussianMethod()

