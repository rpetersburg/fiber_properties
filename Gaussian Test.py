import FiberProperties as FP
import ImageAnalysis as IA

in_images = ["../Alignment Images/2016-06-30/in_sn86_0.17ms.tif",
             "../Alignment Images/2016-06-30/in_sn86_1.00ms.tif",
             "../Alignment Images/2016-06-30/in_sn86_5.00ms.tif",
             "../Alignment Images/2016-06-30/in_sn86_10.0ms.tif",
             "../Alignment Images/2016-06-30/in_sn86_20.0ms.tif"]

in_dark_images = ["../Alignment Images/2016-06-30/dark_sn86_100.0ms.tif"]
in_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn86_1.00ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn86_1.00ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn86_1.00ms_3.tif"]

for in_image in in_images:
    print in_image

    in_image_analysis = IA.ImageAnalysis(in_image, in_dark_images, in_flatfield_images)
    in_image_analysis.showImageArray()

    center_y, center_x = in_image_analysis.getFiberCenter()
    radius = in_image_analysis.getFiberRadius()
    print "Edge:"
    print "Center row:", center_y, "Center column:", center_x
    print "Radius:", radius

    center_y, center_x = in_image_analysis.getFiberCenterGaussianMethod()
    radius = in_image_analysis.getFiberRadius()
    print "Gaussian:"
    print "Center row:", center_y, "Center column:", center_x
    print "Radius:", radius