from FiberProperties import ImageAnalysis

#nf_images = ["../Alignment Images/2016-06-30/nf_sn67_0.30ms.tif",
#             "../Alignment Images/2016-06-30/nf_sn67_0.09ms.tif",
#             "../Alignment Images/2016-06-30/nf_sn67_2.00ms.tif",
#             "../Alignment Images/2016-06-30/nf_sn67_5.00ms.tif",
#             "../Alignment Images/2016-06-30/nf_sn67_10.0ms.tif"]

nf_images = ["../Alignment Images/2016-06-28/Near field normal camera distance zoomed out.tif"]

nf_dark_images = ["../Alignment Images/2016-06-30/dark_sn67_100.0ms.tif"]
nf_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn67_0.98ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_3.tif"]

for nf_image in nf_images:
    print nf_image

    nf_image_analysis = ImageAnalysis(nf_image, nf_dark_images, nf_flatfield_images)

    center_y, center_x = nf_image_analysis.getFiberCenter()
    diameter = nf_image_analysis.getFiberDiameter()
    print "Edge:"
    print "Center row:", center_y, "Center column:", center_x
    print "Diameter:", diameter

    center_y, center_x = nf_image_analysis.getFiberCenterTophatMethod()
    diameter = nf_image_analysis.getFiberDiameter()
    print "Tophat:"
    print "Center row:", center_y, "Center column:", center_x
    print "Diameter:", diameter