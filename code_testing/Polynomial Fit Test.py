from FiberProperties import ImageAnalysis

ff_image = "../Alignment Images/2016-06-30/FF_sn19_1.00ms.tif"
ff_dark_images = ["../Alignment Images/2016-06-30/dark_sn19_100.0ms.tif"]
ff_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn19_1.00ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn19_1.00ms_3.tif"]

nf_image = "../Alignment Images/2016-03-15/NF_67_fiber.tif"
nf_dark_images = ["../Alignment Images/2016-06-30/dark_sn67_100.0ms.tif"]
nf_flatfield_images = ["../Alignment Images/2016-06-30/flat_sn67_0.98ms_1.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_2.tif",
                       "../Alignment Images/2016-06-30/flat_sn67_0.98ms_3.tif"]

nf_image_analysis = ImageAnalysis(nf_image, nf_dark_images, nf_flatfield_images)
nf_image_analysis.showImageArray()
nf_image_analysis.polynomialFit()


"""
ff_image_analysis = IA.ImageAnalysis(ff_image, ff_dark_images, ff_flatfield_images)
ff_image_analysis.showImageArray()
ff_image_analysis.polynomialFit()
"""