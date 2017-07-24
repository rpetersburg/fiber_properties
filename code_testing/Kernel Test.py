from FiberProperties import ImageAnalysis, Calibration

base_folder = 'Stability Measurements/2016-08-15 Stability Test Unagitated/'
ambient_folder = base_folder + 'Ambient/'
dark_folder = base_folder + 'Dark/'
flat_folder = base_folder + 'Flat/'
folder = base_folder + 'Images/'
ext = '.fit'

in_calibration = Calibration(dark=[dark_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(10)],
                             ambient=[ambient_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(10)],
                             flat=[flat_folder + 'in_' + str(i).zfill(3) + ext for i in xrange(8)])

for num_images in [1, 10, 100]:
    for kernel in [1, 9, 19]:
        in_images = [folder + 'in_' + str(i).zfill(3) + ext for i in xrange(num_images)]
        im_obj = ImageAnalysis(in_images, in_calibration, kernel_size=kernel)
        print 'num_images:', num_images, 'kernel:', kernel
        print 'diameter:', im_obj.getFiberDiameter(method='edge')
        print 'center:', im_obj.getFiberCenter(method='edge')
        print
        # im_obj.saveData('Kernel Test/', str(num_images) + '_im_' + str(kernel) + '_kern')