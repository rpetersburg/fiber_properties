import FiberProperties as FP
import ImageAnalysis as IA
import os as os

calibration_folder = "Calibration/TIF/"
image_folder = '../Alignment Images/2016-07-12/'

nf_led_images = []
for i in xrange(3):
    nf_led_images.append(image_folder + 'nf_led_0um_' + str(i) + '_80ms.tif')

nf_laser_images = []
nf_laser_images_agitated = []
for i in xrange(10):
    nf_laser_images.append(image_folder + 'nf_laser_noagitation_' + str(i) + '_1.8ms.tif')
    nf_laser_images_agitated.append(image_folder + 'nf_laser_agitation_' + str(i) + '_1.8ms.tif')    

nf_dark_images = []
nf_ambient_images = []
for i in xrange(3):
    nf_dark_images.append(calibration_folder + 'nf_dark_' + str(i) + '_10ms.tif')
    nf_ambient_images.append(image_folder + 'nf_ambient_' + str(i) + '_1.8ms.tif')

nf_flat_images = []
#for i in xrange(8):
#    nf_flat_images.append(calibration_folder + 'nf_flat_' + str(i) + '_1ms.tif')

LED_nf = IA.ImageAnalysis(nf_led_images, nf_dark_images, nf_flat_images, nf_ambient_images, 3.45, 10, None, 16)
print 'LED image initialized'
laser_nf_agitated = IA.ImageAnalysis(nf_laser_images_agitated, nf_dark_images, nf_flat_images, nf_ambient_images, 3.45, 10, None, 16)
print 'Agitated laser image initialized'
laser_nf = IA.ImageAnalysis(nf_laser_images, nf_dark_images, nf_flat_images, nf_ambient_images, 3.45, 10, None, 16)
print 'Unagitated laser image initialized'
print

laser_fp = FP.FiberProperties(None, laser_nf, None)
laser_fp_agitated = FP.FiberProperties(None, laser_nf_agitated, None)
LED_fp = FP.FiberProperties(None, LED_nf, None)

deg = 16
radius_factor = 0.95
#os.system('mkdir "Modal Noise Data"')
for method in ['fft', 'gaussian', 'polynomial', 'gradient', 'gini', 'tophat', 'contrast']:
    print 'Method:', method

    led_mn = LED_fp.getModalNoise('near', method, radius_factor=radius_factor, deg=deg)
    print 'LED:', led_mn

    laser_ag_mn = laser_fp_agitated.getModalNoise('near', method, radius_factor=radius_factor, deg=deg)
    print 'Agitated Laser:', laser_ag_mn

    laser_unag_mn = laser_fp.getModalNoise('near', method, radius_factor=radius_factor, deg=deg)
    print 'Unagitated Laser:', laser_unag_mn

    #new_file = open('./Modal Noise Data/gaussianimage_' + method + '.txt','w')
    #data = {'LED': led_mn, 'Agitated Laser': laser_ag_mn, 'Unagitated Laser': laser_unag_mn}
    #new_file.write(str(data))
    #new_file.close()
    print
