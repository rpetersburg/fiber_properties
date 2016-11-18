from FiberProperties import ImageAnalysis, Calibration, modalNoise, plotFFT, showPlot, savePlot, showImageArray
import numpy as np

if __name__ == '__main__':
    base_folder = '../data/modal_noise/2016-11-11/'
    ext = '.fit'
    cameras = ['nf', 'ff']
    tests = ['unagitated', 'agitated_5volts_40mm', 'agitated_30volts_40mm', 'agitated_5volts_160mm', 'agitated_30volts_160mm', 'baseline']
    fft = dict.fromkeys(cameras, {})
    freq = dict.fromkeys(cameras, {})

    for camera in cameras:
        camera_folder = base_folder + camera + '/'
        calibration = Calibration(dark=[camera_folder + 'dark_' + str(i).zfill(3) + ext for i in xrange(10)],
                                  ambient=[camera_folder + 'ambient_' + str(i).zfill(3) + ext for i in xrange(10)])
        print camera + ' calibration initialized'

        for test in tests:
            if test == 'baseline':
                continue
            images = [camera_folder + test + '_' + str(i).zfill(3) + ext for i in xrange(10)]
            im_obj = ImageAnalysis(images, calibration, camera=camera)
            im_obj.saveImages(file_name=test, folder=camera_folder)
            fft[camera][test], freq[camera][test] = modalNoise(im_obj,
                                                               method='fft',
                                                               output='array',
                                                               radius_factor=1.0,
                                                               show_image=False)

        if 'baseline' in tests:
            if camera == 'nf':
                perfect_image = im_obj.getTophatFit()
            elif camera == 'ff':
                perfect_image = im_obj.getGaussianFit()
                perfect_image *= (perfect_image > 0.0).astype('float64')

            baseline_image = np.zeros_like(perfect_image)
            for i in xrange(10):
                baseline_image += np.random.normal(perfect_image, np.sqrt(perfect_image+0.001)) / 10

            baseline_obj = ImageAnalysis(baseline_image,
                                         pixel_size=im_obj.getPixelSize(),
                                         camera=camera)
            baseline_obj.saveImages(file_name='baseline', folder=camera_folder)
            fft[camera]['baseline'], freq[camera]['baseline'] = modalNoise(baseline_obj,
                                                                           method='fft',
                                                                           output='array',
                                                                           radius_factor=1.0,
                                                                           show_image=False)


        plotFFT([freq[camera][test] for test in tests],
                [fft[camera][test] for test in tests],
                labels=tests,
                title=camera.upper() + ' Modal Noise Comparison (600um Fiber)')
        savePlot(base_folder + camera + '_modal_noise')
        showPlot()

