from FiberProperties import ImageAnalysis, Calibration, modalNoise
import os as os
import matplotlib.pyplot as plt
from copy import deepcopy
plt.rc('font', size=14, family='sans-serif')
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)
plt.rc('lines', lw=2)


if __name__ == '__main__':
    base_folder = '../data/modal_noise/2016-11-10/'
    ext = '.fit'

    data = {}

    for camera in ['nf','ff']:
        calibration = Calibration(dark=[base_folder + camera + '/dark_' + str(i).zfill(3) + ext for i in xrange(10)],
                                  ambient=[base_folder + camera + '/ambient_' + str(i).zfill(3) + ext for i in xrange(10)])
        print camera + ' calibration initialized'

        empty_data = {'images': [], 'fft': [], 'freq': []}
        for test in ['agitated', 'unagitated']:
            images = [base_folder + camera + '/' + test + '_' + str(i).zfill(3) + ext for i in xrange(20)]

            im_obj = ImageAnalysis(images, calibration)

            fft, freq = modalNoise(im_obj, method='fft', output='array', radius_factor=1.0)

        nf_test_obj = nf['agitated']['obj']
        nf['baseline']['obj'] = ImageAnalysis(nf_test_obj.getTophatFit(),
                                              pixel_size=nf_test_obj.getPixelSize(),
                                              threshold = 0.1,
                                              camera='nf')

        for test in ['agitated', 'unagitated']:
            ff[test]['obj'] = ImageAnalysis(ff[test]['images'], ff['calibration'])
        ff_test_obj = ff['agitated']['obj']
        ff['baseline']['obj'] = ImageAnalysis(ff_test_obj.getGaussianFit(),
                                              pixel_size=ff_test_obj.getPixelSize(),
                                              magnification=1,
                                              camera='ff')

        print

        for cam in [nf, ff]:
            for test in ['agitated', 'unagitated', 'baseline']:
                print test + ':'
                cam[test]['fft'], cam[test]['freq'] = modalNoise(cam[test]['obj'], method='fft', output='array', radius_factor=1.0)

            plotFFT([cam[test]['freq'] for test in ['agitated', 'unagitated', 'baseline']],
                    [cam[test]['fft'] for test in ['agitated', 'unagitated', 'baseline']],
                    labels=['Agitated laser', 'Unagitated laser', 'Baseline'],
                    title=cam[test]['obj'].getCamera().upper() + ' Modal Noise Comparison (600um Fiber)')

