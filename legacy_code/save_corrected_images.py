from FiberProperties import ImageAnalysis, Calibration

if __name__ == "__main__":
    folder = '../data/Modal Noise Measurements/2016-10-19/'

    data = {'nf': {}, 'ff': {}}

    for cam in ['nf']:
        calibration = Calibration(ambient = [folder + cam + '_ambient_' + str(i).zfill(3) + '.fit' for i in xrange(10)])
        for method in ['agitated', 'unagitated']:
            images = [folder + 'oct_60_' + cam + '_' + method + '_' + str(i).zfill(3) + '.fit' for i in xrange(10)]
            im_obj = ImageAnalysis(images, calibration, camera=cam)
            im_obj.saveImages(file_name='oct_60_'+cam+'_'+method, folder=folder)

