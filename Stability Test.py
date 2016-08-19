import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
from ImageAnalysis import ImageAnalysis
from datetime import datetime

def plotDiameterStability(info_dict, title):
    plt.figure()
    plt.title(title)
    plt.plot(info_dict['agitated']['time'], info_dict['agitated']['diameter'], label='agitated')
    plt.plot(info_dict['unagitated']['time'], info_dict['unagitated']['diameter'], label='unagitated')
    plt.xlabel('Time')
    plt.ylabel('Diameter [um]')
    plt.legend()
    plt.show()

def plotCenterStability(info_dict, title):
    plt.figure()
    plt.subplot(211)    
    plt.title(title)
    plt.plot(info_dict['agitated']['time'], info_dict['agitated']['x0'], label='agitated')
    plt.plot(info_dict['unagitated']['time'], info_dict['unagitated']['x0'], label='unagitated')
    plt.xlabel('Time')
    plt.ylabel('Position [um]')
    plt.legend(title='Center X')
    plt.subplot(212)
    plt.plot(info_dict['agitated']['time'], info_dict['agitated']['y0'], label='agitated')
    plt.plot(info_dict['unagitated']['time'], info_dict['unagitated']['y0'], label='unagitated')
    plt.xlabel('Time')
    plt.ylabel('Position [um]')
    plt.legend(title='Center Y')
    plt.show()

if __name__ == "__main__":
    folder = {}
    folder['unagitated'] = 'Stability Measurements/2016-08-15 Stability Test Unagitated/Data2/'
    folder['agitated'] = 'Stability Measurements/2016-08-16 Stability Test Agitated/Data/'

    data_dict = {'x0': [], 'y0': [], 'diameter': [], 'time': []}
    nf_dict = {'agitated': deepcopy(data_dict), 'unagitated': deepcopy(data_dict)}
    ff_dict = deepcopy(nf_dict)
    in_dict = deepcopy(nf_dict)

    for test in ['unagitated', 'agitated']:
        for i in xrange(100):
            nf_data = folder[test] + 'nf_' + str(i).zfill(3) + '_data.p'
            nf_obj = ImageAnalysis(image_input=None, image_data=nf_data)
            y0, x0, diameter = nf_obj.getFiberData(method='radius', units='microns')
            nf_dict[test]['x0'].append(x0)
            nf_dict[test]['y0'].append(y0)
            nf_dict[test]['diameter'].append(diameter)
            nf_dict[test]['time'].append(nf_obj.getImageInfo('date_time'))

            ff_data = folder[test] + 'ff_' + str(i).zfill(3) + '_data.p'
            ff_obj = ImageAnalysis(image_input=None, image_data=ff_data)
            y0, x0, diameter = ff_obj.getFiberData(method='gaussian', units='microns')
            ff_dict[test]['x0'].append(x0)
            ff_dict[test]['y0'].append(y0)
            ff_dict[test]['diameter'].append(diameter)
            ff_dict[test]['time'].append(ff_obj.getImageInfo('date_time'))

            in_data = folder[test] + 'in_' + str(i).zfill(3) + '_data.p'
            in_obj = ImageAnalysis(image_input=None, image_data=in_data)
            y0, x0, diameter = in_obj.getFiberData(method='gaussian', units='microns')
            in_dict[test]['x0'].append(x0)
            in_dict[test]['y0'].append(y0)
            in_dict[test]['diameter'].append(diameter)
            in_dict[test]['time'].append(in_obj.getImageInfo('date_time'))

        for cam_dict in [nf_dict, ff_dict, in_dict]:
            for prop in ['x0', 'y0', 'diameter', 'time']:
                cam_dict[test][prop] = np.array(cam_dict[test][prop]) - cam_dict[test][prop][0]
            for i, time in enumerate(cam_dict[test]['time']):
                cam_dict[test]['time'][i] = time.total_seconds()

    print nf_dict['unagitated']['time']

    plotDiameterStability(nf_dict, 'Near Field Diameter Stability')
    plotDiameterStability(ff_dict, 'Far Field Diameter Stability')
    plotDiameterStability(in_dict, 'Fiber Input Diameter Stability')

    plotCenterStability(nf_dict, 'Near Field Center Stability')
    plotCenterStability(ff_dict, 'Far Field Center Stability')
    plotCenterStability(in_dict, 'Fiber Input Center Stability')
