from FiberProperties import FRD, FRD_Input, imageList, ImageAnalysis
import matplotlib.pyplot as plt
import numpy as np
from sys import platform
import cPickle as pickle
from multiprocessing import Pool

NEW_DATA = True
FRD_CALIBRATION_THRESHOLD = 1500

def imageListFRD(image_name, f_ratios, **kwargs):
    return [imageList(image_name+str(f)+'/im_', **kwargs) for f in f_ratios]

def darkFiles(folder):
    return imageList(folder+'Dark/im_')

def ambientFiles(folder):
    return imageList(folder+'Ambient/im_')

def inputFiles(folder, f):
    return imageList(folder+'Input '+str(f)+'/im_')

def outputFiles(folder, f):
    return imageList(folder+'Output '+str(f)+'/im_')

def inputObjects(folder, in_f):
    if NEW_DATA:
        return [ImageAnalysis(inputFiles(folder, f),
                              dark=darkFiles(folder),
                              ambient=ambientFiles(folder),
                              input_fnum=f,
                              threshold=FRD_CALIBRATION_THRESHOLD
                              ) for f in in_f]
    return [folder+'Input '+str(f)+'/im_data.pkl' for f in in_f]

def outputObjects(folder, out_f):
    if NEW_DATA:
        return [ImageAnalysis(inputFiles(folder, f),
                              dark=darkFiles(folder),
                              ambient=ambientFiles(folder),
                              output_fnum=f,
                              threshold=FRD_CALIBRATION_THRESHOLD
                              ) for f in out_f]
    return [folder+'Output '+str(f)+'/im_data.pkl' for f in out_f]

class Container(object):
    def __init__(name, folder, in_f, out_f):
        self.name = name
        self.folder = folder
        self.in_obj = inputObjects(folder, in_f)
        self.out_obj = outputObject(folder, out_f)
        self.output = None

if __name__ == '__main__':
    if platform == 'darwin':
        folder = '/home/ryanp/Fiber_Characterization/'
    else:
        folder = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/frd/Core Extension/'

    FOCAL_RATIO_DIAMETER = 0.95

    in_f = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    out_f = [3.0, 4.0, 5.0]

    reference = Container('Reference', folder + '2016-08-10 Reference Octagonal/', in_f, out_f)
    prototype_1 = Container('Prototype 1', folder + '2016-08-05 Prototype Core Extension 1/' , in_f[1:], [3.5])
    prototype_2 = Container('Prototype 2', folder + '2016-08-09 Prototype Core Extension 2/', in_f, out_f)
    prototype_A2 = Container('Prototype A2', folder + '2017-01-11 Prototype A2/', in_f, out_f)
    prototype_A3 = Container('Prototype A3', folder + '2017-01-12 Prototype A3/', in_f, out_f)

    tests = [reference, prototype_1, prototype_2, prototype_A2, prototype_A3]

    for test in tests:
        test.output = FRD(test.in_objs, test.out_objs, 'edge', True,
                          fnum_diameter=FOCAL_RATIO_DIAMETER)

    for test in tests:
        frd_info = test.output[0]
        magn = test.output[1]
        magn_list = test.output[2]
        magn_error = test.output[3]

        plt.figure(1)
        for i, f in enumerate(frd_info.input_focal_ratios):
            plt.plot(frd_info.encircled_energy_focal_ratios[i],
                       frd_info.encircled_energy[i],
                       label=str(f))
        plt.xlabel('Output f/#')
        plt.ylabel('Encircled Energy')
        plt.ylim(ymax=1)
        plt.xticks()
        plt.yticks()
        plt.grid()
        plt.legend(loc=3, title='Input f/#')
        plt.title('FRD: ' + test.name)
        plt.savefig(test.folder + test.name + ' FRD.png')

        plt.figure(2)
        plt.errorbar(frd_info.input_fnum,
                     frd_info.energy_loss,
                     xerr=magn_error*np.array(frd_info.input_fnum),
                     label=test.name)

        plt.figure(3)
        plt.errorbar(frd_info.input_fnum,
                     frd_info.output_fnum,
                     xerr=magn_error*np.array(frd_info.input_fnum),
                     yerr=magn_error*np.array(frd_info.input_fnum),
                     label=test.name)

        plt.figure(4)
        for i, f in enumerate([2.5, 3.0, 3.5, 4.0, 4.5, 5.0]):
            if f in frd_info.input_fnum:
                plt.subplot(3, 2, i+1)
                index = frd_info.input_fnum.index(f)
                plt.plot(frd_info.encircled_energy_fnum[index],
                         frd_info.encircled_energy[index],
                         label=test.name,
                         linewidth=1)
                plt.xlabel('Output f/#', fontsize=8)
                plt.ylabel('Encircled Energy', fontsize=8)
                plt.ylim(ymax=1)
                plt.xticks(fontsize=8)
                plt.yticks(fontsize=8)
                plt.grid()
                plt.legend(loc=3, fontsize=8)
                plt.title('Input f/# = ' + str(f), fontsize=8)

    plt.figure(2)
    plt.xlabel('Input f/#')
    plt.ylabel('Energy Loss [%]')
    plt.grid()
    plt.legend(loc=4)
    plt.title('Energy Loss at Constant F/#')
    plt.savefig(folder + 'Energy Loss.png')

    plt.figure(3)
    plt.plot(tests[-1].input_focal_ratios, tests[-1].input_focal_ratios,
             label='Ideal', linestyle='--', color='black')
    plt.xlabel('Input f/#')
    plt.ylabel('Output f/#')
    plt.grid()
    plt.legend(loc=2)
    plt.title('FRD Comparison')
    plt.savefig(folder + 'Input vs Output.png')

    plt.figure(4)
    plt.savefig(folder + 'Encircled Energy Comparison.png')
