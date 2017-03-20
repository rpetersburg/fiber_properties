from fiber_properties import frd, image_list, ImageAnalysis
import matplotlib.pyplot as plt
import numpy as np
from sys import platform
import cPickle as pickle
from multiprocessing import Pool

NEW_DATA = False
FRD_CALIBRATION_THRESHOLD = 1500

def image_list_frd(image_name, f_ratios, **kwargs):
    return [image_list(image_name+str(f)+'/im_', **kwargs) for f in f_ratios]

def dark_files(folder):
    return image_list(folder+'../dark/ff_')

def ambient_files(folder):
    return image_list(folder+'../ambient/ff_')

def input_files(folder, f):
    return image_list(folder+'input_'+str(f)+'/ff_')

def output_files(folder, f):
    return image_list(folder+'output_'+str(f)+'/ff_')

def input_objects(folder, in_f):
    if NEW_DATA:
        output = []
        for f in in_f:
            im_obj = ImageAnalysis(input_files(folder, f),
                                   dark=dark_files(folder),
                                   ambient=ambient_files(folder),
                                   input_fnum=f,
                                   threshold=FRD_CALIBRATION_THRESHOLD,
                                   camera='ff')
            im_obj.save()
            output.append(im_obj)
        return output
    return [folder+'input_'+str(f)+'/ff_object.pkl' for f in in_f]

def output_objects(folder, out_f):
    if NEW_DATA:
        output = []
        for f in out_f:
            im_obj = ImageAnalysis(output_files(folder, f),
                                   dark=dark_files(folder),
                                   ambient=ambient_files(folder),
                                   output_fnum=f,
                                   threshold=FRD_CALIBRATION_THRESHOLD,
                                   camera='ff')
            im_obj.save()
            output.append(im_obj)
        return output
    return [folder+'output_'+str(f)+'/ff_object.pkl' for f in out_f]

class Container(object):
    def __init__(self, name, folder, in_f, out_f):
        self.name = name
        self.folder = folder
        self.in_objs = input_objects(folder, in_f)
        self.out_objs = output_objects(folder, out_f)
        self.output = None

if __name__ == '__main__':
    folder = '../data/EXPRES/rectangular_132/frd2/'

    FOCAL_RATIO_DIAMETER = 0.95

    in_f = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    out_f = [3.0, 4.0, 5.0]

    octagonal = Container('Rectangular', '', in_f, out_f)

    # plt.figure()
    # for diameter in [0.95, 0.98, 0.99]:
    #     print 'FRD for diameter', diameter
    #     frd_info, magn, magn_list, magn_error = frd(octagonal.in_objs,
    #                                                 octagonal.out_objs,
    #                                                 cal_method='edge',
    #                                                 save_objs=True,
    #                                                 fnum_diameter=diameter,
    #                                                 new=True)

    #     plt.errorbar(frd_info.input_fnum,
    #                  frd_info.output_fnum,
    #                  xerr=magn_error*np.array(frd_info.input_fnum),
    #                  yerr=magn_error*np.array(frd_info.input_fnum),
    #                  label=str(diameter*100) + '%')
    # plt.xlabel('Input f/#')
    # plt.ylabel('Output f/#')
    # plt.xticks()
    # plt.yticks()
    # plt.grid()
    # plt.legend(loc='best')
    # plt.title('Octagonal FRD')
    # plt.savefig('Octagonal Input vs. Output.png')

    tests = [octagonal]

    for test in tests:
        print 'Calculating FRD for '+test.name+' Fiber'
        test.output = frd(test.in_objs, test.out_objs, 'edge', True,
                          fnum_diameter=FOCAL_RATIO_DIAMETER, new=True)
    for test in tests:
        frd_info = test.output[0]
        magn = test.output[1]
        magn_list = test.output[2]
        magn_error = test.output[3]

        plt.figure(1)
        for i, f in enumerate(frd_info.input_fnum):
            plt.plot(frd_info.encircled_energy_fnum[i],
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
    plt.plot(tests[-1].output[0].input_fnum, tests[-1].output[0].input_fnum,
             label='Ideal', linestyle='--', color='black')
    plt.xlabel('Input f/#')
    plt.ylabel('Output f/#')
    plt.grid()
    plt.legend(loc=2)
    plt.title('FRD Comparison')
    plt.savefig(folder + 'Input vs Output.png')

    plt.figure(4)
    plt.savefig(folder + 'Encircled Energy Comparison.png')
