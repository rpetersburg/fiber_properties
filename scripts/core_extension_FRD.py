from fiber_properties import (frd, image_list, FiberImage, save_plot,
                              plot_frd_encircled_energy,
                              plot_frd_encircled_energy_comparison,
                              plot_frd_input_output,
                              plot_frd_energy_loss)
import numpy as np
from sys import platform
import cPickle as pickle
from multiprocessing import Pool

NEW_OBJECTS = False
NEW_DATA = False
FRD_CALIBRATION_THRESHOLD = 1500

def image_list_frd(image_name, f_ratios, **kwargs):
    return [image_list(image_name+str(f)+'/im_', **kwargs) for f in f_ratios]

def dark_files(folder):
    return image_list(folder+'Dark/im_')

def ambient_files(folder):
    return image_list(folder+'Ambient/im_')

def input_files(folder, f):
    return image_list(folder+'Input '+str(f)+'/im_')

def output_files(folder, f):
    return image_list(folder+'Output '+str(f)+'/im_')

def input_objects(folder, in_f):
    if NEW_OBJECTS:
        for f in in_f:
            print 'Saving ' + folder + 'Input ' + str(f)
            im_obj = FiberImage(input_files(folder, f),
                                dark=dark_files(folder),
                                ambient=ambient_files(folder),
                                input_fnum=f,
                                threshold=FRD_CALIBRATION_THRESHOLD,
                                camera='ff')
            im_obj.save_object(folder+'Input '+str(f)+'/ff_object.pkl')
            im_obj.save_image(folder+'Input '+str(f)+'/ff_corrected.fit')
    return [folder+'Input '+str(f)+'/ff_object.pkl' for f in in_f]

def output_objects(folder, out_f):
    if NEW_OBJECTS:
        for f in out_f:
            print 'Saving ' + folder + 'Output ' + str(f)
            im_obj = FiberImage(output_files(folder, f),
                                dark=dark_files(folder),
                                ambient=ambient_files(folder),
                                output_fnum=f,
                                threshold=FRD_CALIBRATION_THRESHOLD,
                                camera='ff')
            im_obj.save_object(folder+'Output '+str(f)+'/ff_object.pkl')
            im_obj.save_image(folder+'Output '+str(f)+'/ff_corrected.fit')
    return [folder+'Output '+str(f)+'/ff_object.pkl' for f in out_f]

class Container(object):
    def __init__(self, name, folder, in_f, out_f):
        self.name = name
        self.folder = folder
        self.in_objs = input_objects(folder, in_f)
        self.out_objs = output_objects(folder, out_f)
        self.output = None

if __name__ == '__main__':
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
        print 'Calculating FRD for ' + test.name + ' Fiber'
        test.output = frd(test.in_objs, test.out_objs,
                          cal_method='edge', save_objs=True,
                          fnum_diameter=FOCAL_RATIO_DIAMETER, new=NEW_DATA)
        plot_frd_encircled_energy(test.output)
        save_plot(test.folder + test.name + ' FRD.png')

    frd_outputs = [test.output for test in tests]
    labels = [test.name for test in tests]

    plot_frd_energy_loss(frd_outputs, labels)
    save_plot(folder + 'Energy Loss.png')

    plot_frd_input_output(frd_outputs, labels)
    save_plot(folder + 'Input vs Output.png')

    plot_frd_encircled_energy_comparison(frd_outputs, labels)
    save_plot(folder + 'Encircled Energy Comparison.png')
