from fiber_properties import (frd, image_list, FiberImage, save_plot,
                              plot_frd_encircled_energy,
                              plot_frd_encircled_energy_comparison,
                              plot_frd_input_output,
                              plot_frd_energy_loss)

NEW_OBJECTS = False
NEW_DATA = False
FOCAL_RATIO_DIAMETER = 0.95
FRD_CALIBRATION_THRESHOLD = 1500

class Container(object):
    def __init__(self, name, folder, in_f, out_f):
        self.name = name
        self.folder = folder
        self.in_objs = input_objects(folder, in_f)
        self.out_objs = output_objects(folder, out_f)
        self.output = None

def image_list_frd(image_name, f_ratios, **kwargs):
    return [image_list(image_name+str(f)+'/im_', **kwargs) for f in f_ratios]

def dark_files(folder):
    return image_list(folder+'dark/ff_')

def ambient_files(folder):
    return image_list(folder+'ambient/ff_')

def input_files(folder, f):
    return image_list(folder+'input_'+str(f)+'/ff_')

def output_files(folder, f):
    return image_list(folder+'output_'+str(f)+'/ff_')

def input_objects(folder, in_f):
    if NEW_OBJECTS:
        for f in in_f:
            print('Saving ' + folder + 'input_' + str(f))
            im_obj = FiberImage(input_files(folder, f),
                                ambient=ambient_files(folder),
                                input_fnum=f,
                                threshold=FRD_CALIBRATION_THRESHOLD,
                                camera='ff')
            im_obj.save_object(folder+'input_'+str(f)+'/ff_object.pkl')
            im_obj.save_image(folder+'input_'+str(f)+'/ff_corrected.fit')
    return [folder+'input_'+str(f)+'/ff_object.pkl' for f in in_f]

def output_objects(folder, out_f):
    if NEW_OBJECTS:
        for f in out_f:
            print('Saving ' + folder + 'output_' + str(f))
            im_obj = FiberImage(output_files(folder, f),
                                ambient=ambient_files(folder),
                                output_fnum=f,
                                threshold=FRD_CALIBRATION_THRESHOLD,
                                camera='ff')
            im_obj.save_object(folder+'output_'+str(f)+'/ff_object.pkl')
            im_obj.save_image(folder+'output_'+str(f)+'/ff_corrected.fit')
    return [folder+'output_'+str(f)+'/ff_object.pkl' for f in out_f]

if __name__ == '__main__':
    CASE = 4
    in_f = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    out_f = [3.0, 4.0, 5.0]

    if CASE == 1:
        TITLE = 'Rectangular'
        FOLDER = '../data/EXPRES/rectangular_132/frd2/'
        TESTS = [Container(TITLE, FOLDER, in_f, out_f)]
    if CASE == 2:
        TITLE = 'Octagonal'
        FOLDER = '../data/EXPRES/bare_octagonal/frd/'
        TESTS = [Container(TITLE, FOLDER, in_f, out_f)]
    if CASE == 3:
        TITLE = 'Rect vs Oct'
        FOLDER = '../data/EXPRES/'
        octagonal = Container('Octagonal', FOLDER+'bare_octagonal/frd/', in_f, out_f)
        rectangular = Container('Rectangular', FOLDER+'rectangular_132/frd2/', in_f, out_f)
        TESTS = [octagonal, rectangular]
    if CASE == 4:
        TITLE = 'Rec_test'
        FOLDER = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/frd/bare_rec_fiber/'
        TESTS = [Container(TITLE, FOLDER, in_f, out_f)]

    # diameters = [0.95, 0.97, 0.99]
    # frd_outputs = []
    # labels = []
    # for test in TESTS:
    #     for diameter in diameters:
    #         print 'FRD for diameter', diameter
    #         frd_outputs.append(frd(TESTS[0].in_objs, TESTS[0].out_objs,
    #                                cal_method='edge', save_objs=True,
    #                                fnum_diameter=diameter, new=True))
    #         labels.append(str(diameter*100) + '%')
    #     plot_frd_input_output(frd_outputs, labels)
    #     save_plot(test.folder + 'Input vs. Output Comparison.png')

    for test in TESTS:
        print('Calculating FRD for '+ test.name + ' Fiber...')
        test.output = frd(test.in_objs, test.out_objs,
                          cal_method='edge', save_objs=True,
                          fnum_diameter=FOCAL_RATIO_DIAMETER, new=NEW_DATA)
        plot_frd_encircled_energy(test.output)
        save_plot(test.folder + test.name + ' FRD.png')

    frd_outputs = [test.output for test in TESTS]
    labels = [test.name for test in TESTS]

    plot_frd_energy_loss(frd_outputs, labels)
    save_plot(FOLDER + 'Energy Loss.png')

    plot_frd_input_output(frd_outputs, labels)
    save_plot(FOLDER + 'Input vs Output.png')

    plot_frd_encircled_energy_comparison(frd_outputs, labels)
    save_plot(FOLDER + 'Encircled Energy Comparison.png')
