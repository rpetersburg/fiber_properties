from fiber_properties import FiberImage, plot_modal_noise, show_plots, save_plot
from modal_noise_script import save_modal_noise_data, object_file
import numpy as np

METHOD = 'filter'
CAMS = ['nf', 'ff']
CASE = 1
FOLDER = "C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/"

if CASE == 1:
    TITLE = 'Coupled Fiber'
    FOLDER += 'coupled_fibers/'
    TITLES = ['200um-200um',
              '100um-200um',
              'Oct-Circ',
              '200um']
    FOLDERS = [FOLDER + '200-200um_test2/',
               FOLDER + '100-200um/',
               FOLDER + 'oct-circ-200um/',
               FOLDER + '../Kris_data/circular_200um/']
    TESTS = [['agitated_first', 'agitated_second', 'agitated_both'],
             ['agitated_first_100um', 'agitated_second_200um', 'agitated_both'],
             ['agitated_oct', 'agitated_circ', 'agitated_both'],
             ['linear_agitation' for i in xrange(3)]]
    LABELS = ['first agitated', 'second agitated', 'both agitated']

if CASE == 2:
    TITLE = 'Fiber Geometry'
    FOLDER += 'Kris_data/'
    TITLES = ['Circular 100um', 'Octagonal 100um', 'Circular 200um', 'Octagonal 200um', 'Rectangular 100x300um']
    FOLDERS = [FOLDER + 'circular_100um/',
               FOLDER + 'octagonal_100um/',
               FOLDER + 'circular_200um/',
               FOLDER + 'octagonal_200um/',
               FOLDER + 'rectangular_100x300um/']
    TESTS = [['linear_agitation', 'circular_agitation', 'coupled_agitation'] for i in xrange(len(FOLDERS))]
    LABELS = ['linear agitation', 'circular agitation', 'coupled agitation']

if CASE == 3:
    TITLE = 'Amplitude vs Frequency'
    FOLDER += 'amp_freq_600um/'
    TITLES = ['']
    FOLDERS = [FOLDER]
    TESTS = [['agitated_5volts_40mm_10s', 'agitated_5volts_160mm_10s_test1', 'agitated_30volts_40mm_10s', 'agitated_30volts_160mm_10s_test1']]
    LABELS = ['0.1Hz 40mm agitation', '0.1Hz 160mm agitation', '1.0Hz 40mm agitation', '1.0Hz 160mm agitation']

if CASE == 4:
    TITLE = 'Amplitude vs Frequency'
    FOLDER += 'amp_freq_200um/'
    TITLES = ['']
    FOLDERS = [FOLDER]
    TESTS = [['agitated_5volts_40mm_8s', 'agitated_30volts_40mm_1s', 'agitated_5volts_160mm_8s', 'agitated_30volts_160mm_1s', 'baseline']]
    LABELS = ['0.1Hz 40mm', '1.0Hz 40mm', '0.1Hz 160mm', '1.0Hz 160mm', 'LED source']

if CASE == 5:
    TITLE = 'All'
    FOLDER += 'amp_freq_200um/'
    TITLES = ['']
    FOLDERS = [FOLDER]
    TESTS = [['unagitated_1s',
              'unagitated_8s',
              'agitated_5volts_40mm_1s',
              'agitated_5volts_40mm_8s',
              'agitated_5volts_40mm_80s',
              'agitated_30volts_40mm_1s',
              'agitated_30volts_40mm_10s',
              'agitated_5volts_160mm_1s',
              'agitated_5volts_160mm_8s',
              'agitated_5volts_160mm_80s',
              'agitated_30volts_160mm_1s',
              'agitated_30volts_160mm_10s',
              'baseline']]
    LABELS = ['unagitated 1s',
              'unagitated 8s',
              '0.1Hz 40mm 1s',
              '0.1Hz 40mm 8s',
              '0.1Hz 40mm 80s',
              '1.0Hz 40mm 1s',
              '1.0Hz 40mm 10s',
              '0.1Hz 160mm 1s',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 80s',
              '1.0Hz 160mm 1s',
              '1.0Hz 160mm 10s',
              'baseline']

if CASE == 6:
    TITLE = 'All'
    FOLDER += 'amp_freq_600um/'
    TITLES = ['']
    FOLDERS = [FOLDER]
    TESTS = [['unagitated_1s',
              'unagitated_8s',
              'unagitated_10s',
              'agitated_5volts_40mm_10s',
              'agitated_5volts_160mm_8s',
              'agitated_5volts_160mm_10s_test1',
              'agitated_5volts_160mm_10s_test2',
              'agitated_5volts_160mm_80s',
              'agitated_30volts_40mm_10s',
              'agitated_30volts_160mm_1s',
              'agitated_30volts_160mm_10s_test1',
              'agitated_30volts_160mm_10s_test2',
              'baseline_amp_freq']]
    LABELS = ['unagitated 1s',
              'unagitated 8s',
              'unagitated 10s',
              '0.1Hz 40mm 10s',
              '0.1Hz 160mm 8s',
              '0.1Hz 160mm 10s test1',
              '0.1Hz 160mm 10s test2',
              '0.1Hz 160mm 80s',
              '1.0Hz 40mm 10s',
              '1.0Hz 160mm 1s',
              '1.0Hz 160mm 10s test1',
              '1.0Hz 160mm 10s test2',
              'baseline']

# def object_file(folder, test, cam):
#     return folder + test + '/' + cam + '_obj.pkl'

if __name__ == '__main__':
    for cam in CAMS:       
        # save_modal_noise_data(FOLDERS, TESTS, cam, labels, METHODS, TITLES) 
        modal_noise = []
        for title, folder, test in zip(TITLES, FOLDERS, TESTS):
            modal_noise.append([])
            for t in test:
                # mn = []
                # for im in xrange(10):
                #     im_obj = FiberImage(object_file(folder + t, cam, 1, im))
                #     mn.append(im_obj.get_modal_noise(method=METHOD))
                # mn = np.array(mn).mean()
                im_obj = FiberImage(object_file(folder + t, cam, 10, 0))
                mn = im_obj.get_modal_noise(method=METHOD)
                modal_noise[-1].append(mn)
                print t, mn

        plot_modal_noise(modal_noise, plot_type='bar', bar_labels=LABELS, labels=TITLES, method=METHOD)
        save_plot(FOLDER + 'analysis/' + TITLE + ' ' + cam.upper() + ' ' + METHOD + ' SNR.png')
        # save_plot(FOLDER + 'analysis/' + TITLE + ' ' + cam.upper() + ' ' + METHOD + ' SNR.pdf')

