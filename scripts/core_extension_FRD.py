from FiberProperties import FRD, FRD_Input
import matplotlib.pyplot as plt
import numpy as np
from sys import platform
import cPickle as pickle
from multiprocessing import Pool

NEW_DATA = False
PARALLELIZATION = True

if __name__ == '__main__':
    if platform == 'darwin':
        base_folder = '/home/ryanp/Fiber_Characterization/'
    else:
        base_folder = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/frd/Core Extension/'

    reference = FRD_Input('Reference Fiber',
                          base_folder+'2016-08-10 Reference Octagonal/')
    prototype_1 = FRD_Input('Prototype Fiber 1',
                            base_folder+'2016-08-05 Prototype Core Extension 1/',
                            input_focal_ratios=[3.0, 3.5, 4.0, 4.5, 5.0],
                            cal_focal_ratios=[3.5])
    prototype_2 = FRD_Input('Prototype Fiber 2',
                            base_folder+'2016-08-09 Prototype Core Extension 2/')
    prototype_A2 = FRD_Input('Prototype Fiber A2',
                             base_folder+'2017-01-11 Prototype A2/')
    prototype_A3 = FRD_Input('Prototype Fiber A3',
                             base_folder+'2017-01-12 Prototype A3/')

    tests = [reference, prototype_1, prototype_2, prototype_A2, prototype_A3]

    if NEW_DATA:
        if PARALLELIZATION:
            pool = Pool(processes=5)
            output = pool.map(FRD, tests)
            for i in xrange(len(output)):
                tests[i].output = output[i]
        else:
            for test in tests:
                test.output = FRD(test, save_images=True)
    else:
        for test in tests:
            with open(test.folder + 'FRD_Output.pkl') as file:
                test.output = pickle.load(file)

    plt.figure()
    for test in tests:
        for i, f in enumerate(test.input_focal_ratios):
            plt.subplot(3, 2, i)
            plt.plot(test.output.encircled_energy_focal_ratios[i],
                     test.output.encircled_energy[i],
                     label=str(f))
        plt.xlabel('Output f/#')
        plt.ylabel('Encircled Energy')
        plt.legend(title='Input f/#', loc=3)
        plt.title('FRD: ' + test.name)
        plt.savefig(test.folder + test.name + ' FRD.png')

    plt.figure()
    for test in tests:
        plt.errorbar(test.input_focal_ratios,
                     test.output.energy_loss,
                     xerr=test.output.magn_error*np.array(test.input_focal_ratios),
                     label=test.name)
    plt.xlabel('Input f/#')
    plt.ylabel('Energy Loss [%]')
    plt.grid()
    plt.legend(loc=4)
    plt.title('Energy Loss at Constant F/#')
    plt.savefig(base_folder + 'Energy Loss.png')

    plt.figure()
    for test in tests:
        output = test.output
        plt.errorbar(test.input_focal_ratios,
                     test.output.output_focal_ratios,
                     xerr=test.output.magn_error*np.array(test.input_focal_ratios),
                     yerr=test.output.magn_error*np.array(test.output.output_focal_ratios),
                     label=test.name)
    plt.plot(tests[-1].input_focal_ratios, tests[-1].input_focal_ratios,
             label='Ideal', linestyle='--', color='black')
    plt.xlabel('Input f/#')
    plt.ylabel('Output f/#')
    plt.grid()
    plt.legend(loc=2)
    plt.title('FRD Comparison')
    plt.savefig(base_folder + 'Input vs Output.png')

    plt.figure()
    for test in tests:
        for i, f in enumerate([2.5, 3.0, 3.5, 4.0, 4.5, 5.0]):
            if f in test.input_focal_ratios:
                plt.subplot(3, 2, i+1)
                index = test.input_focal_ratios.index(f)
                plt.plot(test.output.encircled_energy_focal_ratios[index],
                         test.output.encircled_energy[index],
                         label=test.name,
                         linewidth=1)
                plt.xlabel('Output f/#', fontsize=8)
                plt.ylabel('Encircled Energy', fontsize=8)
                plt.xticks(fontsize=8)
                plt.yticks(fontsize=8)
                plt.grid()
                plt.legend(loc=3, fontsize=8)
                plt.title('Input f/# = ' + str(f), fontsize=8)
    plt.savefig(base_folder + 'Encircled Energy Comparison.png')




