from FiberProperties import frdFromTestData, FRD_Test
import matplotlib.pyplot as plt
import numpy as np
from ast import literal_eval
from sys import platform

if __name__ == '__main__':
    if platform == 'darwin':
        base_folder = '/home/ryanp/Fiber_Characterization/'
    else:
        base_folder = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/frd/Core Extension/'

    new = False

    ref = FRD_Test('Reference Fiber',
                   base_folder+'2016-08-10 Reference Octagonal/',
                   new=new)
    prototype_1 = FRD_Test('Prototype Fiber 1',
                       base_folder+'2016-08-05 Prototype Core Extension 1/',
                       new=new,
                       input_focal_ratios=[3.0, 3.5, 4.0, 4.5, 5.0],
                       cal_focal_ratios=[3.5])
    prototype_2 = FRD_Test('Prototype Fiber 2',
                        base_folder+'2016-08-09 Prototype Core Extension 2/',
                        new=new)
    prototype_A2 = FRD_Test('Prototype Fiber A2',
                        base_folder+'2017-01-11 Prototype A2/',
                        new=new)
    prototype_A3 = FRD_Test('Prototype Fiber A3',
                        base_folder+'2017-01-12 Prototype A3/',
                        new=new)

    tests = [ref, prototype_1, prototype_2, prototype_A2, prototype_A3]

    for test in tests:
        print test.name + ':'
        if test.new:
            test.info = frdFromTestData(test)
        else:
            with open(test.folder + 'info.txt') as file:
                test.info = literal_eval(file.read())

    plt.figure()
    for test in tests:
        output = test.info
        plt.errorbar(output['input_focal_ratios'],
                     output['energy_loss'],
                     xerr=output['error']*np.array(output['input_focal_ratios']),
                     label=test.name)
    plt.xlabel('Input f/#')
    plt.ylabel('Energy Loss [%]')
    plt.grid()
    plt.legend(loc=4)
    plt.title('Energy Loss at Constant F/#')
    plt.savefig(base_folder + 'Energy Loss.png')

    plt.figure()
    for test in tests:
        output = test.info
        plt.errorbar(output['input_focal_ratios'],
                     output['output_focal_ratios'],
                     xerr=output['error']*np.array(output['input_focal_ratios']),
                     yerr=output['error']*np.array(output['output_focal_ratios']),
                     label=test.name)
    plt.plot(tests[-1].info['input_focal_ratios'], tests[-1].info['input_focal_ratios'], label='Ideal', linestyle='--', color='black')
    plt.xlabel('Input f/#')
    plt.ylabel('Output f/#')
    plt.grid()
    plt.legend(loc=2)
    plt.title('FRD Comparison')
    plt.savefig(base_folder + 'Input vs Output.png')

    plt.figure()
    for tests in tests:
        output = test.info
        for i, focal_ratio in enumerate([2.5, 3.0, 3.5, 4.0, 4.5, 5.0]):
            if focal_ratio in output['input_focal_ratios']:
                plt.subplot(3, 2, i+1)
                index = output['input_focal_ratios'].index(focal_ratio)
                plt.plot(output['encircled_energy_focal_ratios'][index],
                         output['encircled_energy'][index],
                         label=test.name,
                         linewidth=1)
                plt.xlabel('Output f/#', fontsize=8)
                plt.ylabel('Encircled Energy', fontsize=8)
                plt.xticks(fontsize=8)
                plt.yticks(fontsize=8)
                plt.grid()
                plt.legend(loc=3, fontsize=8)
                plt.title('Input f/# = ' + str(focal_ratio), fontsize=8)
    plt.savefig(base_folder + 'Encircled Energy Comparison.png')

