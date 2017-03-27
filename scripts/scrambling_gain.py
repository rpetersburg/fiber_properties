from fiber_properties import (scrambling_gain, image_list,
                              plot_scrambling_gain_input_output,
                              plot_scrambling_gain, save_plot, show_plots,
                              load_image_object, FiberImage)

if __name__ == '__main__':

    NEW_DATA = False
    folder = '../data/scrambling/2016-08-05 Prototype Core Extension 1/'

    if NEW_DATA:
        from fiber_properties import FiberImage

        in_dark = image_list(folder + 'Dark/in_')
        in_ambient = image_list(folder + 'Ambient/in_')
        nf_dark = image_list(folder + 'Dark/nf_')
        nf_ambient = image_list(folder + 'Ambient/nf_')
        ff_dark = image_list(folder + 'Dark/ff_')
        ff_ambient = image_list(folder + 'Ambient/ff_')

        for shift in ['00', '05', '10', '15', '20', '25', '30']:
            print 'Initializing Shift ' + shift
            in_images = image_list(folder + 'Shift_' + shift + '/in_')
            nf_images = image_list(folder + 'Shift_' + shift + '/nf_')
            ff_images = image_list(folder + 'Shift_' + shift + '/ff_')

            FiberImage(in_images, in_dark, in_ambient, camera='in').save()
            FiberImage(nf_images, nf_dark, nf_ambient, camera='nf').save()
            FiberImage(ff_images, ff_dark, ff_ambient, camera='ff').save()

    shifts = ['00', '05', '10', '15', '20', '25', '30']
    in_objs = [folder + 'Shift_' + shift + '/in_object.pkl' for shift in shifts]
    a = load_image_object(in_objs[0])

    nf_objs = [folder + 'Shift_' + shift + '/nf_object.pkl' for shift in shifts]
    nf_scrambling = scrambling_gain(in_objs, nf_objs, input_method='radius', output_method='radius')

    print 'Minimum NF scrambling: ' + min(nf_scrambling.scrambling_gain)
    print 'Maximum NF scrambling: ' + max(nf_scrambling.scrambling_gain)

    plot_scrambling_gain_input_output(nf_scrambling, 'NF Centroid Shift')
    save_plot(folder + 'Near Field Shift.png')
    plot_scrambling_gain(nf_scrambling, 'NF Scrambling Gains')
    save_plot(folder + 'Near Field SG.png')
    show_plots()

    ff_objs = [folder + 'Shift_' + shift + '/ff_object.pkl' for shift in shifts]
    ff_scrambling = scrambling_gain(in_objs, ff_objs, input_method='radius', output_method='gaussian')

    print 'Minimum FF scrambling: ' + min(ff_scrambling.scrambling_gain)
    print 'Maximum FF scrambling: ' + max(ff_scrambling.scrambling_gain)

    plot_scrambling_gain_input_output(ff_scrambling, 'FF Centroid Shift')
    save_plot(folder + 'Far Field Shift.png')
    plot_scrambling_gain(ff_scrambling, 'FF Scrambling Gains')
    save_plot(folder + 'Far Field SG.png')
    show_plots()