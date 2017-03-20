from fiber_properties import (scrambling_gain, image_list,
                              plot_scrambling_gain_input_output,
                              plot_scrambling_gain, save_plot, show_plots,
                              load_image_object, ImageAnalysis)

NEW_DATA = True
FOLDER = '../data/EXPRES/rectangular_132/scrambling/'
POSITIONS = ['pos_1', 'pos_2', 'pos_3', 'pos_4', 'pos_5']

if __name__ == '__main__':

    if NEW_DATA:
        from fiber_properties import ImageAnalysis

        in_dark = image_list(FOLDER + '../dark/in_')
        in_ambient = image_list(FOLDER + '../ambient/in_')
        nf_dark = image_list(FOLDER + '../dark/nf_')
        nf_ambient = image_list(FOLDER + '../ambient/nf_')
        ff_dark = image_list(FOLDER + '../dark/ff_')
        ff_ambient = image_list(FOLDER + '../ambient/ff_')

        for pos in POSITIONS:
            print 'Initializing ' + pos
            in_images = image_list(FOLDER + pos + '/in_')
            nf_images = image_list(FOLDER + pos + '/nf_')
            ff_images = image_list(FOLDER + pos + '/ff_')

            in_obj = ImageAnalysis(in_images, in_dark, in_ambient, camera='in')
            nf_obj = ImageAnalysis(nf_images, nf_dark, nf_ambient, camera='nf')
            ff_obj = ImageAnalysis(ff_images, ff_dark, ff_ambient, camera='ff')

            in_obj.

            ImageAnalysis(in_images, in_dark, in_ambient, camera='in').save()
            ImageAnalysis(nf_images, nf_dark, nf_ambient, camera='nf').save()
            ImageAnalysis(ff_images, ff_dark, ff_ambient, camera='ff').save()

    in_objs = [FOLDER + pos + '/in_object.pkl' for pos in POSITIONS]

    nf_objs = [FOLDER + pos + '/nf_object.pkl' for pos in POSITIONS]
    nf_scrambling = scrambling_gain(in_objs, nf_objs, input_method='radius', output_method='radius')

    print 'Minimum NF scrambling: ' + min(nf_scrambling.scrambling_gain)
    print 'Maximum NF scrambling: ' + max(nf_scrambling.scrambling_gain)

    plot_scrambling_gain_input_output(nf_scrambling, 'NF Centroid Shift')
    save_plot(FOLDER + 'Near Field Shift.png')
    plot_scrambling_gain(nf_scrambling, 'NF Scrambling Gains')
    save_plot(FOLDER + 'Near Field SG.png')
    show_plots()

    ff_objs = [FOLDER + pos + '/ff_object.pkl' for pos in POSITIONS]
    ff_scrambling = scrambling_gain(in_objs, ff_objs, input_method='radius', output_method='gaussian')

    print 'Minimum FF scrambling: ' + min(ff_scrambling.scrambling_gain)
    print 'Maximum FF scrambling: ' + max(ff_scrambling.scrambling_gain)

    plot_scrambling_gain_input_output(ff_scrambling, 'FF Centroid Shift')
    save_plot(FOLDER + 'Far Field Shift.png')
    plot_scrambling_gain(ff_scrambling, 'FF Scrambling Gains')
    save_plot(FOLDER + 'Far Field SG.png')
    show_plots()