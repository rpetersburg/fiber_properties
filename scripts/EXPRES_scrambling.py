from fiber_properties import (scrambling_gain, image_list,
                              plot_scrambling_gain_input_output,
                              plot_scrambling_gain, save_plot, show_plots,
                              load_image_object, FiberImage)

NEW_DATA = True
FOLDER = '/Users/Dominic/Box Sync/Fiber_Characterization/Image Analysis/data/EXPRES/science_rec/scrambling/'
POSITIONS = ['pos_cent', 'pos_left', 'pos_right']

if __name__ == '__main__':

    if NEW_DATA:
        from fiber_properties import FiberImage

        #in_dark = image_list(FOLDER + '../dark/in_')
        in_ambient = image_list(FOLDER + 'ambient/in_')
        #nf_dark = image_list(FOLDER + '../dark/nf_')
        nf_ambient = image_list(FOLDER + 'ambient/nf_')
        #ff_dark = image_list(FOLDER + '../dark/ff_')
        ff_ambient = image_list(FOLDER + 'ambient/ff_')

        for pos in POSITIONS:
            print('Initializing ' + pos)
            in_images = image_list(FOLDER + pos + '/in_')
            nf_images = image_list(FOLDER + pos + '/nf_')
            ff_images = image_list(FOLDER + pos + '/ff_')

            FiberImage(in_images, ambient=in_ambient, camera='in', threshold=180).save_object(FOLDER + pos + '/in_object.pkl')
            FiberImage(nf_images, ambient=nf_ambient, camera='nf').save_object(FOLDER + pos + '/nf_object.pkl')
            FiberImage(ff_images, ambient=ff_ambient, camera='ff').save_object(FOLDER + pos + '/ff_object.pkl')

    in_objs = [FOLDER + pos + '/in_object.pkl' for pos in POSITIONS]

    nf_objs = [FOLDER + pos + '/nf_object.pkl' for pos in POSITIONS]
    nf_scrambling = scrambling_gain(in_objs, nf_objs, input_method='full', output_method='full')
    print(nf_scrambling.in_d)
    print(nf_scrambling.out_d)
    print(nf_scrambling.scrambling_gain)
    
    print('Minimum NF scrambling:', min(nf_scrambling.scrambling_gain))
    print('Maximum NF scrambling:', max(nf_scrambling.scrambling_gain))

    plot_scrambling_gain_input_output(nf_scrambling)
    save_plot(FOLDER + 'Near Field Shift.png')
    plot_scrambling_gain(nf_scrambling)
    save_plot(FOLDER + 'Near Field SG.png')
    show_plots()

    ff_objs = [FOLDER + pos + '/ff_object.pkl' for pos in POSITIONS]
    ff_scrambling = scrambling_gain(in_objs, ff_objs, input_method='full', output_method='full')

    print('Minimum FF scrambling:', min(ff_scrambling.scrambling_gain))
    print('Maximum FF scrambling:', max(ff_scrambling.scrambling_gain))

    plot_scrambling_gain_input_output(ff_scrambling)
    save_plot(FOLDER + 'Far Field Shift.png')
    plot_scrambling_gain(ff_scrambling)
    save_plot(FOLDER + 'Far Field SG.png')
    show_plots()