from fiber_properties import FiberImage

# image = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/' \
#         + 'modal_noise/coupled_fibers/200-200um_test2/agitated_both/nf_obj.pkl'
image = 'C:/Libraries/Box Sync/ExoLab/Fiber_Characterization/Image Analysis/data/modal_noise/Kris_data/rectangular_100x300um/coupled_agitation/nf_obj.pkl'


im_obj = FiberImage(image)
# im_obj.get_modal_noise(method='fft', new=True, show_image=True)
im_obj.set_modal_noise(method='filter', new=True, show_image=False, kernel_size=31)
# print im_obj.get_modal_noise(method='filter')
im_obj.save_object(image)