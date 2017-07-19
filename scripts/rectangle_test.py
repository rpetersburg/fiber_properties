from fiber_properties import FiberImage, image_list

obj = FiberImage(image_list('pos_1/nf_'),
                    ambient=image_list('../ambient/nf_'),
                    camera='ff')

print obj.get_fiber_center(method='rectangle', show_image=True)