# from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume
# from pathlib import Path
import numpy as np
import skimage.filters
from matplotlib import pyplot as plt
from scipy import ndimage, datasets

def generate_2D_grid_image_for_testing(size=100, grid_spacing=20):
    """
    generates an artificial 2D image (or array) of size [size, size] being [rows, cols].
    the image is nearly all zeros, but has some rows of ones defined by grid spacing

    :param size:
    :param grid_spacing:
    :return:
    """

    image = np.zeros([size, size])  # square image of all zeros

    for i in range(size):  # loop over rows
        if i == 0:
            continue  # skip first row
        if i % grid_spacing == 0:  #nb % in python is modulus
            image[i, :] = 1  # whenever we hit grid_spacing interval overwrite zeros with ones
    for j in range(size):
        if j == 0:
            continue
        if j % grid_spacing == 0:
            image[:, j] = 1

    return image


'''
I've made the below far simpler: removed any inputs you're not going to need, and removed any image processing
not specific to grid. Basically now, this stops at the point we had our debug statement. This would be a good point 
for you to start trying to write some code.
For now, you can just work directly in this example. we can integrate into the main code later.
mr_volume.InputVolume
'''

# data_loc = Path(r'C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data')
# mr_volume = MarkerVolume(data_loc)

'''
BUT maybe what's going to be easier for you, is to start with something simpler. So, I've
put a function in here which will generate a very simple 2D image:
'''
artificial_2D_grid_image = generate_2D_grid_image_for_testing()
# print(artificial_2D_grid_image)
plt.figure()
plt.grid(False)
v_edge_map = skimage.filters.prewitt_v(artificial_2D_grid_image)
intersect_map = skimage.filters.prewitt_h(v_edge_map)
plt.imshow(intersect_map)
# plt.imshow(v_edge_map)
# plt.imshow(artificial_2D_grid_image)
plt.show()
plt.grid(False)
'''
previous authors have mentioned 'prewit operators' to process grid search.
I would google those if I were you.
skimage seems to have a nice implementation
'''

# generated gradient map for vert and horiz
# abs value
# blur the image, smooth points into blob. e.g., gaussian
# set points with value << X to zero
# blob search or label function -- look at skimage

# once settled with generated image - use a single dicom image (slice) - test in dicom_to_numpy

'''
'''


# want to return [x_coordinates, y_coordinates] of each intersection point

# distance between x coordinates: x_spacing = width of dicom image / len(artificial_2D_grid_image[0])
# distance between y coordinates: y_spacing = length of dicom image / len(artificial_2D_grid_image)

#for i in len(artificial_2D_grid_image[0]): # cycle through positions in row
    # save c_coord if value = 1
#    for j in len(artificial_2D_grid_image):  # cycle through positions in column
        # save y coord if value = 1
        # check if x coord is also = 1
        # if yes, save coordinate as grid intersection point