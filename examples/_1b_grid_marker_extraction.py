from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt

def generate_2D_grid_image_for_testing(size=100, grid_spacing=20):
    """
    generates an artifical 2D image (or array) of size [size, size].
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

# data_loc = Path(r'\where\is\your\data')
# mr_volume = MarkerVolume(data_loc)


'''
BUT maybe what's going to be easier for you, is to start with something simpler. So, I've
put a function in here which will generate a very simple 2D image:
'''
artificial_2D_grid_image = generate_2D_grid_image_for_testing()
plt.figure()
plt.grid(False)
plt.imshow(artificial_2D_grid_image)
plt.show()
plt.grid(False)
'''
previous authors have mentioned 'prewit operators' to process grid search.
I would google those if I were you.
skimage seems to have a nice implementation
'''