# from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume
# from pathlib import Path
import numpy as np
import skimage.filters
from matplotlib import pyplot as plt
import scipy.ndimage
import pydicom

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
# generate a 2D numpy array of pixel values representing a slice of a grid-based phantom dicom folder
# pixel_array = generate_2D_grid_image_for_testing()

# OR call on a single dicom image from the grid-based MRI folder
dicom_path = r"C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\MR\1.3.46.670589.11.79127.5.0.6984.2022112517535313493.dcm"
# dicom_path = r"C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\MR\1.3.46.670589.11.79127.5.0.6984.2022112517535367593.dcm"
dicom_data = pydicom.dcmread(dicom_path)

# don't know how to get coordinates

# convert the Dicom image to a 2D numpy array containing the pixel data
if 'PixelData' in dicom_data:
    pixel_array = dicom_data.pixel_array
else:
    raise ValueError("DICOM file does not contain pixel data.")

# plt.figure()
# plt.grid(False)
# plt.imshow(pixel_array)
# plt.show()

# highlight the gradient change of the vertical lines on the slice
v_edge_map = skimage.filters.prewitt_v(pixel_array)

# plt.figure()
# plt.grid(False)
# plt.imshow(v_edge_map)
# plt.show()

# highlight the gradient change down the vertical lines i.e., where the horizontal lines intersect
intersect_map = skimage.filters.prewitt_h(v_edge_map)

# plt.figure()
# plt.grid(False)
# plt.imshow(intersect_map)
# plt.show()

# make all values positive to highlight the intersection points
intersect_map = abs(intersect_map)

# blurring the image points using a Guassian filter
blurred_map = scipy.ndimage.gaussian_filter(intersect_map, sigma=1)

# plt.figure()
# plt.grid(False)
# plt.imshow(blurred_map)
# plt.show()

# clear image by setting values < cut-off (0.03 for generated array and 0.0014 for Dicom image) to zero and all others to 1
cut_off = 0.0014
for i in range(len(blurred_map)): # iterating rows
    for j in range(len(blurred_map[0])): # iterating cols
        if blurred_map[i][j] > cut_off: # if the value is greater than the cut-off value, set it to 1
            blurred_map[i][j] = int(1)
            continue
        blurred_map[i][j] = int(0) # otherwise, set the value to 0

# plt.figure()
# plt.grid(False)
# plt.imshow(blurred_map)
# plt.show()
# plt.grid(False)

# blob search or label function

from math import sqrt
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage.color import rgb2gray


image_gray = rgb2gray(blurred_map)

blobs_log = blob_log(image_gray, max_sigma=30, num_sigma=10, threshold=.1)

# Compute radii in the 3rd column.
blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)

blobs_dog = blob_dog(image_gray, max_sigma=30, threshold=.1)
blobs_dog[:, 2] = blobs_dog[:, 2] * sqrt(2)

blobs_doh = blob_doh(image_gray, max_sigma=30, threshold=.01)

blobs_list = [blobs_log, blobs_dog, blobs_doh]
colors = ['yellow', 'lime', 'red']
titles = ['Laplacian of Gaussian', 'Difference of Gaussian',
          'Determinant of Hessian']
sequence = zip(blobs_list, colors, titles)

fig, axes = plt.subplots(1, 3, figsize=(9, 3), sharex=True, sharey=True)
ax = axes.ravel()

for idx, (blobs, color, title) in enumerate(sequence):
    ax[idx].set_title(title)
    ax[idx].imshow(image)
    for blob in blobs:
        y, x, r = blob
        c = plt.Circle((x, y), r, color=color, linewidth=2, fill=False)
        ax[idx].add_patch(c)
    ax[idx].set_axis_off()

plt.tight_layout()
plt.show()


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
# PATH OF GRIP DICOM "C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\MR\1.3.46.670589.11.79127.5.0.6984.2022112517535313493.dcm"

# want to return [x_coordinates, y_coordinates] of each intersection point

# distance between x coordinates: x_spacing = width of dicom image / len(artificial_2D_grid_image[0])
# distance between y coordinates: y_spacing = length of dicom image / len(artificial_2D_grid_image)

#for i in len(artificial_2D_grid_image[0]): # cycle through positions in row
    # save c_coord if value = 1
#    for j in len(artificial_2D_grid_image):  # cycle through positions in column
        # save y coord if value = 1
        # check if x coord is also = 1
        # if yes, save coordinate as grid intersection point