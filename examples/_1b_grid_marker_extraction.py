from pathlib import Path
from mri_distortion_toolkit.utilities import dicom_to_numpy
import numpy as np
import skimage.filters
from matplotlib import pyplot as plt
import scipy.ndimage  # BW: prefer the below style of input
from scipy.ndimage import gaussian_filter
import pydicom
import sys
from skimage.measure import label

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

def get_marker_max_min_volume(labels, marker_size_lower_tol=1, marker_size_upper_tol=20):
    """
    figures out the range of voxels that a marker should have
    """
    n_voxels = []  # going to keep track of this so we can remove any very small regions if needed
    unique_labels = np.unique(labels)[1:]  # first label is background so skip
    for label_level in unique_labels:  # first label is background so skip
        # extract x,y,z of each connected region
        RegionInd = labels == label_level
        n_voxels.append(np.count_nonzero(RegionInd))
    # Set up min and max marker volumes
    n_voxels_median = np.median(np.array(n_voxels))
    voxel_min = (1 - marker_size_lower_tol) * n_voxels_median
    voxel_max = (1 + marker_size_upper_tol) * n_voxels_median

    return voxel_min, voxel_max

def plot_2D_array(image_to_plot, title=None):
    plt.figure()
    plt.grid(False)
    plt.imshow(image_to_plot)
    if title:
        plt.title(title)
    plt.show()



'''
BUT maybe what's going to be easier for you, is to start with something simpler. So, I've
put a function in here which will generate a very simple 2D image:
'''
# generate a 2D numpy array of pixel values representing a slice of a grid-based phantom dicom folder
# pixel_array = generate_2D_grid_image_for_testing()

# OR call on a single dicom image from the grid-based MRI folder
dicom_path = r"C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\MR\1.3.46.670589.11.79127.5.0.6984.2022112517535313493.dcm"
# BW: I created a folder and copied just one image into it:
# 1.3.46.670589.11.79127.5.0.6984.2022112517535358579.dcm
dicom_path = Path(r'/home/brendan/Downloads/simple_data')

InputVolume, dicom_affine, (X, Y, Z) = dicom_to_numpy(dicom_path,
                                                     file_extension='dcm',
                                                     return_XYZ=True)

# BW fin; you can delete the following line, for demo only:
print(f'the shape of the input volume is {InputVolume.shape}')
# BW note that the last dimension is 1; this is because the code above reads a volume, but we only gave it one image.
# BW we can remove this singleton dimension like this:
InputSlice = InputVolume.squeeze()  # dont delete this line!
print(f'the shape of the input slice is {InputSlice.shape}')
# BW in case you want to plot I added a function:
plot_2D_array(InputSlice)


# use prewit operators to find intersections:
v_edge_map = skimage.filters.prewitt_v(InputSlice)

# highlight the gradient change down the vertical lines i.e., where the horizontal lines intersect
intersect_map = skimage.filters.prewitt_h(v_edge_map)
# take absolute values
intersect_map = abs(intersect_map)

# blurring the image points using a Guassian filter
blurred_map = gaussian_filter(intersect_map, sigma=1)



plt.show()

# clear image by setting values < cut-off (0.03 for generated array and 0.0014 for Dicom image) to zero and all others to 1
cut_off = 0.0014
# bw: I suggest for determinging the cut off, you use otsu's method:
# http://devdoc.net/python/scikit-image-doc-0.13.1/auto_examples/xx_applications/plot_thresholding.html
# BW: then to threshold the image you just need to do this:
binary_image = blurred_map > cut_off

# BW: create a plot of our images so far:
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10,10])
axs[0, 0].imshow(InputSlice); axs[0, 0].grid(False); axs[0, 0].set_title('original')
axs[0, 1].imshow(intersect_map); axs[0, 1].grid(False); axs[0, 1].set_title('prewit operators')
axs[1, 0].imshow(blurred_map); axs[1, 0].grid(False); axs[1, 0].set_title('blurred')
axs[1, 1].imshow(binary_image); axs[1, 1].grid(False); axs[1, 1].set_title('otsu threshold')



# blob search or label function
# BW: label seems to work fine.

from math import sqrt
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage.color import rgb2gray


# image_gray = rgb2gray(blurred_map)
# BW the reason you get an error is your image is already greyscale.
# rgb implies three seperate channels, i.e. the image would have size [x_pixel, y_pixels, 3]
# if you are concerned that imshow shows iamges in color instead of grey, you need to look up
# colormaps instead!

blobs_log = blob_log(binary_image, max_sigma=30, num_sigma=10, threshold=.1)

# Compute radii in the 3rd column.
blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)

blobs_dog = blob_dog(binary_image, max_sigma=30, threshold=.1)
blobs_dog[:, 2] = blobs_dog[:, 2] * sqrt(2)

blobs_doh = blob_doh(binary_image, max_sigma=30, threshold=.01)

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