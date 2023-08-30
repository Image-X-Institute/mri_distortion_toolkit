from pathlib import Path
from mri_distortion_toolkit.utilities import dicom_to_numpy
import numpy as np
import skimage.filters
from matplotlib import pyplot as plt
# import scipy.ndimage  # BW: prefer the below style of input
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
# dicom_path = r"C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\MR\1.3.46.670589.11.79127.5.0.6984.2022112517535313493.dcm"
# BW: I created a folder and copied just one image into it:
# 1.3.46.670589.11.79127.5.0.6984.2022112517535358579.dcm
dicom_path = Path(r"C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\Test Slice")

InputVolume, dicom_affine, (X, Y, Z) = dicom_to_numpy(dicom_path, file_extension='dcm',return_XYZ=True)

# BW fin; you can delete the following line, for demo only:
print(f'the shape of the input volume is {InputVolume.shape}')
# BW note that the last dimension is 1; this is because the code above reads a volume, but we only gave it one image.
# BW we can remove this singleton dimension like this:
InputSlice = InputVolume.squeeze()  # dont delete this line!
print(f'the shape of the input slice is {InputSlice.shape}')
# BW in case you want to plot I added a function:
# plot_2D_array(InputSlice)

# use prewit operators to find intersections:
v_edge_map = skimage.filters.prewitt_v(InputSlice)

# highlight the gradient change down the vertical lines i.e., where the horizontal lines intersect
intersect_map = skimage.filters.prewitt_h(v_edge_map)
# take absolute values
intersect_map = abs(intersect_map)

# blurring the image points using a Guassian filter
blurred_map = gaussian_filter(intersect_map, sigma=1)

# clear image by setting values < cut-off (0.03 for generated array and 0.0014 for Dicom image) to zero and all others to 1
# bw: I suggest for determining the cut off, you use otsu's method:
# http://devdoc.net/python/scikit-image-doc-0.13.1/auto_examples/xx_applications/plot_thresholding.html
# calculate optimal threshold using Otsu's method
threshold = skimage.filters.threshold_otsu(blurred_map)

# BW: then to threshold the image you just need to do this:
binary_image = blurred_map > threshold

# BW: create a plot of our images so far:
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10,10])
axs[0, 0].imshow(InputSlice); axs[0, 0].grid(False); axs[0, 0].set_title('original')
axs[0, 1].imshow(intersect_map); axs[0, 1].grid(False); axs[0, 1].set_title('prewit operators')
axs[1, 0].imshow(blurred_map); axs[1, 0].grid(False); axs[1, 0].set_title('blurred')
axs[1, 1].imshow(binary_image); axs[1, 1].grid(False); axs[1, 1].set_title('otsu threshold')
plt.show()
# blob search or label function
# BW: label seems to work fine.
label_image = label(binary_image)

# plot_2D_array(label_image)

voxel_min, voxel_max = get_marker_max_min_volume(label_image)

# find contour centroids
x_centroids = []
y_centroids = []
skipped = 0

unique_labels = np.unique(label_image)[1:]  # first label is background so skip

# extract x,y of each connected region
for i, label_level in enumerate(unique_labels):
    RegionInd = np.equal(label_image, label_level)
    voxels = np.count_nonzero(RegionInd)
    if voxels < voxel_min or voxels > voxel_max:
        skipped += 1
        continue # skip outliers

    region_sum = np.sum(InputSlice[RegionInd])
    weighted_x = np.sum(np.multiply(X[RegionInd], InputSlice[RegionInd]))
    weighted_y = np.sum(np.multiply(Y[RegionInd], InputSlice[RegionInd]))
    x_centroids.append(weighted_x / region_sum)
    y_centroids.append(weighted_y / region_sum)

centroids = np.array([x_centroids, y_centroids]).T
print(skipped)
print(centroids)

# scatter plot of the centroids
# overlay by updating