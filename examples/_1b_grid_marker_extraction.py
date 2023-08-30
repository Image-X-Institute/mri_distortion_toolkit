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

def plot_2D_array(image_to_plot, title=None, extent=None):
    # BW2: I added extent as a default argument here, which will enable you to pass it in.
    # notes on calculating it below
    plt.figure()
    plt.grid(False)
    plt.imshow(image_to_plot, extent=extent)
    if title:
        plt.title(title)
    plt.show()

def plot_2D_scatter(x_centroids, y_centroids):
    # BW2 this is a scatter plot
    # you will have to figure out how to overlay this onto an image!
    plt.figure()
    plt.scatter(x_centroids, y_centroids)
    plt.axis('equal')
    plt.show()

def detect_if_edge_case(x_center, y_center,X, Y, InputSlice):
    """
    BW2: fin, note what we are doing here when we split this into functions.
    because all this function does is return True or False, it gices us a lot of
    freedom down the track to swap it out for something different.
    This is formally referred to as 'modular' programming and is one of the easiest thigns
    you can do to make your life easier down the track!

    there are actually two distinct tasks here:
    1. get a cluster of pixels
    2. run a test on that cluster.
    If you want some coding practice, you should try splitting this function up into two different functions
    """
    # get the cluster of pixels around the centroid
    r_roi = 7 # we should be more intelligent about how to choose this but for now we can hard code it
    # you might have to fine tune a bit too.
    # we want to detect pixels where:
    # cenroid_x - r_roi <= centroid_x <= centroid_x <= centroid_x + r_roi
    # and similar for y
    x_ind = np.logical_and((X - r_roi <= x_center), (X+r_roi >= x_center))
    #BW2 super useful when doing logical operation like this:
    print(f'Total elements in X: {x_ind.size} number of True values in x_ind: {np.count_nonzero(x_ind)}')
    y_ind = np.logical_and((Y - r_roi <= y_center), ( Y+ r_roi >= y_center))
    combined_ind = np.logical_and(x_ind, y_ind)  # cases where both x and y meet our critera
    # can sanity check what pixels you've selected:
    # plot_2D_array(combined_ind)
    pixels = InputSlice[combined_ind]


    # fin to do: apply a test to pixels to see if we should count it or not. a few ideas:
    # average pixel value
    # number of peak in a histogram
    # derivitive accross the image
    # some combo of the above?

    return True  # BW2 for now will always be True, you nee to update based on the test you write


# call on a single dicom image from the grid-based MRI folder
# dicom_path = r"C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\MR\1.3.46.670589.11.79127.5.0.6984.2022112517535313493.dcm"
# BW: I created a folder and copied just one image into it:
# 1.3.46.670589.11.79127.5.0.6984.2022112517535358579.dcm
dicom_path = Path(r"C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\Test Slice")
#BW2: comment out the below to revert to your own data locatoin
dicom_path = Path(r'/home/brendan/Downloads/3Done phantom/'
                     r'3DOne_zzphys_MR_2022-11-25_173121_post.upgrade_T1.3D.Tra.HN.L_n208__00000/')
InputVolume, dicom_affine, (X, Y, Z) = dicom_to_numpy(dicom_path,
                                                      FilesToReadIn='1.3.46.670589.11.79127.5.0.6984.2022112517535358579.dcm',
                                                      file_extension='dcm',return_XYZ=True)
#BW2: to calculate extent of image:
extent = (X.min(), X.max(), Y.max(), Y.min())
# you can pass this to imshow (or plot_2D_array) to get the correct coordinates


# BW we can remove this singleton dimension like this:
InputSlice = InputVolume.squeeze()  # dont delete this line!
X = X.squeeze()
Y = Y.squeeze()

# use prewit operators to find intersections:
v_edge_map = skimage.filters.prewitt_v(InputSlice)

# highlight the gradient change down the vertical lines i.e., where the horizontal lines intersect
intersect_map = skimage.filters.prewitt_h(v_edge_map)
# take absolute values
intersect_map = abs(intersect_map)

# blurring the image points using a Guassian filter
blurred_map = gaussian_filter(intersect_map, sigma=1)


threshold = skimage.filters.threshold_otsu(blurred_map)

binary_image = blurred_map > threshold

# BW: create a plot of our images so far:
#BW2 I will update the extent in one
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10,10])
axs[0, 0].imshow(InputSlice, extent=extent); axs[0, 0].grid(False); axs[0, 0].set_title('original')
axs[0, 1].imshow(intersect_map); axs[0, 1].grid(False); axs[0, 1].set_title('prewit operators')
axs[1, 0].imshow(blurred_map); axs[1, 0].grid(False); axs[1, 0].set_title('blurred')
axs[1, 1].imshow(binary_image); axs[1, 1].grid(False); axs[1, 1].set_title('otsu threshold')
plt.show()
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
    # weighted_x = np.sum(np.multiply(X[RegionInd], InputSlice[RegionInd]))
    # weighted_y = np.sum(np.multiply(Y[RegionInd], InputSlice[RegionInd]))
    # x_centroids.append(weighted_x / region_sum)
    # y_centroids.append(weighted_y / region_sum)

    # BW2: the above lines of code are intended to weight the centroid of the detected point towards the brightest part of the
    # image. it wasn't working for you for some reason, I've just deleted and replaced with the simpler logic below:
    # x_centroid = mean of X coordinate of region. you can delete all commented lines here once you read this.

    potential_x_centroid = np.mean(X[RegionInd])
    potential_y_centroid = np.mean(Y[RegionInd])
    # BW2: now, we want to detect and remove edge cases:
    # for debugging could do this:
    show_debug_image = True  # note will make a lot of iamges!!
    if show_debug_image:
        plt.figure()
        plt.imshow(InputSlice, extent=extent)
        plt.scatter(potential_x_centroid, potential_y_centroid)
        plt.show()
    proper_centroid = detect_if_edge_case(potential_x_centroid, potential_y_centroid, X, Y, InputSlice)

    if proper_centroid:
        # only count the ones which aren't detected on the edge
        x_centroids.append(potential_x_centroid)
        y_centroids.append(potential_y_centroid)

plt.figure()
extent = (X.min(), X.max(), Y.max(), Y.min())
plt.imshow(InputSlice, extent=extent)
plt.grid(False)
plt.scatter(x_centroids, y_centroids, marker='x')
plt.show()






centroids = np.array([x_centroids, y_centroids]).T
print(skipped)
print(centroids)

# scatter plot of the centroids
plot_2D_scatter(x_centroids, y_centroids)
# overlay by updating