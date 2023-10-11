from pathlib import Path
from mri_distortion_toolkit.utilities import dicom_to_numpy
import numpy as np
import skimage.filters
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks
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
    plt.colorbar()
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
    because all this function does is return True or False, it gives us a lot of
    freedom down the track to swap it out for something different.
    This is formally referred to as 'modular' programming and is one of the easiest things
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
    # centroid_x - r_roi <= centroid_x <= centroid_x <= centroid_x + r_roi
    # and similar for y
    x_ind = np.logical_and((X - r_roi <= x_center), (X+r_roi >= x_center))
    #BW2 super useful when doing logical operation like this:
    # print(f'Total elements in X: {x_ind.size} number of True values in x_ind: {np.count_nonzero(x_ind)}')
    y_ind = np.logical_and((Y - r_roi <= y_center), ( Y+ r_roi >= y_center))
    combined_ind = np.logical_and(x_ind, y_ind)  # cases where both x and y meet our criteria
    # print(f'Total elements: {combined_ind.size} number of True values: {np.count_nonzero(combined_ind)}')
    # can sanity check what pixels you've selected:
    # plot_2D_array(combined_ind, title = 'Potential Centroid Indexes of Interest')

    pixels = InputSlice[combined_ind] # pixels is a 1D array of pixel values within the ROI

    region_of_int = InputSlice*combined_ind # region_of_int is a 2D array of only the pixels within the ROI

    # otsu method to find the optimal threshold value for region of interest
    threshold = skimage.filters.threshold_otsu(region_of_int)

    # apply the threshold
    thresholded_region = region_of_int > threshold

    # fin to do: apply a test to pixels to see if we should count it or not. a few ideas:
    # average pixel value
    proper_centroid1 = mean_detection(pixels, x_center, y_center)

    # percentage of non zero pixels in ROI after thresholding
    #proper_centroid2 = percent_detection(thresholded_region, combined_ind, x_center, y_center)

    # number of peak in a histogram
    #proper_centroid3 = peak_detection(thresholded_region, x_center, y_center)
    # the intersection points should have a peak in each corner whereas edge cases should be missing at least one peak

    # derivitive across the image
    # some combo of the above

    return proper_centroid1 # and proper_centroid2 and proper_centroid3

    # if abs(x_center) < 128 and abs(x_center) > 127 and abs(y_center) < 32 and abs(y_center) > 31:
    #     fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10, 10])
    #     axs[0, 0].imshow(region_of_int); axs[0, 0].grid(False); axs[0, 0].set_title('x -122 y 117')
    #     axs[0, 1].imshow(thresholded_region); axs[0, 1].grid(False); axs[0, 1].set_title('Threshold')
    #     axs[1, 0].imshow(label_region); axs[1, 0].grid(False); axs[1, 0].set_title('Label')
    #     plt.show()
    #
    # if x_center < 125 and x_center > 124 and abs(y_center) < 29 and abs(y_center) > 28:
    #     fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10, 10])
    #     axs[0, 0].imshow(region_of_int); axs[0, 0].grid(False); axs[0, 0].set_title('x 117 y 116')
    #     axs[0, 1].imshow(thresholded_region); axs[0, 1].grid(False); axs[0, 1].set_title('Threshold')
    #     axs[1, 0].imshow(label_region); axs[1, 0].grid(False); axs[1, 0].set_title('Label')
    #     plt.show()

    # if abs(x_center) < 123 and abs(x_center) > 122 and y_center < 118 and y_center > 117:
    #     fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10, 10])
    #     axs[0, 0].imshow(region_of_int); axs[0, 0].grid(False); axs[0, 0].set_title('x -122 y 117')
    #     axs[0, 1].imshow(thresholded_region); axs[0, 1].grid(False); axs[0, 1].set_title('Threshold')
    #     axs[1, 0].imshow(label_region); axs[1, 0].grid(False); axs[1, 0].set_title('Label')
    #     plt.show()
    #
    # if x_center < 118 and x_center > 117 and y_center < 117 and y_center > 116:
    #     fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10, 10])
    #     axs[0, 0].imshow(region_of_int); axs[0, 0].grid(False); axs[0, 0].set_title('x 117 y 116')
    #     axs[0, 1].imshow(thresholded_region); axs[0, 1].grid(False); axs[0, 1].set_title('Threshold')
    #     axs[1, 0].imshow(label_region); axs[1, 0].grid(False); axs[1, 0].set_title('Label')
    #     plt.show()
    #
    # if x_center < 98 and x_center > 97 and y_center < 117 and y_center > 116:
    #     fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10, 10])
    #     axs[0, 0].imshow(region_of_int); axs[0, 0].grid(False); axs[0, 0].set_title('x 97 y 116')
    #     axs[0, 1].imshow(thresholded_region); axs[0, 1].grid(False); axs[0, 1].set_title('Threshold')
    #     axs[1, 0].imshow(label_region); axs[1, 0].grid(False); axs[1, 0].set_title('Label')
    #     plt.show()
    #
    # if abs(x_center) < 103 and abs(x_center) > 102 and y_center < 117 and y_center > 116:
    #     fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10, 10])
    #     axs[0, 0].imshow(region_of_int); axs[0, 0].grid(False); axs[0, 0].set_title('x -102 y 116')
    #     axs[0, 1].imshow(thresholded_region); axs[0, 1].grid(False); axs[0, 1].set_title('Threshold')
    #     axs[1, 0].imshow(label_region); axs[1, 0].grid(False); axs[1, 0].set_title('Label')
    #     plt.show()
    #
    # if abs(x_center) < 83 and abs(x_center) > 82 and y_center < 116 and y_center > 115:
    #     fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10, 10])
    #     axs[0, 0].imshow(region_of_int); axs[0, 0].grid(False); axs[0, 0].set_title('x -82 y 115')
    #     axs[0, 1].imshow(thresholded_region); axs[0, 1].grid(False); axs[0, 1].set_title('Threshold')
    #     axs[1, 0].imshow(label_region); axs[1, 0].grid(False); axs[1, 0].set_title('Label')
    #     plt.show()

def mean_detection(pixels, x_center, y_center):
    # to not hard code this function we need to save mean value of each region and determine a cut-off mean value from the spread
    # by currently returning boolean and not saving mean values, we cannot do this
    pixels_mean = np.mean(pixels)

    print(f"Mean: {pixels_mean:.2f}, X location: {x_center:.2f}, Y location {y_center:.2f}")

    if pixels_mean < 360: # misses two true centroids
        return False
    else:
        return True

def percent_detection(thresholded_region, combined_ind, x_center, y_center):
    thresh_pixels = thresholded_region[combined_ind] # thresh pixels is a 1D array of the pixels within the thresholded ROI
    percent_nonzero = np.count_nonzero(thresh_pixels) / thresh_pixels.size # percent of pixels above the threshold

    print(f'Total: {thresh_pixels.size}, non-zero: {np.count_nonzero(thresh_pixels)}, percent: {percent_nonzero:.4f}, X location: {x_center:.2f}, Y location {y_center:.2f}')

    if percent_nonzero > 0.4:
        return True
    else:
        return False

def peak_detection(thresholded_region, x_center, y_center):
    label_region = label(thresholded_region)
    peaks = np.unique(label_region)[1:]  # disregard background

    number_of_peaks = len(peaks)
    print(f"Peaks: {number_of_peaks}, X location: {x_center:.2f}, Y location {y_center:.2f}")

    if number_of_peaks == 4:
        return True
    else:
        return False



# call on a single dicom image from the grid-based MRI folder
# dicom_path = r"C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\MR\1.3.46.670589.11.79127.5.0.6984.2022112517535313493.dcm"
# BW: I created a folder and copied just one image into it:
# 1.3.46.670589.11.79127.5.0.6984.2022112517535358579.dcm
dicom_path = Path(r"C:\Users\finmu\OneDrive\Documents\2023\Thesis - BMET4111 BMET4112\CODE\Grid-Based Sample Data\Test Slice")
#BW2: comment out the below to revert to your own data location
# dicom_path = Path(r'/home/brendan/Downloads/3Done phantom/'
                    # r'3DOne_zzphys_MR_2022-11-25_173121_post.upgrade_T1.3D.Tra.HN.L_n208__00000/')
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

# plot_2D_array(InputSlice, 'MRI Slice', extent)

# use prewit operators to find intersections:
v_edge_map = skimage.filters.prewitt_v(InputSlice)
# plot_2D_array(v_edge_map, 'Vertical Prewitt Operator', extent)

# highlight the gradient change down the vertical lines i.e., where the horizontal lines intersect
intersect_map = skimage.filters.prewitt_h(v_edge_map)
# plot_2D_array(intersect_map, 'Horizontal Prewitt Operator', extent)

# take absolute values
intersect_map = abs(intersect_map)
# plot_2D_array(intersect_map, 'Magnitude of Gradient Change Along Each Axis', extent)

# blurring the image points using a Guassian filter
blurred_map = gaussian_filter(intersect_map, sigma=1)
# plot_2D_array(blurred_map, 'Gaussian Filter', extent)

# otsu method to find the optimal threshold value
threshold = skimage.filters.threshold_otsu(blurred_map)


# apply the threshold
binary_image = blurred_map > threshold
# plot_2D_array(binary_image, 'Otsu\'s Method Thresholding', extent)

# BW: create a plot of our images so far:
#BW2 I will update the extent in one
# fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[10,10])
# axs[0, 0].imshow(InputSlice, extent=extent); axs[0, 0].grid(False); axs[0, 0].set_title('Dicom Image')
# axs[0, 1].imshow(intersect_map, extent=extent); axs[0, 1].grid(False); axs[0, 1].set_title('Prewit Operators')
# axs[1, 0].imshow(blurred_map, extent=extent); axs[1, 0].grid(False); axs[1, 0].set_title('Blurred')
# axs[1, 1].imshow(binary_image, extent=extent); axs[1, 1].grid(False); axs[1, 1].set_title('Otsu\'s threshold')
# plt.show()

# segment regions using label method
label_image = label(binary_image)
plot_2D_array(label_image, 'Labelled Image', extent)

# calculate the expected min/max volume of the centroids
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
    show_debug_image = False  # note will make a lot of images!!
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
plt.xlabel('X [mm]')
plt.ylabel('Y [mm]')
plt.title('Extracted Grid Marker Positions')
plt.show()



# centroids = np.array([x_centroids, y_centroids]).T


# scatter plot of the centroids
# plot_2D_scatter(x_centroids, y_centroids)
# overlay by updating