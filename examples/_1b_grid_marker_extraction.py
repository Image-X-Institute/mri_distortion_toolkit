import numpy as np
from pathlib import Path
from mri_distortion_toolkit.utilities import dicom_to_numpy
import skimage.filters
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter
from skimage.measure import label
import pandas as pd
from mri_distortion_toolkit.MarkerAnalysis import MarkerVolume, MatchedMarkerVolumes
from mri_distortion_toolkit import utilities as ut

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

def get_marker_max_min_volume(labels, marker_size_lower_tol=0.80, marker_size_upper_tol=5):
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
    plt.figure()
    plt.grid(False)
    plt.imshow(image_to_plot, extent=extent)
    plt.xlim(-180, 175)
    plt.ylim(150, -90)
    plt.xlabel('X [mm]')
    plt.ylabel('Y [mm]')
    if title:
        plt.title(title)
    plt.colorbar()
    plt.show()

def plot_2D_scatter(x_centroids, y_centroids):
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
    # we want to detect pixels where:
    # centroid_x - r_roi <= centroid_x <= centroid_x <= centroid_x + r_roi
    # and similar for y
    x_ind = np.logical_and((X - r_roi <= x_center), (X+r_roi >= x_center))
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

    # apply a test to pixels to see if we should count it or not:
    # average pixel value
    proper_centroid1 = mean_detection(pixels, x_center, y_center)

    # percentage of non zero pixels in ROI after thresholding
    proper_centroid2 = percent_detection(thresholded_region, combined_ind)

    # number of peak in a histogram
    proper_centroid3 = peak_detection(thresholded_region)
    # the intersection points should have a peak in each corner whereas edge cases should be missing at least one peak

    return proper_centroid1 and proper_centroid2 and proper_centroid3


def mean_detection(pixels, x_center, y_center):
    # to not hard code this function we need to save mean value of each region and determine a cut-off mean value from the spread
    # by currently returning boolean and not saving mean values, we cannot do this
    pixels_mean = np.mean(pixels)

    # print(f"Mean: {pixels_mean:.2f}, X location: {x_center:.2f}, Y location {y_center:.2f}")

    if y_center > 110:
        return False

    if pixels_mean < 420:
        return False
    else:
        return True

def percent_detection(thresholded_region, combined_ind):
    thresh_pixels = thresholded_region[combined_ind] # thresh pixels is a 1D array of the pixels within the thresholded ROI
    percent_nonzero = np.count_nonzero(thresh_pixels) / thresh_pixels.size # percent of pixels above the threshold

    # print(f'Total: {thresh_pixels.size}, non-zero: {np.count_nonzero(thresh_pixels)}, percent: {percent_nonzero:.4f}, X location: {x_center:.2f}, Y location {y_center:.2f}')

    if percent_nonzero > 0.4:
        return True
    else:
        return False

def peak_detection(thresholded_region):
    label_region = label(thresholded_region)
    peaks = np.unique(label_region)[1:]  # disregard background

    number_of_peaks = len(peaks)
    # print(f"Peaks: {number_of_peaks}, X location: {x_center:.2f}, Y location {y_center:.2f}")

    if number_of_peaks >= 4:
        return True
    else:
        return False

def generate_ground_truth_marker_volume(x_start=0.0, y_start=0.0, z_start=0.0, spacing=20,
                                        n_markers_x=25, n_markers_y=25, n_markers_z=1):
    """
    generate a ground truth marker volume based on generating points around a single
    'seed' point (x_start, y_start, z_start). This point should be close to the center of the volume.
    Note that you can safely generate more ground truth markers than distorted markers, because the matching
    code is smart enough to discard ground truth that doesn't match to any marker

    :param x_start: x coodinate of seed point
    :param y_start: y coodinate of seed point
    :param z_start: z coodinate of seed point
    :param spacing: spacing between markers, assumed to be equal in all directions
    :param n_markers_x: number of markers in x direction.
    :param n_markers_y: number of markers in y direction.
    :param n_markers_z: number of markers in z direction.
    :return: gt_volume
    """
    # generate start and stop point in each direction
    start = []
    stop = []
    for n_markers in [n_markers_x, n_markers_y, n_markers_z]:
        if n_markers % 2 == 0:
            start.append(-1 * spacing * np.floor((n_markers / 2) - 1))
            stop.append(spacing * int(n_markers / 2))
        else:
            start.append(-1 * spacing * np.floor((n_markers / 2)))
            stop.append(spacing * np.floor((n_markers / 2)))

    # first will generate the points centered on zero
    x_points = np.linspace(start[0], stop[0], n_markers_x) + x_start
    y_points = np.linspace(start[1], stop[1], n_markers_y) + y_start
    z_points = np.linspace(start[2], stop[2], n_markers_z) + z_start
    # now use meshgrid to create an array with all the combined points
    [xx, yy, zz] = np.meshgrid(x_points, y_points, z_points)
    points_dataframe = pd.DataFrame({'x': xx.flatten(), 'y': yy.flatten(), 'z': zz.flatten()})
    gt_volume = MarkerVolume(points_dataframe)

    return gt_volume

# call on a single dicom image from the grid-based MRI folder
dicom_path = Path(r"Your_Folder")
InputVolume, dicom_affine, (X, Y, Z) = dicom_to_numpy(dicom_path,
                                                      FilesToReadIn= '1.3.46.670589.11.79127.5.0.6984.2022112517535397646.dcm',
                                                      file_extension='dcm',return_XYZ=True)
# to calculate extent of image:
extent = (X.min(), X.max(), Y.max(), Y.min())
# we can remove this singleton dimension like this:
InputSlice = InputVolume.squeeze()  # dont delete this line!
X = X.squeeze()
Y = Y.squeeze()
# plot_2D_array(InputSlice, 'MRI Slice', extent)

# use prewit operators to find intersections:
v_edge_map = skimage.filters.prewitt_v(InputSlice)
#plot_2D_array(v_edge_map, 'Vertical Prewitt Operator', extent)

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

# segment regions using label method
label_image = label(binary_image)
# plot_2D_array(label_image, 'Labelled Image', extent)

# calculate the expected min/max volume of the centroids
voxel_min, voxel_max = get_marker_max_min_volume(label_image)

# find centroids
x_centroids = []
y_centroids = []
skipped = 0

unique_labels = np.unique(label_image)[1:]  # first label is background so skip
labels_to_remove = [] # for labelled sections we do not want to extract
# print('Number of labelled regions: ' + str(len(unique_labels)))

# extract x,y of each connected region
for i, label_level in enumerate(unique_labels):
    RegionInd = np.equal(label_image, label_level)
    voxels = np.count_nonzero(RegionInd)
    if voxels < voxel_min or voxels > voxel_max:
        skipped += 1
        labels_to_remove.append(label_level) # exclude labelled regions outside expected sizes
        continue # skip outliers

    # centroid = mean of coordinate of region.
    potential_x_centroid = np.mean(X[RegionInd])
    potential_y_centroid = np.mean(Y[RegionInd])
    # detect and remove edge cases:
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

    if not proper_centroid:
        labels_to_remove.append(label_level)

# iterate through the labeled_image and set specified labels to zero
for label in labels_to_remove:
    label_image[label_image == label] = 0

# plot thresholded image without edge cases
refined_binary_image = label_image > 0
# plot_2D_array(refined_binary_image, 'Removed Edge-Cases', extent)

# generate'ground truth' marker volume
gt_volume = generate_ground_truth_marker_volume(x_start=-2.8, y_start=28.7, z_start=np.mean(Z),
                                                n_markers_x=25, n_markers_y=25, n_markers_z=1,
                                                spacing=20)
# generate 'distorted marker volume
distorted_data_frame = pd.DataFrame({'x': x_centroids, 'y': y_centroids, 'z': np.mean(Z)})
distorted_volume = MarkerVolume(distorted_data_frame)

# match them together
matched_volume = MatchedMarkerVolumes(gt_volume, distorted_volume)
# a few plots
matched_volume.plot_compressed_markers(z_min=-90, z_max= -70)
matched_volume.plot_3D_markers()
ut.plot_distortion_xyz_hist(matched_volume) # histogram of distibution of distortion in each axis
ut.plot_matched_volume_hist([matched_volume]) # histogram of distribution of distortion magnitude

df = matched_volume.MatchedCentroids.match_distance
# Calculate max, min, mean, and standard deviation of the total geometric distortion in the image
df = df[df <= 10]

max_value = df.max()
min_value = df.min()
mean_value = df.mean()
std_deviation = df.std()

# Print or use the calculated values
print("Max:", max_value)
print("Min:", min_value)
print("Mean:", mean_value)
print("Standard Deviation:", std_deviation)