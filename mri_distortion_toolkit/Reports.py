import sys
import os
from pathlib import Path
import numpy as np
import logging
import plotly.express as px
import pandas as pd
from .utilities import get_harmonics, convert_cartesian_to_spherical, reconstruct_Bz
from .utilities import convert_spherical_to_cartesian
import plotly.graph_objects as go
from distutils.dir_util import copy_tree
from datetime import datetime
from jinja2 import Template
from .utilities import get_dicom_data

ch = logging.StreamHandler()
formatter = logging.Formatter('[%(filename)s: line %(lineno)d %(levelname)8s] %(message)s')
ch.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.addHandler(ch)
logger.setLevel(logging.INFO)  # This toggles all the logging in your app
logger.propagate = False


class TG284_Tests:  # pragma: no cover
    """
    tests defined here
    https://aapm.onlinelibrary.wiley.com/doi/full/10.1002/mp.14695
    """
    def test_distortion_less_than_1mm_within_r_100(self):
        """
        distortion test to be run by MRI_QA_Reporter
        """
        test_data = self._extract_data_from_MatchedMarkerVolume(r_max=100)
        if test_data.abs_dis.max() < 1:
            return True
        else:
            return False

    def test_distortion_less_than_2mm_within_r_200(self):
        """
        distortion test to be run by MRI_QA_Reporter
        """
        test_data = self._extract_data_from_MatchedMarkerVolume(r_max=200)
        if test_data.abs_dis.max() < 1:
            return True
        else:
            return False


class DefaultTestSuite:  # pragma: no cover
    """
    these are the tests which are run if no others are specified
    """

    def test_distortion_less_than_2mm_within_r_100(self):
        """
        distortion test to be run by MRI_QA_Reporter
        """
        test_data = self._extract_data_from_MatchedMarkerVolume(r_max=100)
        if test_data.abs_dis.max() < 2:
            return True
        else:
            return False

    def test_distortion_less_than_4mm_within_r_150(self):
        """
        distortion test to be run by MRI_QA_Reporter
        """
        test_data = self._extract_data_from_MatchedMarkerVolume(r_max=150)
        if test_data.abs_dis.max() < 4:
            return True
        else:
            return False


class Elekta_Distortion_tests:  # pragma: no cover
    """
    Geometric tests for Elekta system as described here:
    https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1002/mp.14764
    """

    def test_distortion_less_than_1_mm_within_r_100_mm(self):
        """
        distortion test to be run by MRI_QA_Reporter
        """
        test_data = self._extract_data_from_MatchedMarkerVolume(r_max=100)
        if test_data.abs_dis.max() < 1:
            return True
        else:
            return False

    def test_distortion_less_than_2_mm_within_r_150_mm(self):
        """
        distortion test to be run by MRI_QA_Reporter
        """
        test_data = self._extract_data_from_MatchedMarkerVolume(r_max=150)
        if test_data.abs_dis.max() < 2:
            return True
        else:
            return False

    def test_distortion_less_than_4_mm_within_r_200_mm(self):
        """
        distortion test to be run by MRI_QA_Reporter
        """
        test_data = self._extract_data_from_MatchedMarkerVolume(r_max=200)
        if test_data.abs_dis.max() < 4:
            return True
        else:
            return False


class MRI_QA_Reporter:
    """
    generate a report based on either some matched data, or harmonics. In the latter case, report data is reconstructed.

    :param MatchedMarkerVolume: pandas dataframe containing  ['x_gt', 'y_gt', 'z_gt', 'x_gnl', 'y_gnl', 'z_gnl']
        You must enter one of MatchedMarkerVolume or gradient_harmonics
    :type MatchedMarkerVolume: pandas.DataFrame, optional
    :param gradient_harmonics: list of either three pandas series or three paths to csv files
        from SphericalHarmonicFit.
        List elements should correspond to [G_x_harmonics, G_y_harmonics, G_z_harmonics].
        You must enter one of MatchedMarkerVolume or gradient_harmonics
    :type gradient_harmonics: list, optional
    :param B0_harmonics: pandas series from SphericalHarmonicFit or path to saved csv file. If entered, is used
        to report and plot B0 homogeneity
    :type B0_harmonics: pandas.Series or path, optional
    :param recon_coords: The coordinates you want to calculate the fields at. a pandas Dataframe which can contain
        whatever you want, but MUST contain columns called: x, y, z. If left as None a default is generated.
    :type recon_coords: pandas.DataFrame, optional
    :param tests_to_run: a class containing any tests to run. The class in this instance serves only as a name
        space. All methods inside the class starting with 'test' are identified as tests. The test methods will
        have access to everything inside this class.
    :type tests_to_run: class, optional
    :param dicom_data: either a dictionary or path to a saved json file of the dicom_data object generated by
        MarkerVolume. This data is required for harmonic reconstruction. For MatchedMarkerVolume it is not required,
        but it does contain useful information about the acquisition so still a good idea to use it.
    :type dicom_data: dict or path to json, optional
    :param r_outer: for the harmonic approach, r_outer is used to determine where the data is reconstructed. It is
        also used to control the limits in which data is plotted
    :type r_outer: float, optional
    :param show_plots: if True, all plots will be shown in your browser. Otherwise they are only saved
    :type show_plots: bool, optional
    :param style: control the report style. At the moment can only be 'light' or 'dark'
    :type style: str
    """

    def __init__(self, MatchedMarkerVolume=None, gradient_harmonics=None, B0_harmonics=None,
                 recon_coords_cartesian=None, tests_to_run=DefaultTestSuite, dicom_data=None,
                 r_outer=None, show_plots=False, style='dark'):

        self._html_template_loc = Path(__file__).parent.resolve() / 'jinja_templates' / 'MR_report.html'
        self._show_plots = show_plots
        _available_plotly_themes = ["plotly", "plotly_white", "plotly_dark", "ggplot2",
                                    "seaborn", "simple_white", "none"]
        self._style = style
        self._set_styles()
        self._jinja_dict = {}  # this will hold all the info jinja needs to render
        self._test_class = tests_to_run
        self.r_outer = r_outer
        self.recon_coords_cartesian = recon_coords_cartesian
        self.dicom_data = get_dicom_data(dicom_data)
        self._build_acquisition_dict()
        self._get_analysis_data(MatchedMarkerVolume, gradient_harmonics, B0_harmonics)
        self._update_MatchedMarkerVolume()

        # put plots here
        self._plot_distortion_v_r()
        self._plot_3D_cutplanes()
        self._plot_B0_surface()
        self._build_homogeneity_report_table()
        # self._plot_gradients_surface()
        self._get_test_results()

    def _build_acquisition_dict(self):
        """
        will try and extract acquisition data from dicom_data, and put it in a new dictionary
        if this new dictionary exists, we will include the acquisition data in the report
        """
        if self.dicom_data is None:
            return
        self._jinja_dict['acquisition_data'] = {}
        _fields_of_interest = ['acquisition_date', 'magnetic_field_strength', 'bandwidth', 'pixel_spacing',
                               'freq_encode_direction', 'manufacturer','sequence_name', 'imaging_frequency']
        for key in self.dicom_data:
            if key in _fields_of_interest:
                new_key = key.replace('_', ' ')
                self._jinja_dict['acquisition_data'][new_key] = self.dicom_data[key]

    def _set_styles(self):
        """
        sets the report style by changing the plotly theme and chanding the html css file
        """
        if self._style == 'dark':
            self._plotly_theme = 'plotly_dark'
            self._html_theme = "themes/d42ker-github.css"
        elif self._style == 'light':
            self._plotly_theme = 'plotly_white'
            self._html_theme = "themes/whitey.css"
        else:
            raise AttributeError(f'unknown style entered: {self._style}. allowed options are "light" or "dark"')

    def _unique_name_generator(self, save_loc, input_name, rel_path=True):
        """
        take input name, append the current date to it
        then check if this new name already exists in save_loc; if it does add integers to it (_1, _2, etc.) until it doesn't
        """
        if not save_loc.is_dir():
            # in this case, the code will likely fail downstream but we don't have to do anything here
            return input_name
        input_name, extension = os.path.splitext(input_name)
        now = datetime.now()  # current date and time
        date_string = now.strftime("%d_%m_%Y")
        new_name = input_name + '_' + date_string
        appended_int = 0
        while (save_loc / (new_name + extension)).is_file():
            appended_int = appended_int + 1
            new_name = new_name + '_' + str(appended_int)
        new_name = new_name + extension
        if rel_path:
            new_name = (save_loc  / new_name).relative_to(save_loc.parent)
        else:
            new_name = save_loc / new_name
        return new_name

    def _get_test_results(self):
        """
        runs through all tests found in the _test_class and keeps all the results in _jinja_dict
        """

        self._jinja_dict['test_collection'] = {}
        self._jinja_dict['TestSuiteName'] = self._test_class.__name__
        method_list = [method for method in dir(self._test_class) if method.startswith('__') is False]
        # this will list all methods except 'dunder' methods e.g. __init__
        for test in method_list:
            if test[:4] == 'test':
                # then we found a true test
                test_name = test.replace('_', ' ') # under scores to space
                test_name = test_name.replace('test ','' )  # remove test from the name
                test_string = getattr(self._test_class, test)(self)
                if isinstance(test_string, bool):
                    if test_string:
                        test_string = "<span style='color:green'>Pass</span>"
                    else:
                        test_string = "<span style='color:red'>Fail</span>"
                elif not isinstance(test_string, str):
                    raise TypeError('please return test results as either booleans or strings')
                self._jinja_dict['test_collection'][test_name] = test_string

    def _get_analysis_data(self, MatchedMarkerVolume, gradient_harmonics, B0_harmonics):
        """
        check what the user has input, and attach to self
        if necessary, convert gradient_harmonics to distortion map
        """

        if MatchedMarkerVolume is not None and gradient_harmonics is not None:
            raise AttributeError(
                f'you cannot enter both MatchedMarkerVolume and gradient_harmonics at the same time...')

        if MatchedMarkerVolume is not None:
            self._check_marker_volume(MatchedMarkerVolume)
            self._MatchedMarkerVolume = MatchedMarkerVolume.copy()

        if gradient_harmonics is not None:
            self._check_dicom_data()
            self._Gx_Harmonics, self._Gy_Harmonics, self._Gz_Harmonics, self._B0_harmonics = \
                get_harmonics(gradient_harmonics[0], gradient_harmonics[1], gradient_harmonics[2])
            self._Gx_Harmonics = self._Gx_Harmonics * self.dicom_data['gradient_strength'][0]
            self._Gy_Harmonics = self._Gy_Harmonics * self.dicom_data['gradient_strength'][1]
            self._Gz_Harmonics = self._Gz_Harmonics * self.dicom_data['gradient_strength'][2]
            if self.r_outer is None:
                logger.warning('no r_outer value entered, using 150 mm')
                self.r_outer = 150
            if self.recon_coords_cartesian is None:
                self.recon_coords_cartesian = self._generate_recon_coords()
            self.recon_coords = convert_cartesian_to_spherical(self.recon_coords_cartesian)
            self._reconstruct_gradient_fields()
            self._convert_magnetic_field_to_distortion()

        if B0_harmonics is not None:
            if isinstance(B0_harmonics, pd.Series):
                self.B0_harmonics = B0_harmonics
            elif isinstance(B0_harmonics, (str, Path)):
                self.B0_harmonics = pd.read_csv(B0_harmonics, index_col=0).squeeze("columns")
            else:
                raise AttributeError('could not read in B0_harmonics...please input either a series or a '
                                     'path to a saved csv file')
        else:
            self.B0_harmonics = None

    def _convert_magnetic_field_to_distortion(self):
        """
        converts reconstructed gradient fields into geometric distortion.
        Essentially this process is doing the opposite of FieldCalculation.ConvertMatchedMarkersToBz
        """
        self._MatchedMarkerVolume = pd.DataFrame()
        self._MatchedMarkerVolume = self._MatchedMarkerVolume.assign(x_gt=self.recon_coords_cartesian.x)
        self._MatchedMarkerVolume = self._MatchedMarkerVolume.assign(y_gt=self.recon_coords_cartesian.y)
        self._MatchedMarkerVolume = self._MatchedMarkerVolume.assign(z_gt=self.recon_coords_cartesian.z)
        bandwidth = np.array(self.dicom_data['bandwidth'])
        image_size = np.array(self.dicom_data['image_size'])
        gama = np.array(self.dicom_data['gama'])
        FOV = np.array(self.dicom_data['FOV'])

        gradient_strength = bandwidth * image_size / (gama * 1e6 * FOV * 1e-3)  # unit(T / m)
        # ^ this is a vector [gx, gy, gz]
        # factor of 1e3 is for m to mm
        self._MatchedMarkerVolume = \
            self._MatchedMarkerVolume.assign(x_gnl=self.Gx_Bfield * 1e3 / (gradient_strength[0]))
        self._MatchedMarkerVolume = \
            self._MatchedMarkerVolume.assign(y_gnl=self.Gy_Bfield * 1e3 / (gradient_strength[1]))
        self._MatchedMarkerVolume = \
            self._MatchedMarkerVolume.assign(z_gnl=self.Gz_Bfield * 1e3 / (gradient_strength[2]))

    def _reconstruct_gradient_fields(self):
        """
        utilise the utilities.reconstruct_Bz function to reconsutrct each of the gradient fields at recon_coords
        """
        self.Gx_Bfield = reconstruct_Bz(self._Gx_Harmonics, self.recon_coords, quantity='T')
        self.Gy_Bfield = reconstruct_Bz(self._Gy_Harmonics, self.recon_coords, quantity='T')
        self.Gz_Bfield = reconstruct_Bz(self._Gz_Harmonics, self.recon_coords, quantity='T')

    def _build_homogeneity_report_table(self):
        """
        Generate a table of predicted pk-pk homogeneity at different r
        :return:
        """
        if self.B0_harmonics is None:
            return
        radii_of_interest = [50, 100, 150, 200]
        pk_pk = []
        for DSV_radius in radii_of_interest:
            azimuth = np.linspace(0, 2 * np.pi, 100)
            elevation = np.linspace(0, np.pi, 100)
            [AZ, EL, R] = np.meshgrid(azimuth, elevation, DSV_radius, indexing='ij')
            surface_coords = pd.DataFrame({'r': R.flatten(), 'azimuth': AZ.flatten(), 'elevation': EL.flatten()})
            surface_coords = convert_spherical_to_cartesian(surface_coords)
            B0 = reconstruct_Bz(self.B0_harmonics, surface_coords, quantity='T', r_outer=DSV_radius)
            pk_pk.append(int((B0.max() - B0.min()) * 1e6))
        _B0_stats = pd.DataFrame({'r': radii_of_interest, 'pk_pk [\u03BCT]': pk_pk} )

        self._B0_table = go.Figure(data=[go.Table(
            header=dict(values=list(_B0_stats.columns),
                        align='left'),
            cells=dict(values=[_B0_stats.r, _B0_stats['pk_pk [\u03BCT]']],
                       align='left'))
        ])
        self._B0_table.update_layout(template=self._plotly_theme)

    def _check_dicom_data(self):
        """
        check that if input dicom data is required, it confirms to minimum standard
        """
        _required_fields = ['FOV', 'bandwidth', 'gama', 'image_size', 'gradient_strength']

        if not isinstance(self.dicom_data, dict) or not all(
                name in self.dicom_data.keys() for name in _required_fields):
            raise AttributeError(f'dicom_data must be a dict with at least the following keys: {_required_fields}')

    def _check_marker_volume(self, MatchedMarkerVolume):
        """
        check that if a MatchedMarkerVolume is entered, is meets expected format
        """
        _required_cols = ['x_gt', 'y_gt', 'z_gt', 'x_gnl', 'y_gnl', 'z_gnl']
        if not isinstance(MatchedMarkerVolume, pd.DataFrame) or \
                not all(name in MatchedMarkerVolume.columns for name in _required_cols):
            raise AttributeError(f'MatchedMarkerVolume must be a dataframe containing columns'
                                 f'[x_gt, y_gt, z_gt, x_gnl, y_gnl, z_gnl]')

    def _generate_recon_coords(self):
        """
        Generate a grid of coordinates to perform reconstruction
        """
        x = np.linspace(-self.r_outer, self.r_outer, 50)
        y = np.linspace(-self.r_outer, self.r_outer, 50)
        z = np.linspace(-self.r_outer, self.r_outer, 50)
        [x_recon, y_recon, z_recon] = np.meshgrid(x, y, z, indexing='ij')
        recon_coords = pd.DataFrame({'x': x_recon.flatten(), 'y': y_recon.flatten(), 'z': z_recon.flatten()})
        return recon_coords

    def _plot_distortion_v_r(self):
        """
        generate a histogram as a function of r
        """

        # data over different planes
        xy_plot_data = self._extract_data_from_MatchedMarkerVolume(r_max=self.r_outer, z_select=0)
        xy_plot_data['plane'] = 'XY'
        zx_plot_data = self._extract_data_from_MatchedMarkerVolume(r_max=self.r_outer, y_select=0)
        zx_plot_data['plane'] = 'ZX'
        zy_plot_data = self._extract_data_from_MatchedMarkerVolume(r_max=self.r_outer, x_select=0)
        zy_plot_data['plane'] = 'ZY'
        # create a new data frame combined from these for easy plotting
        all_dsv_plot_data = self._extract_data_from_MatchedMarkerVolume(r_max=self.r_outer)
        all_dsv_plot_data['plane'] = 'Entire DSV'


        plot_data = pd.concat([all_dsv_plot_data, xy_plot_data, zx_plot_data, zy_plot_data], ignore_index=True, axis=0)
        self._fig_distortion_v_r = px.scatter(plot_data, x="r_gt", y="abs_dis",color="plane",
                                              labels={
                                                  "r_gt": "Distance from the center (mm)", #rename the x axis
                                                  "abs_dis": " Absolute Distortion (mm)", #rename the y axis
                                              },
                                              title=f"Cardinal plane data, {self.r_outer} DSV")


        self._fig_distortion_v_r.update_layout(template=self._plotly_theme)
        if self._show_plots:
            self._fig_distortion_v_r.show()

    def _plot_3D_cutplanes(self):
        """
        create a 3D cone plot in each of the three cardinal directions
        """
        xy_plot_data = self._extract_data_from_MatchedMarkerVolume(r_max=self.r_outer, z_select=0)
        zx_plot_data = self._extract_data_from_MatchedMarkerVolume(r_max=self.r_outer, y_select=0)
        zy_plot_data = self._extract_data_from_MatchedMarkerVolume(r_max=self.r_outer, x_select=0)
        plot_data = pd.concat([xy_plot_data, zx_plot_data, zy_plot_data], ignore_index=True, axis=0)

        self._fig_3D_planes = go.Figure(data=go.Cone(
            x=plot_data['x_gt'],
            y=plot_data['y_gt'],
            z=plot_data['z_gt'],
            u=plot_data['x_dis'],
            v=plot_data['y_dis'],
            w=plot_data['z_dis'],
            colorbar=dict(title= 'mm'),
            sizemode="absolute",
            sizeref=1))


        self._fig_3D_planes.update_layout(scene=dict(aspectratio=dict(x=1, y=1, z=0.8),
                                                     camera_eye=dict(x=1.2, y=1.2, z=0.6)))
        self._fig_3D_planes.update_layout(template=self._plotly_theme, scene = dict(
                    xaxis_title='X Axis (mm)',
                    yaxis_title='Y Axis (mm)',
                    zaxis_title='Z Axis (mm)'))
        if self._show_plots:
            self._fig_3D_planes.show()

    def _plot_B0_surface(self):
        """
        if self.B0_harmonics exists, reconstruct B0 on the surface of r_outer and make a pretty 3D plot
        """

        if self.B0_harmonics is None:
            return
        # generate surface recon coordinates
        if self.r_outer is not None:
            r = self.r_outer
        else:
            r = 150
        azimuth = np.linspace(0, 2 * np.pi, 100)
        elevation = np.linspace(0, np.pi, 100)
        [AZ, EL, R] = np.meshgrid(azimuth, elevation, r, indexing='ij')
        surface_coords = pd.DataFrame({'r': R.flatten(), 'azimuth': AZ.flatten(), 'elevation': EL.flatten()})
        surface_coords = convert_spherical_to_cartesian(surface_coords)

        B0 = reconstruct_Bz(self.B0_harmonics, surface_coords, quantity='uT', r_outer=r)
        x = surface_coords.x.to_numpy().reshape(100, 100)
        y = surface_coords.y.to_numpy().reshape(100, 100)
        z = surface_coords.z.to_numpy().reshape(100, 100)
        B0 = B0.to_numpy().reshape(100, 100)

        self._fig_DSV_surface = go.Figure(data=[go.Surface(z=z, x=x, y=y,colorbar=dict(title= 'B(uT)'), surfacecolor=B0)])

        self._fig_DSV_surface.update_layout(title='DSV surface [uT]', autosize=True,
                          template=self._plotly_theme, scene = dict(
                    xaxis_title='X Axis (mm)',
                    yaxis_title='Y Axis (mm)',
                    zaxis_title='Z Axis (mm)') )

        if self._show_plots:
            self._fig_DSV_surface.show()

    def _update_MatchedMarkerVolume(self):
        """
        adds in 'r_gt', 'r_gnl' and 'dis_x,y,z' to _MatchedMarkerVolume input
        """
        self._MatchedMarkerVolume['r_gt'] = np.sqrt(self._MatchedMarkerVolume.x_gt ** 2
                                                    + self._MatchedMarkerVolume.y_gt ** 2
                                                    + self._MatchedMarkerVolume.z_gt ** 2)
        self._MatchedMarkerVolume['r_gnl'] = np.sqrt(self._MatchedMarkerVolume.x_gnl ** 2
                                                     + self._MatchedMarkerVolume.y_gnl ** 2
                                                     + self._MatchedMarkerVolume.z_gnl ** 2)
        self._MatchedMarkerVolume['x_dis'] = self._MatchedMarkerVolume.x_gnl - self._MatchedMarkerVolume.x_gt
        self._MatchedMarkerVolume['y_dis'] = self._MatchedMarkerVolume.y_gnl - self._MatchedMarkerVolume.y_gt
        self._MatchedMarkerVolume['z_dis'] = self._MatchedMarkerVolume.z_gnl - self._MatchedMarkerVolume.z_gt
        self._MatchedMarkerVolume['abs_dis'] = np.sqrt(self._MatchedMarkerVolume.x_dis ** 2
                                                       + self._MatchedMarkerVolume.y_dis ** 2
                                                       + self._MatchedMarkerVolume.z_dis ** 2)
        # get resolution in each direction
        self._x_res = np.mean(np.diff(np.sort(self._MatchedMarkerVolume.x_gt.unique())))
        self._y_res = np.mean(np.diff(np.sort(self._MatchedMarkerVolume.y_gt.unique())))
        self._z_res = np.mean(np.diff(np.sort(self._MatchedMarkerVolume.z_gt.unique())))

    def _extract_data_from_MatchedMarkerVolume(self, x_select=None, y_select=None, z_select=None,
                                               max_num_points=None, r_max=None):
        """
        Extract data from the matched volume for use in direct plotting

        uses the values of _x_res etc. to define a tolerance
        """
        if (x_select is None) and (y_select is None) and (z_select is None) and (r_max is None):
            logger.warning('you should set a value for at least one of x_ind, y_ind, and z_ind or r_max')
        # start with everything selected, overwrite below
        x_ind = np.ones(self._MatchedMarkerVolume.shape[0])
        y_ind = np.ones(self._MatchedMarkerVolume.shape[0])
        z_ind = np.ones(self._MatchedMarkerVolume.shape[0])
        r_ind = np.ones(self._MatchedMarkerVolume.shape[0])

        if x_select is not None:
            x_ind = np.logical_and(self._MatchedMarkerVolume.x_gt <= (x_select + self._x_res / 2),
                                   self._MatchedMarkerVolume.x_gt >= (x_select - self._x_res / 2))
        if y_select is not None:
            y_ind = np.logical_and(self._MatchedMarkerVolume.y_gt <= (y_select + self._y_res / 2),
                                   self._MatchedMarkerVolume.y_gt >= (y_select - self._y_res / 2))
        if z_select is not None:
            z_ind = np.logical_and(self._MatchedMarkerVolume.z_gt <= (z_select + self._z_res / 2),
                                   self._MatchedMarkerVolume.z_gt >= (z_select - self._z_res / 2))
        if r_max is not None:
            r_ind = self._MatchedMarkerVolume.r_gt <= r_max

        select_ind = np.logical_and(np.logical_and(x_ind, y_ind), z_ind)
        select_ind = np.logical_and(select_ind, r_ind)
        plot_data = self._MatchedMarkerVolume[select_ind]
        plot_data = plot_data.reset_index(drop=True)
        if (max_num_points is not None) and plot_data.shape[0] > max_num_points:
            every_nth_point = int(plot_data.shape[0] / max_num_points)
            plot_data = plot_data[plot_data.reset_index().index % every_nth_point == 0]

        return plot_data

    def _set_up_directory_structure(self):
        """
        ensure we have a place to save the reports.
        will set up a directory in
        """
        if not os.path.isdir(self.output_folder):
            os.makedirs(self.output_folder)
        if not os.path.isdir(self.output_folder / 'plots'):
            os.mkdir(self.output_folder / 'plots')

        # copy themes
        theme_dir = (Path(__file__).parent / 'jinja_templates' / 'themes').resolve()
        copy_tree(str(theme_dir), str(self.output_folder / 'themes'))

    def _get_template(self):
        """
        loads the jinja template
        """

        with open(self._html_template_loc) as f:
            template_text = f.read()
        j2_template = Template(template_text)
        return j2_template

    # public methods

    def write_html_report(self, output_folder=None, report_name=None):
        """
        Generates a html report encompassing available acquisition information, test results, and intercative
        plotly figures

        :param output_folder: folder to write report. Defaults to ~/Documents/MR_QA_Reports
        :type output_folder: Path or string, optional
        """
        if output_folder is None:
            self.output_folder = Path(os.path.expanduser('~')) / 'Documents' / 'MR_QA_Reports'
        else:
            self.output_folder = Path(output_folder)
        self._set_up_directory_structure()

        if report_name is None:
            report_name = 'MR_QA_report.html'
        else:
            report_name = str(Path(report_name).with_suffix('.html'))

        # save plots and update jinja_dict

        distortion_v_r_save_name = self._unique_name_generator(self.output_folder / 'plots', 'distortion_v_r.html')
        self._fig_distortion_v_r.write_html(self.output_folder / distortion_v_r_save_name, full_html=False, include_plotlyjs='cdn')
        self._jinja_dict['dist_v_r_source'] = distortion_v_r_save_name

        threeD_plane_save_name = self._unique_name_generator(self.output_folder / 'plots', '3D_planes.html')
        self._fig_3D_planes.write_html(self.output_folder / threeD_plane_save_name, full_html=False, include_plotlyjs='cdn')
        self._jinja_dict['cutplanes_source'] = threeD_plane_save_name

        if self.B0_harmonics is not None:
            DSV_surface_save_name = self._unique_name_generator(self.output_folder / 'plots', 'DSV_surface.html',)
            self._fig_DSV_surface.write_html(self.output_folder / DSV_surface_save_name, full_html=False, include_plotlyjs='cdn')
            self._jinja_dict['dsv_surf_source'] = DSV_surface_save_name

            B0_table_save_name = self._unique_name_generator(self.output_folder / 'plots', 'B0_table.html')
            self._B0_table.write_html(self.output_folder / B0_table_save_name, full_html=False, include_plotlyjs='cdn')
            self._jinja_dict['B0_table_source'] = B0_table_save_name
        else:
            self._jinja_dict['dsv_surf_source'] = None

        # set up template
        self._jinja_dict['html_theme_loc'] = self._html_theme
        j2_template = self._get_template()
        report_name = self._unique_name_generator(self.output_folder, report_name, rel_path=False)
        report_code = os.path.split(report_name)[1]  # for printing
        report_code = os.path.splitext(report_code)[0]
        report_code = report_code.replace('MR_QA_report_','')
        self._jinja_dict['report_name'] = report_code
        report_string = j2_template.render(self._jinja_dict)
        with open(report_name, 'w') as f:
            f.write(report_string)

        print(f'The report has been compiled and can be found at {self.output_folder}')
