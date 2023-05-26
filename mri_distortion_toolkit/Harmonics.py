from pathlib import Path
import sys, os
import numpy as np
import pandas as pd
import logging
from matplotlib import pyplot as plt
from .utilities import bcolors
from .utilities import convert_cartesian_to_spherical, generate_legendre_basis
from .utilities import reconstruct_Bz
from .utilities import generate_harmonic_names
import seaborn as sns

sns.set_theme(style="whitegrid")
ch = logging.StreamHandler()
formatter = logging.Formatter('[%(filename)s: line %(lineno)d %(levelname)8s] %(message)s')
ch.setFormatter(formatter)
logger = logging.getLogger(__name__)
logger.addHandler(ch)
logger.setLevel(logging.INFO)  # This toggles all the logging in your app
logger.propagate = False


class SphericalHarmonicFit:
    """
    Uses single value decomposition to fit spherical harmonics to InputData.

    Let the legendre basis be called L, the Harmonics be H and the magnetic field B.
    L will have size [n_coords, n_harmonics], Harmonics has size [n_harmonics, 1] and B has size [n_coords, 1]

    Our proposition is that there exists some set of harmonics such that

    .. math:: L.H = B

    We can find H by inverting L:

    .. math::

        L^{-1}.L.H = L^{-1}.B \\
        H = L^{-1}.B

    This task is performed using numpys pseudo_inverse functionality to invert L

    :param _input_Bz_data: A pandas dataframe with columns x, y, z, Bz. Data should be in mm and T
    :type _input_Bz_data: Pandas.DataFrame
    :param r_outer: radius of sphere of interest. If AssessHarmonicPk_Pk=True, data is reconstucted on r_outer. if
        TrimDataBy_r_outer=True, data lying outside r_outer will be deleted. Finally, r_outer is used to set the limits
        on plotting data.
    :type r_outer: float, optional
    :param n_order: order of harmonic fit. higher order means more terms to fit, which may result in a better fit
        but also is slower. If the data you input doesn't satisfy some nyquist criterion (I need to figure this out!)
        then you will get a spurious result. For most MRI application 8 should be a reasonable number.
    :type n_order: int, optional
    :param AssessHarmonicPk_Pk: if True, data is reconstructed on the surface of r_outer, and the pk-pk contribution
        from each harmonic calculated. This is a very useful way to interpret harmonics which are otherwise quite
        abstract.
    :type AssessHarmonicPk_Pk: bool, optional
    :param QuantifyFit: if True, reconstruct data at all input data points, and compare input data point to
        reconstructed data point. This is useful to check how well the fit has worked.
    :type QuantifyFit: bool, optional
    :param TrimDataBy_r_outer: if True, any input data outside r_outer is deleted
    :type TrimDataBy_r_outer: bool, optional
    :param scale: Value to scale the harmonics by. The main intention of this parameter is to normalise gradient
        field harmonics to gradient strength. harmonics will be multiplied by the value of scale
    :type scale: float, optional
    """

    def __init__(self, input_Bz_data, r_outer=150, n_order=5, AssessHarmonicPk_Pk=True, QuantifyFit=True,
                 TrimDataBy_r_outer=False, scale=1):

        # attributes:
        self._tol = 1e-3  # tolerance for coordinate filtering in mm
        self.scale = scale
        self.TrimDataBy_r_outer = TrimDataBy_r_outer
        self.r_outer = r_outer
        self.QuantifyFit = QuantifyFit
        self._AssessHarmonicPk_Pk = AssessHarmonicPk_Pk
        self.n_order = n_order
        self.n_harmonics = (self.n_order + 1) ** 2
        self._input_Bz_data = input_Bz_data
        self.__recon_warning_thrown = False
        self._data_type = None
        self._determine_Bz_data_type()
        # basic fitting routine:
        if self._data_type == 'field_map':
            self._input_Bz_data = convert_cartesian_to_spherical(self._input_Bz_data)
            self._check_data_input()
            if self.TrimDataBy_r_outer:
                self._filter_data()
            self._coeff_names = generate_harmonic_names(self.n_order)
            self.legendre_basis = generate_legendre_basis(self._input_Bz_data, self.n_order)
            self._svd_fit()
            self.harmonics = self.harmonics * scale
        elif self._data_type == 'harmonics':
            r = self.r_outer
            azimuth = np.linspace(0, 2 * np.pi, 100)
            elevation = np.linspace(0, np.pi, 100)
            [r, azimuth, elevation] = np.meshgrid(r, azimuth, elevation)
            coords = pd.DataFrame({'r': r.flatten(), 'azimuth': azimuth.flatten(), 'elevation': elevation.flatten()})
            self.legendre_basis = generate_legendre_basis(coords, self.n_order)
            self.harmonics = self._input_Bz_data
        else:
            raise AttributeError('invalid input data')
        if self.QuantifyFit and self._data_type == 'field_map':
            self._quantify_fit()
        if self._AssessHarmonicPk_Pk:
            self._assess_harmonic_pk_pk()

    def _determine_Bz_data_type(self):
        """
        check whether we have field data, or harmonics
        """
        if isinstance(self._input_Bz_data, pd.Series):
            if not len(self._input_Bz_data) == (self.n_order + 1) ** 2:
                raise ValueError('length of series does not match input n_order')
            else:
                self._data_type = 'harmonics'
        elif isinstance(self._input_Bz_data, pd.DataFrame):
            self._data_type = 'field_map'  # not error checks in place downstream

    def _check_data_input(self):
        """
        - Make sure the input data actually covers a sphere. If it doesn't, it normally means something has or will go wrong
            so a warning is triggered.
        - Make sure that input data is in mm, which is an implicit assumption
        """

        if not isinstance(self.n_order, int):
            logger.warning(f'n_order should be integer; round {self.n_order} to {int(self.n_order)}')
            self.n_order = int(self.n_order)

        if self._data_type == 'field_map':
            if self._input_Bz_data.r.mean() < 10:
                logger.warning('it appears that your input data is in m, not mm - please use mm!')

            if (self._input_Bz_data.elevation.max() - self._input_Bz_data.elevation.min()) < 0.75 * np.pi or \
                    (self._input_Bz_data.azimuth.max() - self._input_Bz_data.azimuth.min()) < 1.8 * np.pi:
                logger.warning('input sample points do not appear to cover a full sphere')
            field_map_headings = ['x', 'y', 'z', 'Bz']
            for col in field_map_headings:
                if not col in self._input_Bz_data.columns:
                    raise AttributeError(f'Field map data must contain {field_map_headings}')

    def _filter_data(self):
        """
        filter the Bz data by radial coordinate, removing any entries that fall outside self.r_outer
        """

        data_to_delete_ind = self._input_Bz_data.r > (self.r_outer + self._tol)
        self._input_Bz_data = self._input_Bz_data.drop(self._input_Bz_data[data_to_delete_ind].index)
        self._input_Bz_data.reset_index(inplace=True)
        n_deleted = np.count_nonzero(data_to_delete_ind)
        if n_deleted > 0:
            logger.warning(f'deleting {n_deleted} of the original {data_to_delete_ind.size} data points because'
                           f'they fall outside the r_outer value of {self.r_outer} and TrimDataBy_r_outer=True')

        if self._input_Bz_data.shape[0] == 0:
            logger.error(f'After filtering for data insdide {self.r_outer} mm, there is no data left!')
            sys.exit(1)

    def _svd_fit(self):
        """
        calculate harmonics using single value decomposition of legendre_basis
        """

        inverseLegendre = np.linalg.pinv(self.legendre_basis)

        if not np.allclose(inverseLegendre @ self.legendre_basis,
                           np.eye(inverseLegendre.shape[0], self.legendre_basis.shape[1]),
                           rtol=5e-3, atol=5e-4):
            # check that the inverse times itself returns an idendity matrix as it should.
            logger.error('it appears that single value decomposition has failed. This occurs for a few reasons:'
                         '\n-  Very small or very large data in the input coordinates'
                         '\n-  The data does not cover a full sphere: hence the matrix is poorly conditioned and has no inverse'
                         '\n-  Continuing but you should be careful...')

        harmonics = inverseLegendre @ self._input_Bz_data.Bz
        self.harmonics = pd.Series(harmonics, index=self._coeff_names)

    def _quantify_fit(self):
        """
        Compare the reconstructed field to the original field and calculate some basic goodness of fit metrics
        """
        if len(self._input_Bz_data.index) > 1e4:
            logger.warning('you are reconstructing a lot of points and it might be a bit slow.'
                           'I should write a down sampling routine here...')

        Bz_recon = self.legendre_basis @ (self.harmonics / self.scale)
        Residual = np.subtract(self._input_Bz_data.Bz, Bz_recon)
        self._residual_pk_pk = float(abs(Residual.max() - Residual.min()) * 1e6)

        initial_pk_pk = float(abs(self._input_Bz_data.Bz.max() - self._input_Bz_data.Bz.min()) * 1e6)
        recon_pk_pk = float(abs(Bz_recon.max() - Bz_recon.min()) * 1e6)
        residual_percentage = abs(self._residual_pk_pk) * 100 / initial_pk_pk

        try:
            print(f'Initial pk-pk:       {initial_pk_pk: 1.3e} \u03BCT')
            print(f'Reconstructed pk-pk: {recon_pk_pk: 1.3e} \u03BCT')
            print(f'Residual pk-pk:      {self._residual_pk_pk: 1.3e} \u03BCT ({residual_percentage: 1.1f}%)')
        except UnicodeError:
            print(f'Initial pk-pk:       {initial_pk_pk: 1.3e} uT')
            print(f'Reconstructed pk-pk: {recon_pk_pk: 1.3e} uT')
            print(f'Residual pk-pk:      {self._residual_pk_pk: 1.3e} uT ({residual_percentage: 1.1f}%)')

        if residual_percentage > 2:
            logger.warning('residual_pk_pk is greater than 2 %. This may indicate that the order is not high enough,'
                           'or that the data is non physical (but it heavily depends on your use case!')

    def _assess_harmonic_pk_pk(self):
        """
        generate the peak-peak perturbation caused over the surface of the sphere due to each harmonic
        """
        if len(self.legendre_basis.index) < 1e4:
            # generate some points to ensure we cover a good portion of the sphere
            r = self.r_outer
            azimuth = np.linspace(0, 2 * np.pi, 100)
            elevation = np.linspace(0, np.pi, 100)
            [r, azimuth, elevation] = np.meshgrid(r, azimuth, elevation)
            coords = pd.DataFrame({'r': r.flatten(), 'azimuth': azimuth.flatten(), 'elevation': elevation.flatten()})
            legendre_basis_to_use = generate_legendre_basis(coords, self.n_order)
        else:
            legendre_basis_to_use = self.legendre_basis

        BasisRange = legendre_basis_to_use.apply(lambda x: abs(np.max(x) - np.min(x)))
        self.HarmonicsPk_Pk = (self.harmonics / self.scale) * BasisRange * 1e6

    # Public Methods

    def plot_harmonics_pk_pk(self, cut_off=.1, title=None, return_axs=False, plot_percentage_of_dominant=False,
                             drop_dominant_harmonic=False):  # pragma: no cover
        """
                produces a barplot of harmonics.

        :param cut_off: cutoff point relative to highest harmonic. e.g. cut_off=.1 means that harmonics which produce
            less than 10% the pk-pk perturbation of the dominant harmonic are ignored.
        :type cut_off: float, optional
        :param title: title of plot
        :param return_axs: if True, will return the plotting axs for further manipulation
        :param plot_percentage_of_dominant: if True, switches from absolute pk-pk to percentage
        :param drop_dominant_harmonic: if True, drops dominant harmonic before plotting
        """

        if not hasattr(self, 'HarmonicsPk_Pk'):
            self._assess_harmonic_pk_pk()
        plt.figure(figsize=[8, 8])
        # ax = sns.barplot(self.HarmonicsPk_Pk, y='pk-pk [\u03BCT]', x=)
        CutOffInd = abs(self.HarmonicsPk_Pk) < cut_off * abs(self.HarmonicsPk_Pk).max()
        KeyHarmonics = self.HarmonicsPk_Pk.drop(self.HarmonicsPk_Pk[CutOffInd].index)
        HarmonicsToPlot = KeyHarmonics

        _ylabel = 'pk-pk [\u03BCT]'
        if plot_percentage_of_dominant:
            _ylabel = 'pk-pk [%]'
            HarmonicsToPlot = HarmonicsToPlot * 100 / HarmonicsToPlot.abs().max()
        dominant_ind = np.argwhere(np.max(HarmonicsToPlot))
        if drop_dominant_harmonic:
            HarmonicsToPlot = HarmonicsToPlot.drop(HarmonicsToPlot.index[HarmonicsToPlot.abs().argmax()])
        sns.set_theme(font_scale=2)
        axs = sns.barplot(x=HarmonicsToPlot.index, y=HarmonicsToPlot.values, palette="Blues_d")
        
        if title is None:
            axs.set_title(f'Principle Harmonics pk-pk (>{cut_off * 100: 1.0f}% of max)')
        else:
            axs.set_title(title)

        axs.set_ylabel(_ylabel)
        for item in axs.get_xticklabels():
            item.set_rotation(45)
        plt.tight_layout()
        if not return_axs:
            plt.show()
        else:
            return axs

    def plot_cut_planes(self, resolution=2.5, AddColorBar=True, quantity='uT', vmin=None,
                        vmax=None):  # pragma: no cover
        """
        Reconstruct the Bz field at the cardinal planes.
        Note this is basically the same code copied three times in  a row, one for each plane.

        :param resolution: the resolution of reconstructed data. Default is 2.5 mm
        :type resolution: float, optional
        :param AddColorBar: Add a color bar to each plot
        :type AddColorBar: Boolean, optional
        :param quantity: quantity to use in reconstruct_Bz; must match the options for that function
        :type quantity: string, optional
        :param vmin: same as matplotlib.pyplot.imshow; vmin and vmax control the color map bounds.
        :type vmin: float, optional
        :param vmax: same as matplotlib.pyplot.imshow; vmin and vmax control the color map bounds.
        :type vmax: float, optional
        """

        fig, (axs1, axs2, axs3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
        extent = [-self.r_outer - resolution, self.r_outer + resolution, -self.r_outer - resolution,
                  self.r_outer + resolution]  # same for all plots

        conversion_factor = 1
        # XY plane
        x = np.linspace(-self.r_outer, self.r_outer, int((2 * self.r_outer) / resolution)) * conversion_factor
        y = np.linspace(-self.r_outer, self.r_outer, int((2 * self.r_outer) / resolution)) * conversion_factor
        z = 0
        [x_recon, y_recon, z_recon] = np.meshgrid(x, y, z, indexing='ij')
        recon_coords = pd.DataFrame({'x': x_recon.flatten(), 'y': y_recon.flatten(), 'z': z_recon.flatten()})
        recon_coords = convert_cartesian_to_spherical(recon_coords)
        Bz_recon = reconstruct_Bz(harmonics=self.harmonics / self.scale, coords=recon_coords, quantity=quantity)
        Bz_recon = Bz_recon.to_numpy().reshape(np.squeeze(x_recon).shape)
        XYimage = axs1.imshow(Bz_recon.T, extent=extent, vmin=vmin, vmax=vmax)
        axs1.set_xlabel('x [mm]')
        axs1.set_ylabel('y [mm]')
        axs1.set_title('a) XY plane')
        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs1.grid(False)
        axs1.add_artist(draw_circle)

        # ZX plane
        x = np.linspace(-self.r_outer, self.r_outer, int((2 * self.r_outer) / resolution)) * conversion_factor
        y = 0
        z = np.linspace(-self.r_outer, self.r_outer, int((2 * self.r_outer) / resolution)) * conversion_factor
        [x_recon, y_recon, z_recon] = np.meshgrid(x, y, z, indexing='ij')
        recon_coords = pd.DataFrame({'x': x_recon.flatten(), 'y': y_recon.flatten(), 'z': z_recon.flatten()})
        recon_coords = convert_cartesian_to_spherical(recon_coords)
        Bz_recon = reconstruct_Bz(harmonics=self.harmonics / self.scale, coords=recon_coords, quantity=quantity)
        Bz_recon = Bz_recon.to_numpy().reshape(np.squeeze(y_recon).shape)
        ZXimage = axs2.imshow(Bz_recon.T, extent=extent, vmin=vmin, vmax=vmax)
        axs2.set_xlabel('x [mm]')
        axs2.set_ylabel('z [mm]')
        axs2.set_title('b) ZX plane')
        axs2.grid(False)
        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs2.add_artist(draw_circle)

        # ZY plane
        x = 0
        y = np.linspace(-self.r_outer, self.r_outer, int((2 * self.r_outer) / resolution)) * conversion_factor
        z = np.linspace(-self.r_outer, self.r_outer, int((2 * self.r_outer) / resolution)) * conversion_factor
        [x_recon, y_recon, z_recon] = np.meshgrid(x, y, z, indexing='ij')
        recon_coords = pd.DataFrame({'x': x_recon.flatten(), 'y': y_recon.flatten(), 'z': z_recon.flatten()})
        recon_coords = convert_cartesian_to_spherical(recon_coords)
        Bz_recon = reconstruct_Bz(harmonics=self.harmonics / self.scale, coords=recon_coords, quantity=quantity)
        Bz_recon = Bz_recon.to_numpy().reshape(np.squeeze(z_recon).shape)
        ZYimage = axs3.imshow(Bz_recon.T, extent=extent, vmin=vmin, vmax=vmax)
        axs3.set_xlabel('y [mm]')
        axs3.set_ylabel('z [mm]')
        axs3.grid(False)
        axs3.set_title('c) ZY plane')
        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs3.add_artist(draw_circle)

        if AddColorBar:
            if quantity == 'uT':
                cb_title = '\u03BCT'
            elif quantity == 'T':
                cb_title = 'T'
            cbar = fig.colorbar(XYimage, ax=axs1, fraction=0.046, pad=0.04)
            cbar.set_label(cb_title)
            cbar = fig.colorbar(ZXimage, ax=axs2, fraction=0.046, pad=0.04)
            cbar.set_label(cb_title)
            cbar = fig.colorbar(ZYimage, ax=axs3, fraction=0.046, pad=0.04)
            cbar.set_label(cb_title)

        plt.tight_layout()
        plt.show()

    def print_key_harmonics(self, cut_off=.01):
        """
        print the harmonics with value > cut_off to the terminal in pk-pk.

        :param cut_off: cutoff point relative to highest harmonic. e.g. cut_off=.1 means that harmonics which produce
            less than 10% the pk-pk perturbation of the dominant harmonic are ignored.
        :type cut_off: float, optional
        """
        if not hasattr(self, 'harmonics_pk_pk'):
            self._assess_harmonic_pk_pk()

        CutOffInd = abs(self.HarmonicsPk_Pk) < cut_off * abs(self.HarmonicsPk_Pk).max()
        KeyHarmonics = self.HarmonicsPk_Pk.drop(self.HarmonicsPk_Pk[CutOffInd].index)
        pd.set_option('display.float_format', lambda x: '%1.4ef' % x)
        print(f'\nOnly displaying values >= {cut_off * 100: 1.0f}% of the peak harmonic.')
        try:
            print('Values are in pk=pk [\u03BCT]')
        except UnicodeError:
            print('Values are in pk=pk uT]')
        print(f'{bcolors.OKBLUE}{KeyHarmonics.to_string()}{bcolors.ENDC}')
