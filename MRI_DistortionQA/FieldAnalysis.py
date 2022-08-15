from pathlib import Path
import sys, os
import numpy as np
import pandas as pd
import logging
from matplotlib import pyplot as plt
from .utilities import bcolors
from .utilities import convert_cartesian_to_spherical, generate_legendre_basis
from .utilities import reconstruct_Bz
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

    :param InputData: A pandas dataframe with columns x, y, z, Bz. Data should be in mm and T
    :type InputData: Pandas.DataFrame
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
        field harmonics to gradient strength. To work with our distortion correction code, you should scale each
        gradient by 1/gradient_strength in mT/m; e.g. if the gradient strength is 10 mT/m then scale should be 1/10
    :type scale: float, optional
    """
    def __init__(self, input_Bz_data, r_outer=150, n_order=8, AssessHarmonicPk_Pk=True, QuantifyFit=True,
                 TrimDataBy_r_outer=False, scale=1):

        # attributes:
        self._tol = 1e-3  # tolerance for coordinate filtering in mm
        self.scale = scale
        self.TrimDataBy_r_outer = TrimDataBy_r_outer
        self.r_outer = r_outer
        self.QuantifyFit = QuantifyFit
        self.AssessHarmonicPk_Pk = AssessHarmonicPk_Pk
        self.n_order = n_order
        self.n_harmonics = (self.n_order + 1) ** 2
        self.input_Bz_data = input_Bz_data.copy()
        self.__recon_warning_thrown = False
        # basic fitting routine:
        self.input_Bz_data = convert_cartesian_to_spherical(self.input_Bz_data)
        self._check_data_input()
        if self.TrimDataBy_r_outer:
            self._filter_data()
        self._generate_harmonic_names()
        self.legendre_basis = generate_legendre_basis(self.input_Bz_data, self.n_order)
        self._svd_fit()
        # optional methods:
        if self.QuantifyFit:
            self._quantify_fit()

        if self.AssessHarmonicPk_Pk:
            self._assess_harmonic_pk_pk()
        
        if not self.scale == 1:
            logger.warning(f'scaling harmonics by {self.scale}')
            self.harmonics = self.harmonics*scale
        
    def _check_data_input(self):
        """
        - Make sure the input data actually covers a sphere. If it doesn't, it normally means something has or will go wrong
            so a warning is triggered.
        - Make sure that input data is in mm, which is an implicit assumption
        """

        if self.input_Bz_data.r.mean() < 10:
            logger.warning('it appears that your input data is in m, not mm - please use mm!')

        if (self.input_Bz_data.elevation.max() - self.input_Bz_data.elevation.min()) < 0.75 * np.pi or \
                (self.input_Bz_data.azimuth.max() - self.input_Bz_data.azimuth.min()) < 1.8 * np.pi:
            logger.warning('input sample points do not appear to cover a full sphere')

    def _filter_data(self):
        """
        filter the Bz data by radial coordinate, removing any entries that fall outside self.r_outer
        """

        data_to_delete_ind = self.input_Bz_data.r > (self.r_outer + self._tol)
        self.input_Bz_data = self.input_Bz_data.drop(self.input_Bz_data[data_to_delete_ind].index)
        self.input_Bz_data.reset_index(inplace=True)
        n_deleted = np.count_nonzero(data_to_delete_ind)
        if n_deleted > 0:
            logger.warning(f'deleting {n_deleted} of the original {data_to_delete_ind.size} data points because'
                           f'they fall outside the r_outer value of {self.r_outer} and TrimDataBy_r_outer=True')

        if self.input_Bz_data.shape[0] == 0:
            logger.error(f'After filtering for data insdide {self.r_outer} mm, there is no data left!')
            sys.exit(1)

    def _generate_harmonic_names(self):
        """
        generate the names of each harmonic. Used to label the columns in the harmonics dataframe
        """
        k = 0
        self._coeff_names = []
        for n in range(0, self.n_order + 1):  # the plus 1 is because range stops at -1 for some reason
            self._coeff_names.append(f'A_{n}_0')
            k = k + 1
            for m in range(0, n):
                self._coeff_names.append(f'A_{n}_{m+1}')
                k = k + 1
                self._coeff_names.append(f'B_{n}_{m + 1}')
                k = k + 1

    def _svd_fit(self):
        """
        calculate harmonics using single value decomposition of legendre_basis

        Formalism:

        Let the legendre basis be called L, the Harmonics be H and the B field B.
        L will have size [n_coords, n_harmonics], Harmonics has size [n_harmonics, 1] and B has size [n_coords, 1]

        Our proposition is that there exists some set of harmonics such that

        .. math:: L.H = B

        We can find H by inverting L:

        .. math::

            L^{-1}.L.H = L^{-1}.B \\
            H = L^{-1}.B

        However, direct matrix inversion is tricky. We therefore use single value decomposition to calculate the pseudo inverse
        """

        inverseLegendre = np.linalg.pinv(self.legendre_basis)

        if not np.allclose(inverseLegendre @ self.legendre_basis, np.eye(inverseLegendre.shape[0], self.legendre_basis.shape[1]),
                           rtol=5e-3, atol=5e-4):
            # check that the inverse times itself returns an idendity matrix as it should.
            logger.error('it appears that single value decomposition has failed. This occurs for a few reasons:'
                         '\n-  Very small or very large data in the input coordinates'
                         '\n-  The data does not cover a full sphere: hence the matrix is poorly conditioned and has no inverse'
                         '\n-  Continuing but you should be careful...')

        harmonics = inverseLegendre @ self.input_Bz_data.Bz
        self.harmonics = pd.Series(harmonics, index=self._coeff_names)

    def _quantify_fit(self):
        """
        Compare the reconstructed field to the original field and calculate some basic goodness of fit metrics
        """
        if len(self.input_Bz_data.index) > 1e4:
            logger.warning('you are reconstructing a lot of points and it might be a bit slow.'
                           'I should write a down sampling routine here...')

        Bz_recon = self.legendre_basis @ self.harmonics
        Residual = np.subtract(self.input_Bz_data.Bz, Bz_recon)

        initial_pk_pk = float(abs(Bz_recon.max() - Bz_recon.min()) * 1e6)
        recon_pk_pk = float(abs(self.input_Bz_data.Bz.max() - self.input_Bz_data.Bz.min()) * 1e6)
        self._residual_pk_pk = float(abs(Residual.max() - Residual.min()) * 1e6)

        print(f'Initial pk-pk:       {initial_pk_pk: 1.3f} \u03BCT')
        print(f'Reconstructed pk-pk: {recon_pk_pk: 1.3f} \u03BCT')
        print(f'Residual pk-pk:      {self._residual_pk_pk: 1.3f} \u03BCT')

        residual_percentage = abs(self._residual_pk_pk)*100/initial_pk_pk
        if residual_percentage > 5:
            logger.warning('residual_pk_pk is greater than 5%. This may indicate that the order is not high enough,'
                           'or that the data is non physical (but it heavily depends on your use case!')

    def _assess_harmonic_pk_pk(self):
        """
        generate the peak-peak perturbation caused over the surface of the sphere due to each harmonic
        """
        if len(self.legendre_basis.index) < 1e4:
            # generate some points to ensure we cover a good portion of the sphere
            r = self.r_outer
            azimuth = np.linspace(0, 2*np.pi, 100)
            elevation = np.linspace(0, np.pi, 100)
            [r, azimuth, elevation] = np.meshgrid(r, azimuth, elevation)
            coords = pd.DataFrame({'r': r.flatten(), 'azimuth': azimuth.flatten(), 'elevation': elevation.flatten()})
            legendre_basis_to_use = generate_legendre_basis(coords, self.n_order)
        else:
            legendre_basis_to_use = self.legendre_basis

        BasisRange = legendre_basis_to_use.apply(lambda x: abs(np.max(x) - np.min(x)))
        self.HarmonicsPk_Pk = self.harmonics * BasisRange * 1e6

    # Public Methods

    def plot_harmonics_pk_pk(self, cut_off=.1):  # pragma: no cover
        """
        produces a barplot of harmonics.

        :param cut_off: cutoff point relative to highest harmonic. e.g. cut_off=.1 means that harmonics which produce
            less than 10% the pk-pk perturbation of the dominant harmonic are ignored.
        :type cut_off: float, optional
        """

        if not hasattr(self,'HarmonicsPk_Pk'):
            self._assess_harmonic_pk_pk()
        plt.figure(figsize=[10, 5])
        # ax = sns.barplot(self.HarmonicsPk_Pk, y='pk-pk [\u03BCT]', x=)
        CutOffInd = abs(self.HarmonicsPk_Pk) < cut_off * abs(self.HarmonicsPk_Pk).max()
        KeyHarmonics = self.HarmonicsPk_Pk.drop(self.HarmonicsPk_Pk[CutOffInd].index)
        HarmonicsToPlot = KeyHarmonics

        axs = sns.barplot(x=HarmonicsToPlot.index, y=HarmonicsToPlot.values,  palette="Blues_d")
        axs.set_title(f'Principle Harmonics pk-pk (>{cut_off*100: 1.0f}% of max)')

        axs.set_ylabel('pk-pk [\u03BCT]')
        for item in axs.get_xticklabels():
            item.set_rotation(45)
        plt.show()

    def plot_cut_planes(self, resolution=2.5, AddColorBar=True, quantity='uT', vmin=None, vmax=None):  # pragma: no cover
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
        x = np.linspace(-self.r_outer, self.r_outer, int((2*self.r_outer)/resolution)) * conversion_factor
        y = np.linspace(-self.r_outer, self.r_outer, int((2*self.r_outer)/resolution)) * conversion_factor
        z = 0
        extent = [-self.r_outer - resolution, self.r_outer + resolution, -self.r_outer - resolution,
                  self.r_outer + resolution]  # same for all plots
        [x_recon, y_recon, z_recon] = np.meshgrid(x, y, z, indexing='ij')
        recon_coords = pd.DataFrame({'x': x_recon.flatten(), 'y': y_recon.flatten(), 'z': z_recon.flatten()})
        recon_coords = convert_cartesian_to_spherical(recon_coords)
        Bz_recon = reconstruct_Bz(harmonics=self.harmonics, coords=recon_coords, quantity=quantity)
        Bz_recon = Bz_recon.to_numpy().reshape(np.squeeze(x_recon).shape)
        XYimage = axs1.imshow(Bz_recon.T, extent=extent, vmin=vmin, vmax=vmax)
        axs1.set_xlabel('x [mm]')
        axs1.set_ylabel('y [mm]')
        axs1.set_title('a) XY plane')
        draw_circle = plt.Circle((0, 0), 150, fill=False)
        axs1.grid(False)
        axs1.add_artist(draw_circle)

        # ZX plane
        x = np.linspace(-self.r_outer,self.r_outer,int((2*self.r_outer)/resolution)) * conversion_factor
        y = 0
        z = np.linspace(-self.r_outer, self.r_outer, int((2*self.r_outer)/resolution)) * conversion_factor
        [x_recon, y_recon, z_recon] = np.meshgrid(x, y, z, indexing='ij')
        recon_coords = pd.DataFrame({'x': x_recon.flatten(), 'y': y_recon.flatten(), 'z': z_recon.flatten()})
        recon_coords = convert_cartesian_to_spherical(recon_coords)
        Bz_recon = reconstruct_Bz(harmonics=self.harmonics, coords=recon_coords, quantity=quantity)
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
        y = np.linspace(-self.r_outer,self.r_outer,int((2*self.r_outer)/resolution)) * conversion_factor
        z = np.linspace(-self.r_outer,self.r_outer,int((2*self.r_outer)/resolution)) * conversion_factor
        [x_recon, y_recon, z_recon] = np.meshgrid(x, y, z, indexing='ij')
        recon_coords = pd.DataFrame({'x': x_recon.flatten(), 'y': y_recon.flatten(), 'z': z_recon.flatten()})
        recon_coords = convert_cartesian_to_spherical(recon_coords)
        Bz_recon = reconstruct_Bz(harmonics=self.harmonics, coords=recon_coords, quantity=quantity)
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

    def print_key_harmonics(self, cut_off=.1):
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

        print(f'\nOnly displaying values >= {cut_off*100: 1.0f}% of the peak harmonic.')
        print('Values are in pk=pk [\u03BCT]')
        print(f'{bcolors.OKBLUE}{KeyHarmonics.to_string()}{bcolors.ENDC}')


