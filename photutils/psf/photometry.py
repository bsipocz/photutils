# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Module which provides classes to perform PSF Photometry"""

from __future__ import division
import abc
import numpy as np
from astropy.extern import six
from astropy.table import Table, vstack
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.nddata.utils import overlap_slices, NoOverlapError
from astropy.stats import gaussian_sigma_to_fwhm
from ..psf import subtract_psf
from ..aperture import CircularAperture, aperture_photometry


__all__ = ['DAOPhotPSFPhotometry']


class DAOPhotPSFPhotometry(object):
    """
    This class implements the DAOPHOT algorithm proposed by Stetson
    (1987) to perform point spread function photometry in crowded fields,
    which consists basically in applying the loop FIND, GROUP, NSTAR,
    SUBTRACT, FIND until no more stars are detected or a given number of
    iterations is reached.
    """

    def __init__(self, group, bkg, psf, fitshape, find=None,
                 fitter=LevMarLSQFitter(), niters=3, aperture_radius=None):
        """
        Attributes
        ----------
        find : callable or instance of any StarFinderBase subclasses
            ``find`` should be able to identify stars, i.e. compute a rough
            estimate of the centroids, in a given 2D image.
            ``find`` receives as input a 2D image an return an
            `~astropy.table.Table` object which contains columns with names:
            ``id``, ``xcentroid``, ``ycentroid``, and ``flux``. In which
            ``id`` is an interger-valued column starting from ``1``,
            ``xcentroid`` and ``ycentroid`` are center position estimates of
            the sources and ``flux`` contains flux estimates of the sources.
            See, e.g., `~photutils.detection.DAOStarFinder`
        group : callable or instance of any GroupStarsBase subclasses
            ``group`` should be able to decide whether a given star overlaps
            with any other and label them as beloging to the same group.
            ``group`` receives as input an `~astropy.table.Table` object
            with columns named as ``id``, ``x_0``, ``y_0``, in which ``x_0``
            and ``y_0`` have the same meaning of ``xcentroid`` and
            ``ycentroid``. This callable must return an `~astropy.table.Table`
            with columns ``id``, ``x_0``, ``y_0``, and ``group_id``. The
            column ``group_id`` should cotain integers starting from ``1``
            that indicate in which group a given source belongs to.
            See, e.g., `~photutils.psf.DAOGroup`
        bkg : callable or instance of any `~photutils.BackgroundBase`
            subclasses ``bkg`` should be able to compute either a scalar
            background or a 2D background of a given 2D image.
            See, e.g., `~photutils.background.MedianBackground`
        psf : `astropy.modeling.Fittable2DModel` instance
            PSF or PRF model to fit the data. Could be one of the models in
            this package like `~photutils.psf.sandbox.DiscretePRF`,
            `~photutils.psf.IntegratedGaussianPRF`, or any other suitable
            2D model.
            This object needs to identify three parameters (position of
            center in x and y coordinates and the flux) in order to set them
            to suitable starting values for each fit. The names of these
            parameters should be given as ``x_0``, ``y_0`` and ``flux``.

            `~photutils.psf.prepare_psf_model` can be used to prepare any 2D
            model to match this assumption.
        fitshape : array-like
            Rectangular shape around the center of a star which will be used
            to collect the data to do the fitting, e.g. (5, 5) means to take
            the following relative pixel positions: [-2, -1, 0, 1, 2].
            Also, each element of ``fitshape`` must be an odd number.
        fitter : Fitter instance (default=LevMarLSQFitter())
            Fitter object used to compute the optimized centroid positions
            and/or flux of the identified sources. See
            `~astropy.modeling.fitting` for more details on fitters.
        niters : int (default=3)
            Number of iterations to perform the loop FIND, GROUP, SUBTRACT,
            NSTAR.
        aperture_radius : float (default=None)
            The radius (in units of pixels) used to compute initial estimates
            for the fluxes of sources. If ``None``, one fwhm will be used. 

        Notes
        -----
        If there are problems with fitting large groups, change the parameters
        of the grouping algorithm to reduce the number of sources in each
        group or input a ``star_groups`` table that only includes the groups
        that are relevant (e.g. manually remove all entries that coincide with
        artifacts).

        References
        ----------
        [1] Stetson, Astronomical Society of the Pacific, Publications,
            (ISSN 0004-6280), vol. 99, March 1987, p. 191-222.
            Available at: http://adsabs.harvard.edu/abs/1987PASP...99..191S
        """

        self.find = find
        self.group = group
        self.bkg = bkg
        self.psf = psf
        self.fitter = fitter
        self.niters = niters
        self.fitshape = fitshape
        self.aperture_radius = aperture_radius

    @property
    def niters(self):
        return self._niters

    @niters.setter
    def niters(self, value):
        try:
            if value <= 0:
                raise ValueError('niters must be positive.')
            else:
                self._niters = int(value)
        except:
            raise ValueError('niters must be an integer or convertable '
                             'into an integer.')
    
    @property
    def fitshape(self):
        return self._fitshape

    @fitshape.setter
    def fitshape(self, value):
        value = np.asarray(value)
        if value.size == 2:
            if np.all(value) > 0:
                if np.all(value % 2) == 1:
                    self._fitshape = tuple(value)
                else:
                    raise ValueError('fitshape must be odd integer-valued, '
                                     'received fitshape = {}'\
                                     .format(value))
            else:
                raise ValueError('fitshape must have positive elements, '
                                 'received fitshape = {}'\
                                 .format(value))
        else:
            raise ValueError('fitshape must have two dimensions, '
                             'received fitshape = {}'.format(value))

    @property
    def aperture_radius(self):
        return self._aperture_radius

    @aperture_radius.setter
    def aperture_radius(self, value):
        if isinstance(value, (int, float)) and value > 0:
            self._aperture_radius = value
        elif value is None:
            self._aperture_radius = value
        else:
            raise ValueError('aperture_radius must be a real-valued '
                             'number, received aperture_radius = {}'
                             .format(value))


    def __call__(self, image, positions=None):
        """
        Parameters
        ----------
        image : 2D array-like, `~astropy.io.fits.ImageHDU`,
        `~astropy.io.fits.HDUList`
            Image to perform photometry.
        positions : `~astropy.table.Table` (optional)
            Positions, in pixel coordinates, at which stars are located.
            The columns must be named as 'x_0' and 'y_0'. 'flux_0' can also
            be provided to set initial fluxes.

        Returns
        -------
        outtab : `~astropy.table.Table`
            Table with the photometry results, i.e., centroids and fluxes
            estimations.
        residual_image : array-like, `~astropy.io.fits.ImageHDU`,
        `~astropy.io.fits.HDUList`
            Residual image calculated by subtracting the fitted sources
            and the original image.
        """
        
        if positions is None:
            return self.do_photometry(image)
        else:
            return self.do_fixed_photometry(image, positions)

    def do_photometry(self, image):
        """
        Perform PSF photometry in ``image``. This method assumes that
        ``psf`` has centroids and flux parameters which will be fitted to the
        data provided in ``image``. A compound model, in fact a sum of
        ``psf``, will be fitted to groups of starts automatically identified
        by ``group``. Also, ``image`` is not assumed to be background
        subtracted.

        Parameters
        ----------
        image : 2D array-like, `~astropy.io.fits.ImageHDU`,
        `~astropy.io.fits.HDUList`
            Image to perform photometry.
        
        Returns
        -------
        outtab : `~astropy.table.Table`
            Table with the photometry results, i.e., centroids and fluxes
            estimations.
        residual_image : array-like, `~astropy.io.fits.ImageHDU`,
        `~astropy.io.fits.HDUList`
            Residual image calculated by subtracting the fitted sources
            and the original image.
        """

        outtab = Table([[], [], [], [], [], []],
                       names=('id', 'group_id', 'x_fit', 'y_fit', 'flux_fit',
                              'iter_detected'),
                       dtype=('i4', 'i4', 'f8', 'f8', 'f8', 'i4'))

        residual_image = image.copy()
        residual_image = residual_image - self.bkg(image)
        sources = self.find(residual_image)
        
        if self.aperture_radius is None:
            if hasattr(self.psf, 'fwhm'):
                self.aperture_radius = self.psf.fwhm.value
            elif hasattr(self.psf, 'sigma'):
                self.aperture_radius = self.psf.sigma.value*gaussian_sigma_to_fwhm

        apertures = CircularAperture((sources['xcentroid'],
                                      sources['ycentroid']),
                                     r=self.aperture_radius)

        sources['aperture_flux'] = aperture_photometry(residual_image,
                                              apertures)['aperture_sum']
        n = 1
        while(n <= self.niters and len(sources) > 0):
            intab = Table(names=['id', 'x_0', 'y_0', 'flux_0'],
                          data=[sources['id'], sources['xcentroid'],
                          sources['ycentroid'], sources['aperture_flux']])
            star_groups = self.group(intab)
            tab, residual_image = self.nstar(residual_image, star_groups)
            tab['iter_detected'] = n*np.ones(tab['x_fit'].shape, dtype=np.int)
            outtab = vstack([outtab, tab])
            sources = self.find(residual_image)

            if len(sources) > 0:
                apertures = CircularAperture((sources['xcentroid'],
                                              sources['ycentroid']),
                                             r=self.aperture_radius)
                sources['flux'] = aperture_photometry(residual_image,
                                                      apertures)['aperture_sum']
            n += 1
        return outtab, residual_image

    def do_fixed_photometry(self, image, positions):
        """
        Perform PSF photometry for the case that the centroid positions of the
        starts are known with high accuracy. If the centroid positions are set
        as ``fixed`` in the PSF model ``psf``, then the optimizer will only
        consider the flux as a variable. Otherwise, ``positions`` will be used
        as initial guesses for the centroids. Also, ``image`` is not assumed
        to be background subtracted.

        Parameters
        ----------
        image : 2D array-like, `~astropy.io.fits.ImageHDU`,
        `~astropy.io.fits.HDUList`
            Image to perform photometry
        positions: `~astropy.table.Table`
            Positions (in pixel coordinates) at which to *start* the fit for
            each object. Columns 'x_0' and 'y_0' must be present.
            'flux_0' can also be provided to set initial fluxes.

        Returns
        -------
        outtab : `~astropy.table.Table`
            Table with the photometry results, i.e., centroids and fluxes
            estimations.
        residual_image : array-like, `~astropy.io.fits.ImageHDU`,
        `~astropy.io.fits.HDUList`
            Residual image calculated by subtracting the fitted sources
            from the original image.
        """

        residual_image = image.copy()
        residual_image = residual_image - self.bkg(image)

        if 'flux_0' not in positions.colnames:
            if self.aperture_radius is None:
                 self.aperture_radius = self.psf.sigma.value*\
                                        gaussian_sigma_to_fwhm
            apertures = CircularAperture((positions['x_0'], positions['y_0']),
                                         r=self.aperture_radius)

            positions['aperture_flux'] = aperture_photometry(residual_image,\
                                         apertures)['aperture_sum']

        intab = Table(names=['x_0', 'y_0', 'flux_0'],
                      data=[positions['x_0'], positions['y_0'],
                      positions['aperture_flux']])

        star_groups = self.group(intab)
        outtab, residual_image = self.nstar(residual_image, star_groups)
        
        return outtab, residual_image

    def get_uncertainties(self):
        """
        Return the uncertainties on the fitted parameters
        """
        raise NotImplementedError

    def nstar(self, image, star_groups):
        """
        Fit, as appropriate, a compound or single model to the given
        ``star_groups``. Groups are fitted sequentially from the smallest to
        the biggest. In each iteration, ``image`` is subtracted by the
        previous fitted group.
        
        Parameters
        ----------
        image : numpy.ndarray
            Background-subtracted image.
        star_groups : `~astropy.table.Table`
            This table must contain the following columns: ``id``,
            ``group_id``, ``x_0``, ``y_0``, ``flux_0``.
            ``x_0`` and ``y_0`` are initial estimates of the centroids
            and ``flux_0`` is an initial estimate of the flux.

        Returns
        -------
        result_tab : `~astropy.table.Table`
            Astropy table that contains photometry results.
        image : numpy.ndarray
            Residual image.
        """

        result_tab = Table([[], [], [], [], []],
                           names=('id', 'group_id', 'x_fit', 'y_fit',
                                  'flux_fit'),
                           dtype=('i4', 'i4', 'f8', 'f8', 'f8'))
        star_groups = star_groups.group_by('group_id')
        
        y, x = np.indices(image.shape)

        for n in range(len(star_groups.groups)):
            group_psf = self.GroupPSF(self.psf,
                                      star_groups.groups[n]).get_model()
            usepixel = np.zeros_like(image, dtype=np.bool)
            
            for row in star_groups.groups[n]:
                usepixel[overlap_slices(large_array_shape=image.shape,
                                        small_array_shape=self.fitshape,
                                        position=(row['y_0'], row['x_0']),
                                        mode='trim')[0]] = True

            fit_model = self.fitter(group_psf, x[usepixel], y[usepixel],
                                    image[usepixel])
            param_table = self._model_params2table(fit_model,
                                                   star_groups.groups[n])
            result_tab = vstack([result_tab, param_table])
           
            # do not subtract if the fitting did not go well
            try:
                image = subtract_psf(image, self.psf, param_table,
                                     subshape=self.fitshape)
            except NoOverlapError:
                pass

        return result_tab, image

    def _model_params2table(self, fit_model, star_group):
        """
        Place fitted parameters into an astropy table.
        
        Parameters
        ----------
        fit_model : `astropy.modeling.Fittable2DModel` instance
            PSF or PRF model to fit the data. Could be one of the models in
            this package like `~photutils.psf.sandbox.DiscretePRF`,
            `~photutils.psf.IntegratedGaussianPRF`, or any other suitable
            2D model.
        star_group : ~astropy.table.Table
        
        Returns
        -------
        param_tab : ~astropy.table.Table
            Table that contains the fitted parameters.
        """

        param_tab = Table([[], [], [], [], []],
                          names=('id', 'group_id', 'x_fit', 'y_fit',
                                 'flux_fit'),
                          dtype=('i4', 'i4', 'f8', 'f8', 'f8'))
        if np.size(fit_model) == 1:
            param_tab.add_row([[star_group['id'][0]],
                               [star_group['group_id'][0]],
                               [getattr(fit_model,'x_0').value],
                               [getattr(fit_model, 'y_0').value],
                               [getattr(fit_model, 'flux').value]])
        else:
            for i in range(np.size(fit_model)):
                param_tab.add_row([[star_group['id'][i]],
                                   [star_group['group_id'][i]],
                                   [getattr(fit_model,'x_0_'+str(i)).value],
                                   [getattr(fit_model, 'y_0_'+str(i)).value],
                                   [getattr(fit_model, 'flux_'+str(i)).value]])
        return param_tab

    class GroupPSF(object):
        """
        Construct a joint psf model which consists in a sum of `self.psf`
        whose parameters are given in `star_group`.

        Attributes
        ----------
        star_group : `~astropy.table.Table`
            Table from which the compound PSF will be constructed.
            It must have columns named as `x_0`, `y_0`, and `flux_0`.
        psf : `astropy.modeling.Fittable2DModel` instance
        """

        def __init__(self, psf, star_group):
            self.star_group = star_group
            self.psf = psf
        
        def get_model(self):
            """        
            Returns
            -------
            group_psf : CompoundModel
                `CompoundModel` instance which is a sum of the given PSF
                models.
            """
            
            psf_class = type(self.psf)
            group_psf = psf_class(sigma=self.psf.sigma.value,
                                  flux=self.star_group['flux_0'][0],
                                  x_0=self.star_group['x_0'][0],
                                  y_0=self.star_group['y_0'][0],
                                  fixed=self.psf.fixed, tied=self.psf.tied,
                                  bounds=self.psf.bounds)
            for i in range(len(self.star_group) - 1):
                group_psf += psf_class(sigma=self.psf.sigma.value,
                                       flux=self.star_group['flux_0'][i+1],
                                       x_0=self.star_group['x_0'][i+1],
                                       y_0=self.star_group['y_0'][i+1],
                                       fixed=self.psf.fixed,
                                       tied=self.psf.tied,
                                       bounds=self.psf.bounds)
            return group_psf
