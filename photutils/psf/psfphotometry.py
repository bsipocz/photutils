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


__all__ = ['PSFPhotometryBase', 'DAOPhotPSFPhotometry']


@six.add_metaclass(abc.ABCMeta)
class PSFPhotometryBase(object):
    @abc.abstractmethod
    def do_photometry(self, image):
        pass

class DAOPhotPSFPhotometry(PSFPhotometryBase):
    """
    This class implements the DAOPHOT algorithm proposed by Stetson
    (1987) to perform point spread function photometry in crowded fields,
    which consists basically in applying the loop FIND, GROUP, NSTAR,
    SUBTRACT, FIND until no more stars are detected or a given number of
    iterations is reached.
    """

    def __init__(self, group, bkg, psf, fitshape, find=None,
                 fitter=LevMarLSQFitter(), niters=3):
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
            parameters can be given as follows:

            - Set ``psf.psf_xname``, ``psf.psf_yname`` and
              ``psf.psf_fluxname`` to strings with the names of the respective
              psf model parameter.
            - If those attributes are not found, the names ``x_0``, ``y_0``
              and ``flux`` are assumed.

            `~photutils.psf.prepare_psf_model` can be used to prepare any 2D
            model to match these assumptions.
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

        @property
        def niters(self):
            return self._niters

        @niters.setter
        def niters(self, niters):
            if isinstance(niters, int) and niters > 0:
                self._niters = niters
            else:
                raise ValueError('niters must be an interger-valued number, '
                                 'received niters = {}'.format(niters))
        
        @property
        def fitshape(self):
            return self._fitshape

        @fitshape.setter
        def fitshape(self, fitshape):
            fitshape = np.asarray(fitshape)
            if len(fitshape) == 2:
                if np.all(fitshape) > 0:
                    if np.all(fitshape) % 2 == 1:
                        self._fitshape = fitshape
                    else:
                        raise ValueError('fitshape must be odd '
                                         'integer-valued, '
                                         'received fitshape = {}'\
                                         .format(fitshape))
                else:
                    raise ValueError('fitshape must have positive elements, '
                                     'received fitshape = {}'\
                                     .format(fitshape))
            else:
                raise ValueError('fitshape must have two dimensions, '
                                 'received fitshape = {}'.format(fitshape))

    def __call__(self, **kwargs):
        """
        Parameters
        ----------
        image : 2D array-like, `~astropy.io.fits.ImageHDU`,
        `~astropy.io.fits.HDUList`
            Image to perform photometry
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
        
        if len(kwargs) == 1:
            return self.do_photometry(kwargs['image'])
        elif len(kwargs) == 2:
            return self.do_fixed_photometry(kwargs['image'], kwargs['positions'])

    def do_photometry(self, image):
        """
        Parameters
        ----------
        image : 2D array-like, `~astropy.io.fits.ImageHDU`,
        `~astropy.io.fits.HDUList`
            Image to perform photometry
        
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

        apertures = CircularAperture((sources['xcentroid'],
                                      sources['ycentroid']),
                                     r=self.psf.sigma.value*\
                                       gaussian_sigma_to_fwhm)

        sources['flux'] = aperture_photometry(residual_image,
                                              apertures)['aperture_sum']
        n = 1
        while(n <= self.niters and len(sources) > 0):
            intab = Table(names=['id', 'x_0', 'y_0', 'flux_0'],
                          data=[sources['id'], sources['xcentroid'],
                          sources['ycentroid'], sources['flux']])
            star_groups = self.group(intab)
            tab, residual_image = self.nstar(residual_image, star_groups)
            tab['iter_detected'] = n*np.ones(tab['x_fit'].shape, dtype=np.int)
            outtab = vstack([outtab, tab])
            sources = self.find(residual_image)
            n += 1
        return outtab, residual_image

    def do_fixed_photometry(self, image, positions):
        """
        Parameters
        ----------
        image : 2D array-like, `~astropy.io.fits.ImageHDU`,
        `~astropy.io.fits.HDUList`
            Image to perform photometry
        positions: `~astropy.table.Table`
            Positions at which to *start* the fit for each object, in pixel
            coordinates. Columns 'x_0' and 'y_0' must be present.
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
            apertures = CircularAperture((positions['x_0'], positions['y_0']),
                                         r=self.psf.sigma.value*\
                                           gaussian_sigma_to_fwhm)

            positions['flux_0'] = aperture_photometry(residual_image,\
                                        apertures)['aperture_sum']

        intab = Table(names=['x_0', 'y_0', 'flux_0'],
                      data=[positions['x_0'], positions['y_0'],
                      positions['flux_0']])

        star_groups = self.group(intab)
        outtab, residual_image = self.nstar(residual_image, star_groups)
        
        return outtab, residual_image

    def get_uncertainties(self):
        """
        Return the uncertainties on the fitted parameters
        """
        pass

    def nstar(self, image, star_groups):
        """
        Fit, as appropriate, a compound or single model to the given
        `star_groups`. Groups are fitted sequentially from the smallest to
        the biggest. In each iteration, `image` is subtracted by the previous
        fitted group. 
        
        Parameters
        ----------
        image : numpy.ndarray
            Background-subtracted image.
        star_groups : `~astropy.table.Table`
            This table must contain the following columns: ``id``,
            ``group_id``, ``x_0``, ``y_0``, ``flux_0``.
            ``x_0`` and ``y_0`` are initial estimates of the centroids
            and ``flux_0`` is an initial estimate of the flux.

        Return
        ------
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
                                  fixed=self.psf.fixed)
            for i in range(len(self.star_group) - 1):
                group_psf += psf_class(sigma=self.psf.sigma.value,
                                       flux=self.star_group['flux_0'][i+1],
                                       x_0=self.star_group['x_0'][i+1],
                                       y_0=self.star_group['y_0'][i+1],
                                       fixed=self.psf.fixed)
            return group_psf
