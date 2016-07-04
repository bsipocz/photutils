# Licensed under a 3-clause BSD style license - see LICENSE.rst
class PSFPhotometryBase(object):
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod
    def do_photometry(self):
        pass

class NStarPSFPhotometry(PSFPhotometryBase):
    """
    This is an implementation of the NSTAR algorithm proposed by Stetson
    (1987) to perform point spread function photometry in crowded fields.

    This class basically implements the loop FIND, GROUP, NSTAR,
    SUBTRACT, FIND until no more stars are detected.
    """

    def __init__(self, find, group, bkg, psf, fitter, niters, fitshape):
        """
        Attributes
        ----------
        find : an instance of any StarFinderBase subclasses
        group : an instance of any GroupStarsBase subclasses
        bkg : an instance of any BackgroundBase2D (?) subclasses
        psf : Fittable2DModel instance
        fitter : Fitter instance
        niters : int
            number of iterations for the loop FIND, GROUP, SUBTRACT, NSTAR
        fitshape : array-like
            rectangular shape around the center of a star which will be used
            to collect the data to do the fitting
        """

        self.find = find
        self.group = group
        self.bkg = bkg
        self.psf = psf
        self.fitter = fitter
        self.niters = niters
        self.fitshape = fitshape

    @property
    def find(self):
        return self._find
    
    @find.setter
    def find(self, find):
        if isinstance(find, StarFinderBase)
            self._find = find
        else:
            raise ValueError('find is expected to be an instance of '
                             'StarFinderBase, received {}'.format(type(find)))
    @property
    def group(self):
        return self._group

    @group.setter
    def group(self, group):
        if isinstance(group, GroupStarsBase):
            self._group = group
        else:
            raise ValueError('group is expected to be an instance of '
                             'GroupStarsBase, received {}.'.\
                             format(type(group)))
    
    @property
    def bkg(self)
        return self.bkg

    @background.setter
    def bkg(self, bkg):
        if isinstance(bkg, BackgroundBase2D):
            self._background = background
        else:
            raise ValueError('bkg is expected to be an instance of '
                             'BackgroundBase2D, received {}.'\
                             .format(type(bkg)))
    
    @property
    def psf(self):
        return self._psf

    @psf.setter
    def psf(self, psf):
        if isinstance(psf, Fittable2DModel):
            self._psf = psf_model
        else:
            raise ValueError('psf_model is expected to be an instance of '
                             'Fittable2DModel, received {}.'\
                             .format(type(psf)))

    @property
    def fitter(self):
        return self._fitter

    @fitter.setter
    def fitter(self, fitter):
        if isinstance(fitter, Fitter):
            self._fitter = fitter
        else:
            raise ValueError('fitter is not a valid astropy Fitter, '
                             'received {}'.format(type(fitter)))
    
    @property
    def niters(self):
        return self._niters

    @niters.setter
    def niters(self, niters):
        if isinstance(niters, int) and niters > 0:
            self._niters = niters
        else:
            raise ValueError('niters is not defined properly, '
                             'received niters = {}'.format(niters))

    @property
    def fitshape(self):
        return self._fitshape

    @fitshape.setter
    def fitshape(self, fitshape):
        fitshape = np.asarray(fitshape)
        if len(fitshape) == 2 and np.all(fitshape) > 0:
            self._fitshape = fitshape
        else:
            raise ValueError('fitshape is not defined properly, '
                             'received fitshape = {}'.format(fitshape))

    def __call__(self, image):
        """
        Parameters
        ----------
        image : array-like, ImageHDU, HDUList
            image to perform photometry
        
        Returns
        -------
        outtab : astropy.table.Table
            Table with the photometry results, i.e., centroids and flux
            estimations.
        residual_image : array-like, ImageHDU, HDUList
            Residual image calculated by subtracting the fitted sources
            and the original image.
        """

        return self.do_photometry(image)

    def _nstar(self, image, star_groups):
        """
        Fit, as appropriate, a compound or single model to a given `groups` of
        stars. Groups are fitted sequentially from the smallest to the biggest. In
        each iteration, `image` is subtracted by the previous fitted group. 
        
        Parameters
        ----------
        image : numpy.ndarray
            Background-subtracted image.
        star_groups : list of `~astropy.table.Table`
            Each `~astropy.table.Table` in this list corresponds to a group of
            mutually overlapping starts.

        Return
        ------
        result_tab : `~astropy.table.Table`
            Astropy table that contains the results of the photometry.
        image : numpy.ndarray
            Residual image.
        """
        
        result_tab = Table([[], [], [], [], []],
                           names=('id', 'group_id', 'x_fit', 'y_fit',
                                  'flux_fit'),
                           dtype=('i4', 'i4', 'f8', 'f8', 'f8'))

        star_groups = star_groups.group_by('group_id')
        
        groups_order = []
        for g in star_groups.groups:
            groups_order.append(len(g))

        while len(groups_length) > 0:
            curr_order = np.min(groups_order)
            n = 0
            N = len(groups_order)
            while(n < N):
                if curr_order == len(star_groups.groups[n]):
                    group_psf = _create_sum_psf_model(star_groups.groups[n])
                    x, y, data = _extract_shape_and_data(self.fitshape,
                                                         star_groups.groups[n],
                                                         image)
                    fit_model = self.fitter(self.psf, x, y, data)
                    param_table = _model_params2table(fit_model,
                                                      star_groups.groups[n])
                    result_tab = vstack([result_tab, param_table])
                    image = subtract_psf(image, self.psf, param_table)
                    N = N - 1
                n += 1
        return result_tab, image

    def _model_params2table(fit_model, star_group):
        """
        Place fitted parameters into an astropy table.
        
        Parameters
        ----------
        fit_model : Fittable2DModel
        star_group : ~astropy.table.Table
        
        Returns
        -------
        param_tab : ~astropy.table.Table
            Table that contains the fitted parameters.
        """

        param_tab = Table([[], [], [], [], []],
                          names=('id', 'group_id', 'x_fit','y_fit','flux_fit'),
                          dtype=('i4', 'i4', 'f8', 'f8', 'f8'))
        for i in range(np.size(fit_model)):
            param_table.add_row([[star_group['id'][i]],
                                 [star_group['group_id'][i]],
                                 [getattr(fit_model,'x_0_'+str(i)).value],
                                 [getattr(fit_model, 'y_0_'+str(i)).value],
                                 [getattr(fit_model, 'flux_'+str(i)).value]])
        return param_tab

    def _extract_shape_and_data(shape, star_group, image):
        """
        Parameters
        ----------
        shape : tuple
            Shape of a rectangular region around the center of an isolated source.
        star_group : `astropy.table.Table`
            Group of stars
        image : numpy.ndarray

        Returns
        -------
        x, y : numpy.mgrid
            All coordinate pairs (x,y) in a rectangular region which encloses all
            sources of the given group
        image : numpy.ndarray
            Pixel value
        """

        xmin = int(np.around(np.min(star_group['x_0'])) - shape[0])
        xmax = int(np.around(np.max(star_group['x_0'])) + shape[0])
        ymin = int(np.around(np.min(star_group['y_0'])) - shape[1])
        ymax = int(np.around(np.max(star_group['y_0'])) + shape[1])
        y, x = np.mgrid[ymin:ymax+1, xmin:xmax+1]

        return x, y, image[ymin:ymax+1, xmin:xmax+1]

    def _create_sum_psf_model(self, star_group):
        """
        Construct a joint psf model which consists in a sum of `self.psf`
        whose parameters are given in `star_group`.

        Parameters
        ----------
        star_group : `~astropy.table.Table`
            Table from which the compound PSF will be constructed.
            It must have columns named as `x_0`, `y_0`, and `flux_0`.
        
        Returns
        -------
        sum_psf : CompoundModel
            `CompoundModel` instance which is a sum of the given PSF
            models.
        """

        sum_psf = self.psf(sigma=self.psf.sigma.value,
                           flux=star_group['flux_0'][0],
                           x_0=star_group['x_0'][0], y_0=star_group['y_0'][0])
        for i in range(len(group) - 1):
            sum_psf += self.psf(sigma=self.psf.sigma.value,
                                flux=star_group['flux_0'][i+1],
                                x_0=star_group['x_0'][i+1],
                                y_0=star_group['y_0'][i+1])
        return sum_psf
    
    def do_photometry(self, image):
        # prepare output table
        outtab = Table([[], [], [], [], [], []],
                       names=('id', 'group_id', 'x_fit', 'y_fit', 'flux_fit',
                              'iter_detected'),
                       dtype=('i4', 'i4', 'f8', 'f8', 'f8', 'i4'))

        # make a copy of the input image
        residual_image = image.copy()

        # find potential sources on the given image
        sources = self.find(residual_image)

        n = 1
        # iterate until no more sources are found or the number of iterations
        # has been reached
        while(n <= self.niters and len(sources) > 0):
            # prepare input table
            intab = Table(names=['id', 'x_0', 'y_0', 'flux_0'],
                          data=[sources['id'], sources['xcentroid'],
                          sources['ycentroid'], sources['flux']])

            # find groups of overlapping sources
            star_groups = self.group(intab)

            # fit the sources within in each group in a simultaneous manner
            # and get the residual image
            tab, residual_image = self.nstar(residual_image, groups)

            # mark in which iteration those sources were fitted
            tab['iter_detected'] = np.int(n*np.ones(tab['x_fit'].shape))

            # populate output table
            outtab = vstack([outtab, tab])
            
            # find remaining sources in the residual image
            sources = self.find(residual_image)
            n += 1
        return outtab, residual_image
