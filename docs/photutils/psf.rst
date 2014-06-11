PSF photometry
==============

.. warning::
   The psf photometry API is currently *experimental*
   and may change in the future. For example, the functions currently
   accept `~numpy.ndarray` objects for the parameters ``data``. They may be
   changed to accept `astropy.nddata` objects.


A Simple Example
----------------

Suppose there are 4 sources located at (10, 10), (20, 20), and (30,
30), (40, 40), in pixel coordinates.

    >>> import numpy as np
    >>> import photutils
    >>> data = np.ones((100, 100))
    >>> xc = [10., 20., 30., 40.]
    >>> yc = [10., 20., 30., 40.]



Creating PRFs from image data
-----------------------------

   >>> from photutils.psf import create_prf
   >>> prf_discrete = create_prf(data, coords, 7, fluxes=fluxes_catalog,
   ...                           mask=np.logical_not(mask), subsampling=5)



Photometry with DiscretePRF
---------------------------



   >>> from photutils.psf import psf_photometry
   >>> fluxes_photutils = psf_photometry(data_slice, coords_slice, prf_discrete)


Photometry with GaussianPSF
---------------------------

   >>> from photutils.psf import GaussianPSF
   >>> psf_gaussian = GaussianPSF(1)
   >>> fluxes_gaussian = psf_photometry(data_slice, coords_slice, psf_gaussian)


Making residual images
----------------------

   >>> from photutils.psf import remove_prf
   >>> residuals = remove_prf(data_slice.copy(), prf_discrete, coords_slice, fluxes_photutils)


Reference/API
--------------

.. automodapi:: photutils.psf
    :no-heading:
