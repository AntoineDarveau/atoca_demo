#
# Setup
#
# To use ATOCA (Algorithm to Treat Order ContAmination), the jwst dms needs to be properly installed.
# Also, the simulated data set needs to be downloaded. This notebook assumes that every reduction steps required before the 1d extraction have been perform.
#
# Imports
#
# Import os to make the output directory
import os
import numpy as np

from jwst import datamodels
from jwst.datamodels import dqflags
from astropy.nddata.bitmask import bitfield_to_boolean_mask

# Only import the extraction step
from jwst.extract_1d import Extract1dStep

# Import utilities to analyse the outputs
# from analysis_tools import 
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm #for better display of FITS images

#
# Set parameters for extraction
#

# Path and name to the data
filename = 'example_data/IDTSOSS_clear_noisy_single_bkgd.fits'

# Create output directory
output_dir = 'atoca_results'
try:
    os.mkdir(output_dir)
except FileExistsError:
    pass

# Set parameters given to the step.
parameters = dict()

# Override reference files for ATOCA
# If not, it will be downloaded from the crds database
# Here, since these are simulations, we have the specify the proper reference files.
ref_path = 'ref_files'
# SpecTraceModel
fname = os.path.join(ref_path, 'SOSS_ref_trace_table_SUBSTRIP256.fits')
parameters['override_spectrace'] = fname
# WaveMapModel
fname = os.path.join(ref_path, 'SOSS_ref_2D_wave_SUBSTRIP256.fits')
parameters['override_wavemap'] = fname
# SpecProfileModel
fname = os.path.join(ref_path, 'SOSS_ref_2D_profile_SUBSTRIP256.fits')
parameters['override_specprofile'] = fname
# SpecKernelModel
fname = os.path.join(ref_path, 'SOSS_ref_spectral_kernel.fits')
parameters['override_speckernel'] = fname

# The soss_atoca parameters (use atoca decontamination)  
parameters['soss_atoca'] = True

# Set the output directory and output filename
parameters['output_dir'] = output_dir  # Where to save the outputs
parameters['output_file'] = 'out'  # Name of output file
# Note that the suffix `extract1dstep` will be added,
# so the real output will be `out_extract1dstep.fits`

# To check that the decontamination went well,
# it is better to save the model used for decontamination.
# It won't be saved by default.
parameters['soss_modelname'] = 'out'
# Note that the suffix `SossExtractModel` will be added,
# so the real output will be `out_SossExtractModel.fits`

# The background subtraction implemented in Extract1dStep
# for the SOSS mode is not optimal yet.
# The best is to do this step beforehand
parameters['subtract_background'] = False  # So no bkgd substraction

# The width of the box extraction in pixels
parameters['soss_width'] = 40

# The wavelength grid used by ATOCA will eventually be an input,
# so if one has a good estimate of the incident flux, it can be
# directly specified. However, for now, the algorithm uses the
# data itself to estimate the incident flux.
# The three folowing parameters are used to build the grid.
parameters['soss_n_os'] = 2  # Minimum oversampling of the native pixel grid
parameters['soss_rtol'] = 1e-4  # Relative tolerance needed on each pixel
parameters['soss_max_grid_size'] = 20000  # Max size of the wavelength grid

# The Tikhonov regularization factor can be directly specified.
# The higher it is, the more the solution will be regularized (smoothed).
# If it is not specified, the algorithm will try too find the best estimate.
parameters['soss_tikfac'] = None

# For the box extraction at the end, a value needs to be assigned to the bad pixels.
# It can be modeled or masked here. If modeled, the model from ATOCA is used.
# This is the same model that is output by `soss_modelname`
parameters['soss_bad_pix'] = 'model'

# The reference files can be shifted or rotated if necessary.
rotation, column_shift, row_shift = [0, 0, 0]
parameters['soss_transform'] = [rotation, column_shift, row_shift]

#
# Run extraction
#
# Note: This step can be relatively long, few minutes, especially when testing for different factors

# Run extraction
result = Extract1dStep().call(filename, **parameters)

# Note: We can see here that it is relatively long to converge to a good value of the tikhonov factor. However, this step is done only once for the first integration. The same factor is kept for the next integrations.
