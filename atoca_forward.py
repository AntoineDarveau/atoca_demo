import numpy as np

from jwst.extract_1d.soss_extract.soss_extract import get_ref_file_args, get_trace_1d
from jwst.extract_1d.soss_extract.atoca_utils import make_combined_adaptive_grid, grid_from_map
from jwst.extract_1d.soss_extract.atoca import ExtractionEngine


def forward_model(input_flux, masked_pixels, wave_map, trace_profile, throughput,
                  kernel=None, wave_grid=None, rtol=1e-3,
                  max_grid_size=10000, n_os=2, spectral_order=1):
    """Perform the spectral extraction on a single image.

    Parameters
    ----------
    input_flux : callable or 1d array
        flux to be projected on the detector. Can be a function (callable)
        or a 1d array. If a 1d array is given, the corresponding wavelength
        grid, `wave_grid`, must be specified. If callable, it will be projected
        on `wave_grid`.
    masked_pixels : array[bool]
        Pixel mask to apply to detector image. Set to True if pixel is masked.
    wave_map : (n_row, n_col) 2-D array
        2d-array of the central wavelength position for a given
        order on the detector.
    trace_profile : (n_row, n_col) 2-D array
        2d-array of the spatial profile for a given order
        on the detector.
    throughput : ([, N_k]) list of array or callable
        A list of functions or array of the throughput at each order.
        If callable, the functions depend on the wavelength.
        If array, projected on `wave_grid`.
    kernel : array, callable or sparse matrix
        Convolution kernel to be applied on the spectrum (f_k) for a given order.
        Can be array of the shape (N_ker, N_k_c).
        Can be a callable with the form f(x, x0) where x0 is
        the position of the center of the kernel. In this case, it must
        return a 1D array (len(x)), so a kernel value
        for each pairs of (x, x0). If array or callable,
        it will be passed to `atoca_utils.get_c_matrix` function
        and the `c_kwargs` can be passed to this function.
        If sparse, the shape has to be (N_k_c, N_k) and it will
        be used directly. N_ker is the length of the effective kernel
        and N_k_c is the length of the spectrum (f_k) convolved.
    wave_grid : (N_k) array_like, optional
        The wavelength grid on which `input_flux` will be projected. If not given,
        a grid is generated using `soss_extract.make_combined_adaptive_grid`
        with `n_os`, `r_tol`, `max_grid_size` and `input_flux` as input.
    rtol: float, optional
        The desired relative tolerance. Default is 1e-3.
    max_total_size: int, optional
        maximum size of the output grid. Default is 10 000.
    n_os : int, optional
        The minimum oversampling factor to generate the wavelength grid.
        If not specified, defaults to 2.
    spectral_order: int
        The spectral order for which to return the trace parameters.

    Returns
    -------
    tracemodels : 2d array
        Modeled detector images for a given spectral order.
    engine : atoca ExtractionEngine object
        engine used to create image simulations. If nothing has
        changed in the reference files, can be used to project
        flux again on image with `image = engine.rebuild(new_flux)`.
    """

    if kernel is None:
        # No convolution needed (so give equivalent of identity)
        kernel = np.array([1.])
        c_kwargs = None
    else:
        # Set the convolution kwargs using the minimum value of the kernel
        c_kwargs = {'thresh': kernel.min_value}

    # Generate grid based on flux function and tolerance (if grid not given)
    if wave_grid is None:
        wave_grid = grid_from_map(wave_map, trace_profile, n_os=n_os)
        wave_grid = make_combined_adaptive_grid([wave_grid], [input_flux], rtol=1e-3,
                                                max_total_size=max_grid_size)

    # Initialize the Engine.
    ref_file_args = ([wave_map], [trace_profile], [throughput], [kernel])
    engine = ExtractionEngine(*ref_file_args,
                              wave_grid=wave_grid,
                              mask_trace_profile=[masked_pixels],
                              c_kwargs=c_kwargs,
                              orders=[spectral_order])

    # Build image
    tracemodel = engine.rebuild(input_flux, fill_value=0)

    return tracemodel, engine


# ########################################################################
# ################## Example code using the forward model ################
# ########################################################################

# Imports to get the appropriate reference files
from jwst.extract_1d import Extract1dStep
from jwst import datamodels

# Import to generate fake flux
import numpy as np

# Import to generate aperture based mask
from jwst.extract_1d.soss_extract.soss_boxextract import get_box_weights
# Imports to generate the reference pixel mask
from jwst.datamodels.dqflags import pixel
from astropy.nddata.bitmask import bitfield_to_boolean_mask


# ################################
# ######## Parameters ############
# ################################

# Which order is simulated
spectral_order = 1

# Parameter to get appropriate reference files
transform = [0,0,0]

# Need an observation example to get the appropriate reference files
observation_name = 'simulations/timeseries_20210729_forGJ/nint1/IDTSOSS_clear_noisy_rate.fits'

# Define a fake flux function
def fake_flux_fct(wavelength):
    return (10 + np.sin(wavelength / 0.0005)) * 1e5

flux_fct = fake_flux_fct

# ################################
# #### Get reference inputs ######
# ################################

step = Extract1dStep()
ref_files = dict()

# Get the reference files model class
ref_datamodel = {'spectrace': datamodels.SpecTraceModel,
                 'wavemap': datamodels.WaveMapModel,
                 'specprofile': datamodels.SpecProfileModel,
                 'speckernel': datamodels.SpecKernelModel}

# Read the reference files
for ref_type in ['spectrace', 'wavemap', 'specprofile', 'speckernel']:
    ref_filename = step.get_reference_file(observation_name, ref_type)
    ref_files[ref_type] = ref_datamodel[ref_type](ref_filename)

# Prepare the reference file arguments.
wave_maps, trace_profiles, throughputs, kernels = get_ref_file_args(ref_files, transform)

# Take only the order's specific ref_files
wave_map = wave_maps[spectral_order - 1]
trace_profile = trace_profiles[spectral_order - 1]
throughput = throughputs[spectral_order - 1]
kernel = kernels[spectral_order - 1]

# ################################
# ####### Generate mask ##########
# ################################

# Example of a mask of pixel to simulate based on box aperture
width = 40  # In pixels
xtrace, ytrace, _ = get_trace_1d(ref_files, transform, spectral_order)
box_weights = get_box_weights(ytrace, width, trace_profile.shape, cols=xtrace)
masked_pixels = ~(box_weights > 0)

# Add reference pixels to the mask (based on observation)
obs_datamodel = datamodels.open(observation_name)
ref_pix = bitfield_to_boolean_mask(obs_datamodel.dq,
                                   ignore_flags=pixel['REFERENCE_PIXEL'],
                                   flip_bits=True)
masked_pixels |= ref_pix[0]  # Take first integration

# ##############################
# ####### Build image ##########
# ##############################

# Example 1: Generate a single image with input flux as callable
# note: the computation time is sensitive to
#       - the total number of non-masked pixel (~masked_pixels)
#       - the length of the wavelength grid (wave_grid)
#         which is influenced by rtol, max_grid_size and n_os.
image, atoca_engine = forward_model(flux_fct, masked_pixels, wave_map, trace_profile, throughput,
                                    kernel=None, wave_grid=None, rtol=1e-3,
                                    max_grid_size=20000, n_os=2, spectral_order=spectral_order)

# Example 2: Generate image with flux and grid arrays
grid = atoca_engine.wave_grid
image, _ = forward_model(flux_fct(grid), masked_pixels, wave_map, trace_profile, throughput,
                      kernel=None, wave_grid=grid, spectral_order=spectral_order)

# Example 3: Generate second image without re-creating the engine
# note: Useful if the reference files did not change to save time
image = atoca_engine.rebuild(flux_fct)

# Example 4: Convolve the flux to the orders resolution (takes more time)
grid = atoca_engine.wave_grid
image, _ = forward_model(flux_fct(grid), masked_pixels, wave_map, trace_profile, throughput,
                         kernel=kernel, wave_grid=grid, spectral_order=spectral_order)

