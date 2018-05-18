from skimage.io import imread
from skimage.io import imsave
import matplotlib.pyplot as plt
import numpy as np
import os
from BilateralGrid import *
from BilateralSolver import *

data_folder = 'E:\\Project\\FastBilateralSolver\\data'

# The RGB image that whose edges we will respect
reference = imread(os.path.join(data_folder, 'reference.png'))
# The 1D image whose values we would like to filter
target = imread(os.path.join(data_folder, 'target.png'))
# A confidence image, representing how much we trust the values in "target".
# Pixels with zero confidence are ignored.
# Confidence can be set to all (2^16-1)'s to effectively disable it.
confidence = imread(os.path.join(data_folder, 'confidence.png'))

im_shape = reference.shape[:2]
assert(im_shape[0] == target.shape[0])
assert(im_shape[1] == target.shape[1])
assert(im_shape[0] == confidence.shape[0])
assert(im_shape[1] == confidence.shape[1])

grid_params = {
    'sigma_luma' : 4, # Brightness bandwidth
    'sigma_chroma': 4, # Color bandwidth
    'sigma_spatial': 8 # Spatial bandwidth
}

bs_params = {
    'lam': 128, # The strength of the smoothness parameter
    'A_diag_min': 1e-5, # Clamp the diagonal of the A diagonal in the Jacobi preconditioner.
    'cg_tol': 1e-5, # The tolerance on the convergence in PCG
    'cg_maxiter': 25 # The number of PCG iterations
}

grid = BilateralGrid(reference, **grid_params)

t = target.reshape(-1, 1).astype(np.float64) / (pow(2,16)-1)
c = confidence.reshape(-1, 1).astype(np.float64) / (pow(2,16)-1)
tc_filt = grid.filter(t * c)
c_filt = grid.filter(c)
output_filter = (tc_filt / c_filt).reshape(im_shape)

output_solver = BilateralSolver(grid, bs_params).solve(t, c).reshape(im_shape)

imsave('refined_depth.png', output_filter)

plt.figure('bilateral solver')
plt.imshow(output_filter)
plt.show()