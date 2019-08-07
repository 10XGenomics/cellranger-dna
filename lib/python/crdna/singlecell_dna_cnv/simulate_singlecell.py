import numpy as np
import crdna.constants as cnst

#...............................................................................
def is_parabola_nonnegative(linear, quadratic):
    min_gc = cnst.MIN_GC
    max_gc = cnst.MAX_GC
    gc0 = 0.45 #cnst.GC_ORIGIN
    #
    assert(gc0 > min_gc) # Important: >, not >= 
    assert(max_gc > gc0) # Important: >, not >= 
    #
    delta_left = gc0 - min_gc
    delta_right = max_gc - gc0
    quad_1 = -1.0 / (delta_left * delta_left) + linear / delta_left
    quad_2 = -1.0 / (delta_right * delta_right) - linear / delta_right
    #
    return((quadratic >= quad_1) and (quadratic >= quad_2))
# is_parabola_nonnegative

#...............................................................................
def generate_gc_coefficients(
        gc_linear_per_cell_mean, gc_quadratic_per_cell_mean,
        gc_linear_per_cell_sd, gc_quadratic_per_cell_sd,
        gc_linear_quadratic_covar, n_cells):
    #
    gc_bias_mean = np.array([gc_linear_per_cell_mean, gc_quadratic_per_cell_mean])
    gc_bias_covariance = np.array([
        [gc_linear_per_cell_sd * gc_linear_per_cell_sd, gc_linear_quadratic_covar],
        [gc_linear_quadratic_covar, gc_quadratic_per_cell_sd * gc_quadratic_per_cell_sd]])
    sampled_points = np.random.multivariate_normal(
        mean=gc_bias_mean, 
        cov=gc_bias_covariance, 
        size=n_cells)
    #
    # Make sure sampled points satisfy non-negativity condition
    # for the GC parabola on the segment [gc_min, gc_max], anchored at gc0
    for cell in xrange(n_cells):
        linear, quadratic = sampled_points[cell]
        nonnegative = is_parabola_nonnegative(linear, quadratic)
        while not nonnegative:
            trial = np.random.multivariate_normal(
                mean=gc_bias_mean, 
                cov=gc_bias_covariance, 
                size=1)
            linear, quadratic = trial[0]
            nonnegative = is_parabola_nonnegative(linear, quadratic)
            if nonnegative:
                sampled_points[cell, :] = [linear, quadratic]
            # if nonnegative
        # while not nonnegative
    # for cell
    #
    linear = sampled_points[:, 0]
    quadratic = sampled_points[:, 1]
    result = {'Linear': linear, 'Quadratic': quadratic}
    return(result)
# generate_gc_coefficients

