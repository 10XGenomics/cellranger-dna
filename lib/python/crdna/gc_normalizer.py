import time
import numpy as np
from itertools import product
from crdna.constants import (
    MIN_CURVE, MAX_CURVE, MIN_GC, MAX_GC, GC_ORIGIN, BREAKPOINT_READ_THRESHOLD)
from plot_utils import aggregate_counts


def gc_curve(g, g2, linear, quadratic):
    return np.clip(1.0 + np.outer(linear, g) + np.outer(quadratic, g2),
                   MIN_CURVE,
                   MAX_CURVE)


def generate_grid_points(bounds, num):
    l = np.linspace(*bounds[:, 0], num=num+1, dtype=float)
    q = np.linspace(*bounds[:, 1], num=num+1, dtype=float)
    grid = np.array(list(product(l, q)), dtype=float)
    g = np.array([MIN_GC, MAX_GC], dtype=float) - GC_ORIGIN
    g2 = g * g
    mask = np.all(gc_curve(g, g2, grid[:, 0], grid[:, 1]) > 0, axis=1)
    return grid[mask]


def entropy(density, delta_x):
    eps = 1e-7
    e = -np.nansum(density * np.log(density + eps)) * delta_x
    return e


def minimize_entropy(profile, gc):
    count_ceiling = max(700, int(np.percentile(profile, 98)))
    selector, = np.where(
        (MIN_GC <= gc) & (gc < MAX_GC) &
        (profile < count_ceiling))
    profile = profile[selector]
    gc = gc[selector]
    g = gc - GC_ORIGIN
    g2 = g * g
    FLOAT_MAX = np.finfo(float).max
    best_entropy = FLOAT_MAX
    best_linear = 0.0
    best_quadratic = 0.0
    bins = min(count_ceiling, 140)
    bin_width = float(count_ceiling) / bins
    bounds = np.array([[-10, -20], [10, 20]], dtype=float)
    grid_n = 20
    # this controls regularization strength
    lam = 5e-3
    tmp = np.zeros((profile.shape[0],), dtype=bool)

    def f(x):
        regularizer = lam / 2.0 * np.linalg.norm(x)
        linear, quadratic = x
        parabola = np.clip(1.0 + linear*g + quadratic*g2, MIN_CURVE, MAX_CURVE)
        parabola *= bin_width
        norm_profile = profile / parabola

        if (np.isnan(norm_profile, out=tmp).any() or
            np.less(norm_profile, 0, out=tmp).any()):
            return FLOAT_MAX

        norm_profile = norm_profile[norm_profile < bins]
        frequency = np.bincount(norm_profile.astype(int),
                                minlength=bins).astype(float)
        frequency /= frequency.sum() * bin_width
        return entropy(frequency, bin_width) + regularizer

    # recursive grid search
    for _ in xrange(3):
        grid = generate_grid_points(bounds, grid_n)
        for i in xrange(grid.shape[0]):
            e = f(grid[i, :])
            if e < best_entropy:
                best_entropy = e
                best_linear, best_quadratic = grid[i, :]
        # regenerate bounds for next iteration
        l, q = (bounds[1, :] - bounds[0, :]) / [grid_n, grid_n] / 0.5
        bounds = np.array([[best_linear-l, best_quadratic-q],
                           [best_linear+l, best_quadratic+q]], dtype=float)

    return (best_entropy, best_linear, best_quadratic)


def estimate_window_size(profile):
    mean_count = np.nanmean(profile.astype(float))
    target_count = float(BREAKPOINT_READ_THRESHOLD)
    window_size = np.ceil(target_count / mean_count)
    if np.isfinite(window_size):
        return max(50, int(window_size))
    return profile.shape[0]


def get_cv(profile):
    std = np.nanstd(profile)
    mean = np.nanmean(profile)
    if mean <= 0:
        return np.nan
    return std / mean


def get_delta_gc(raw_cv, norm_cv):
    tmp = raw_cv * raw_cv - norm_cv * norm_cv
    delta_gc = np.sign(tmp) * np.sqrt(np.abs(tmp))
    return delta_gc


def get_delta_gc_cv(raw_profile, norm_profile, gc):
    mask = np.logical_and(MIN_GC <= gc, gc <= MAX_GC)
    raw_cv = get_cv(raw_profile[mask])
    norm_cv = get_cv(norm_profile[mask])
    delta_gc_cv = get_delta_gc(raw_cv, norm_cv)
    return delta_gc_cv


def estimate_gc_norm_cell(profile, gc, mask):
    window_size = estimate_window_size(profile[mask])
    profile = aggregate_counts(profile[mask].astype(float), window_size)
    gc = aggregate_counts(gc[mask].astype(float), window_size) / window_size
    _, linear, quadratic = minimize_entropy(profile, gc)
    g = gc - GC_ORIGIN
    g2 = g * g
    parabola = gc_curve(g, g2, linear, quadratic).ravel()
    norm_profile = profile / parabola

    delta_gc_cv = get_delta_gc_cv(profile, norm_profile, gc)
    if delta_gc_cv < 0:
        linear = 0.0
        quadratic = 0.0
        delta_gc_cv = 0.0

    return linear, quadratic


def estimate_gc_normalization(profiles, gc, mask):
    ncells, _ = profiles.shape
    scale = [1.0] * ncells
    linear = [0.0] * ncells
    quadratic = [0.0] * ncells
    for cell in xrange(ncells):
        tick = time.time()
        result = estimate_gc_norm_cell(profiles[cell, :], gc, mask)
        linear[cell], quadratic[cell] = result
        tock = time.time()
        duration = tock - tick
        print "completed cell {}/{}, linear={:+.3f}, quadratic={:+.3f}, time={:.3f}".format(
                cell + 1, ncells, linear[cell], quadratic[cell], duration)
    return scale, linear, quadratic


def gc_normalize_chrom(profiles, gc, linear, quadratic):
    g = gc - GC_ORIGIN
    g2 = g * g
    parabola = gc_curve(g, g2, linear, quadratic)
    assert profiles.shape[0] == parabola.shape[0]
    norm_profile = profiles / parabola
    norm_profile[norm_profile < 0] = 0.0
    return norm_profile


def gc_normalize(profiles, gc, linear, quadratic, chroms):
    linear = np.array(linear, dtype=float)
    quadratic = np.array(quadratic, dtype=float)
    norm_profiles = []
    for i, chrom in enumerate(chroms):
        norm_profiles.append(gc_normalize_chrom(profiles[i], gc[i], linear, quadratic))
    return norm_profiles
