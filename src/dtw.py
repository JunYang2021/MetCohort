from fastdtw import fastdtw
import numpy as np
import statsmodels.api as sm


def dtw(rt1: np.ndarray, spec1, rt2: np.ndarray, spec2):
    # File 1 is reference, file 2 is target file
    _, warp_path = fastdtw(spec1, spec2, dist=2)

    warp_path = np.array(warp_path)
    indices1, indices2 = warp_path[:, 0], warp_path[:, 1]

    aligned_rt1_array = rt1[indices1]
    aligned_rt2_array = rt2[indices2]

    return aligned_rt1_array, aligned_rt2_array


def loess_dtw(aligned_rt1, aligned_rt2):
    # Non-linear fitting: aligned_rt2 -> aligned_rt1
    # 1. solve the repeat values in aligned_rt2
    unique_rt2, indices = np.unique(aligned_rt2, return_inverse=True)
    averaged_rt1 = np.zeros_like(unique_rt2)
    for i in range(len(unique_rt2)):
        averaged_rt1[i] = np.mean(aligned_rt1[indices == i])
    # 2. loess regression (return loess function and MAE)
    lowess = sm.nonparametric.lowess
    loess_result = lowess(averaged_rt1, unique_rt2, frac=0.1)

    fitted_rt1 = loess_result[:, 1]
    loess_function = lambda x: np.interp(x, loess_result[:, 0], fitted_rt1)

    return loess_function


if __name__ == '__main__':
    pass
