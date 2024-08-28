import warnings
import numpy as np
from scipy import signal
from scipy.special import entr
import matplotlib.pyplot as plt
from targeted_extraction import non_mins_1d, non_maxs_1d, find_path_dp_tar_ext


class BatchRoi:
    def __init__(self, rt_array, i_array, mz_array):
        self.rt_array = rt_array  # (#points, )
        self.i_array = i_array  # (#files, #points)
        self.mz_array = mz_array  # (#points, )
        self.peak_pos = None  # a list length is the number of peaks, every element is sample * 2 array. Every line of array is start and end index

    def plot(self, k, ax):
        ax.clear()
        sample_l, points_l = self.i_array.shape
        plot_length = min(50, sample_l)
        colors = plt.cm.get_cmap('tab10', plot_length)
        for i in range(plot_length):
            import random  # plot samples randomly
            random_i = random.randint(0, sample_l - 1)  # plot samples randomly
            peak_start, peak_end = self.peak_pos[k][random_i]  # plot samples randomly
            # peak_start, peak_end = self.peak_pos[k][i]
            start_t, end_t = self.rt_array[peak_start] - 10, self.rt_array[peak_end] + 10
            rt_range = (self.rt_array >= start_t) & (self.rt_array <= end_t)
            color = colors(i)
            ax.plot(self.rt_array[rt_range], self.i_array[random_i][rt_range], linewidth=1, color=color)
            ax.fill_between(self.rt_array[peak_start:peak_end + 1], self.i_array[random_i][peak_start: peak_end + 1],
                            color=color, alpha=0.3)


def safe_mean_variance(arr):
    n = arr.shape[1]
    mean = np.mean(arr, axis=1)
    var = np.zeros(mean.shape, dtype=np.float64)
    for i in range(n):
        var += (arr[:, i] - mean) ** 2
    var /= n
    return mean, var


def row_norm_uv(arr):
    arr = arr.astype(np.float64)  # Convert to float64 to prevent overflow during computations
    row_means, row_vars = safe_mean_variance(arr)
    row_means = row_means[:, np.newaxis]
    row_stds = np.sqrt(row_vars)[:, np.newaxis]
    new_arr = (arr - row_means) / (row_stds + np.finfo(np.float64).eps)
    return new_arr.astype(np.float32)


def gradient_norm(arr):
    Z = np.zeros_like(arr)

    pos_values = arr > 0
    neg_values = arr < 0

    if np.any(pos_values):
        max_val_pos = np.max(arr[pos_values])
        Z[pos_values] = arr[pos_values] / max_val_pos  # positive values scaled to 0-1

    if np.any(neg_values):
        min_val_neg = np.min(arr[neg_values])
        Z[neg_values] = arr[neg_values] / abs(min_val_neg)  # negative values scaled to -1-0

    return Z


def gaussian_kernel(size, sigma):
    size = int(size) // 2
    x, y = np.mgrid[-size:size + 1, -size:size + 1]
    normal = 1 / (2.0 * np.pi * sigma ** 2)
    g = np.exp(-((x ** 2 + y ** 2) / (2.0 * sigma ** 2))) * normal
    return g


def sobel_filter(d2signal):
    Kx = np.array([[1, 0, -1], [2, 0, -2], [1, 0, -1]], np.float32)
    Ix = signal.convolve2d(d2signal, Kx, boundary='symm', mode='same')
    return Ix


def peak_mask(arr, res_mask, d=0):
    mask = np.ones_like(arr, dtype=bool)
    M = arr.shape[0]
    for res in res_mask:
        for i in range(M):
            peak_start, peak_end = res[i]
            mask[i, peak_start - d:peak_end + d + 1] = False
    return mask


def non_ext_suppression(intensity):
    M, N = intensity.shape
    Z = np.zeros((M, N))
    padded_intensity = np.pad(intensity, ((0, 0), (2, 2)))
    left = padded_intensity[:, :-4]
    right = padded_intensity[:, 4:]
    left_1 = padded_intensity[:, 1:-3]
    right_1 = padded_intensity[:, 3:-1]

    # Generate masks for where the center pixel is greater or less than all of its neighbors
    greater = (intensity > left) & (intensity > right) & (intensity > left_1) & (intensity > right_1) & (intensity > 0)
    lesser = (intensity < left) & (intensity < right) & (intensity < left_1) & (intensity < right_1) & (intensity < 0)

    # Apply the masks to the output array
    Z[greater] = intensity[greater]
    Z[lesser] = intensity[lesser]
    return Z


def integration(d2signal, res):
    """
    Perform integration and baseline estimation
    :param d2signal: A M*N numpy array, M is the number of samples, N is the length of signal.
    :param res: Peak position list. Every element in the list is a M*2 array representing a peak. Every row in the array
     represent the start and end point index of the peak.
    :return: baseline (M*1 array, which is baseline of M samples), area_sn (A list with the same length as res, every
    element is an M*2 array, with every line representing area and s/n)
    """
    area_sn = []
    max_indices_list = []
    M, N = d2signal.shape
    masked_array = d2signal.copy()
    for peak in res:
        for i in range(M):
            peak_start, peak_end = peak[i]
            masked_array[i, peak_start:peak_end + 1] = np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        baseline = np.nanmean(masked_array, axis=1)
        noise_level = np.nanstd(masked_array, axis=1)
    baseline[np.isnan(baseline)] = 0.1

    for peak in res:
        peak_area_sn = np.zeros((M, 2))
        peak_max_indices = np.zeros(M, dtype=np.int32)
        for i in range(M):
            peak_start, peak_end = peak[i]
            peak_region = d2signal[i, peak_start:peak_end + 1] - baseline[i]
            # 靶向提取的时候出现了空的peak_region    ValueError: attempt to get argmax of an empty sequence
            peak_max_indices[i] = np.argmax(peak_region) + peak_start
            peak_region[peak_region < 0] = 0  # Ensuring that the values are not less than 0 after baseline subtraction
            peak_area_sn[i, 0] = np.trapz(peak_region)  # Area
            peak_area_sn[i, 1] = np.max(peak_region) / (noise_level[i] + np.finfo(float).eps)  # S/N
        area_sn.append(peak_area_sn)
        max_indices_list.append(peak_max_indices)

    return area_sn, max_indices_list


def compute_entropy(d2array, lamda=0.8):
    masked_array = np.where(d2array >= 0, d2array, 0)
    masked_array = masked_array.sum(axis=0)
    p = masked_array / (masked_array.sum() + np.finfo(float).eps)
    masked_array = np.where(d2array <= 0, d2array, 0)
    masked_array = masked_array.sum(axis=0)
    n_p = masked_array / (masked_array.sum() + np.finfo(float).eps)
    return min(entr(p).sum(), entr(n_p).sum()) / np.log(d2array.shape[1]) * lamda


def find_path_dp(mI, spectral_entropy, already_res=None, already_resarray=None, penalty=0.9, t=0):
    M, N = mI.shape
    # find positive line
    mD_p = np.empty((M, N), dtype=np.int8)
    mS_p = np.empty((M, N), dtype=np.float16)
    mD_p[0, :] = 0
    mS_p[0, :] = mI[0, :]
    # find negative line
    mD_n = np.empty((M, N), dtype=np.int8)
    mS_n = np.empty((M, N), dtype=np.float16)
    mD_n[0, :] = 0
    mS_n[0, :] = mI[0, :]

    for i in range(1, M):
        left_p = np.pad(mS_p[i - 1, :-1], (1, 0), mode='edge') * penalty
        up_p = mS_p[i - 1, :]
        right_p = np.pad(mS_p[i - 1, 1:], (0, 1), mode='edge') * penalty
        mS_p[i, :] = mI[i, :] + np.maximum.reduce([up_p, left_p, right_p])
        directions_p = np.argmax([up_p, left_p, right_p], axis=0)
        mD_p[i, :] = np.array([0, -1, 1])[directions_p]

        left_n = np.pad(mS_n[i - 1, :-1], (1, 0), mode='edge') * penalty
        up_n = mS_n[i - 1, :]
        right_n = np.pad(mS_n[i - 1, 1:], (0, 1), mode='edge') * penalty
        mS_n[i, :] = mI[i, :] + np.minimum.reduce([up_n, left_n, right_n])
        directions_n = np.argmin([up_n, left_n, right_n], axis=0)
        mD_n[i, :] = np.array([0, -1, 1])[directions_n]

    res = []
    res_array = np.zeros((M, N))
    start_id, end_id = None, None
    find_pos = True
    non_maxs_ll = non_maxs_1d(mS_p[-1, :])
    non_mins_ll = non_mins_1d(mS_n[-1, :])
    for col in range(N):
        if find_pos:
            # if mS_p[-1, col] > M * 0.5: # M * 0.5 is threshold of line
            if non_maxs_ll[col] > M * spectral_entropy:
                colIdx = col
                seam = [np.clip(colIdx - 2, 0, N - 1)]
                for i in range(M - 2, -1, -1):
                    colIdx = np.clip(mD_p[i + 1, colIdx] + colIdx, 0, N - 1)
                    # colIdx = np.clip(colIdx, 0, N - 1)
                    seam.append(np.clip(colIdx - 2, 0, N - 1))
                if abs(seam[-1] - seam[0]) > 10:
                    continue
                start_id = np.array(seam[::-1])
                if end_id is None or (start_id > (end_id - 3)).all():
                    find_pos = False

        else:
            # if mS_n[-1, col] < -M * 0.5:
            if non_mins_ll[col] < -M * spectral_entropy:
                colIdx = col
                seam = [np.clip(colIdx + 2, 0, N - 1)]
                for i in range(M - 2, -1, -1):
                    colIdx = np.clip(mD_n[i + 1, colIdx] + colIdx, 0, N - 1)
                    # colIdx = np.clip(colIdx + 2, 0, N - 1)
                    seam.append(np.clip(colIdx + 2, 0, N - 1))
                if abs(seam[-1] - seam[0]) > 10:
                    continue
                end_id = np.array(seam[::-1])
                if (start_id < end_id).all():
                    find_pos = True
                    res_array[np.arange(M), start_id] = 1
                    res_array[np.arange(M), end_id] = -1
                    res.append(np.column_stack((start_id, end_id)))
    if already_res:
        old_mask = peak_mask(mI, already_res, d=0)
        for r in res:
            new_mask = peak_mask(mI, [r], d=0)
            overlap = np.logical_and(np.logical_not(old_mask), np.logical_not(new_mask))
            if not np.any(overlap):
                start_id, end_id = r[:, 0], r[:, 1]
                already_resarray[np.arange(M), start_id] = 1
                already_resarray[np.arange(M), end_id] = -1
                already_res.append(r)
                t += 1
        res = already_res
        res_array = already_resarray
    return res, res_array, t


def UntargetedEdgeDetector(d2signal, gaussian_kernel_size=7, sigma=2, lamda=0.8):
    d2signal_n = row_norm_uv(d2signal)
    d2signal_n = signal.convolve2d(d2signal_n, gaussian_kernel(gaussian_kernel_size, sigma),
                                   boundary='symm', mode='same')
    gradient = sobel_filter(d2signal_n)

    nonMaxIntensity = non_ext_suppression(gradient)
    norm_gradient = gradient_norm(nonMaxIntensity)
    res, res_array, _ = find_path_dp(norm_gradient, compute_entropy(norm_gradient, lamda))
    if res:  # 2ed detect
        mask_g = np.where(peak_mask(nonMaxIntensity, res, d=3), nonMaxIntensity, 0)
        if not np.all(mask_g == 0):

            nonMaxIntensity = gradient_norm(mask_g)
            res, res_array, t = find_path_dp(nonMaxIntensity, compute_entropy(nonMaxIntensity, lamda), res,
                                             res_array, t=1)
            if t >= 2:  # 3rd detect
                mask_g = np.where(peak_mask(mask_g, res, d=3), mask_g, 0)
                nonMaxIntensity = gradient_norm(mask_g)
                res, res_array, _ = find_path_dp(nonMaxIntensity, compute_entropy(nonMaxIntensity, lamda),
                                                 res,
                                                 res_array)

        area_sn, peak_max_indices = integration(d2signal, res)
    else:
        area_sn, peak_max_indices = 0, 0
    return res, area_sn, peak_max_indices


def TargetedEdgeDetector(d2signal, d2signal_rt, cur_tar_cpd, gaussian_kernel_size=7, sigma=2, lamda=0.8):
    d2signal_n = row_norm_uv(d2signal)
    d2signal_n = signal.convolve2d(d2signal_n, gaussian_kernel(gaussian_kernel_size, sigma),
                                   boundary='symm', mode='same')
    gradient = sobel_filter(d2signal_n)

    nonMaxIntensity = non_ext_suppression(gradient)
    norm_gradient = gradient_norm(nonMaxIntensity)
    cur_tar_cpd_rt = cur_tar_cpd['RT(s)']  # retention time of targeted compound
    cur_tar_cpd_start_rt = cur_tar_cpd_rt - 10  # beginning of targeted region
    cur_tar_cpd_end_rt = cur_tar_cpd_rt + 10  # end of targeted region

    cur_tar_index = np.searchsorted(d2signal_rt, cur_tar_cpd_rt)
    cur_tar_start_index = np.searchsorted(d2signal_rt, cur_tar_cpd_start_rt)
    cur_tar_end_index = np.searchsorted(d2signal_rt, cur_tar_cpd_end_rt)

    res, res_array = find_path_dp_tar_ext(norm_gradient,
                                          cur_tar_index,
                                          cur_tar_start_index,
                                          cur_tar_end_index)
    # 2ed detect
    # mask_g = np.where(peak_mask(nonMaxIntensity, res, d=3), nonMaxIntensity, 0)
    # if not np.all(mask_g == 0):
    #     nonMaxIntensity = gradient_norm(mask_g)
    #     res, res_array, t = find_path_dp(nonMaxIntensity, compute_entropy(nonMaxIntensity, lamda), res,
    #                                      res_array, t=1)
    #     if t >= 2:  # 3rd detect
    #         mask_g = np.where(peak_mask(mask_g, res, d=3), mask_g, 0)
    #         nonMaxIntensity = gradient_norm(mask_g)
    #         res, res_array, _ = find_path_dp(nonMaxIntensity, compute_entropy(nonMaxIntensity, lamda),
    #                                          res,
    #                                          res_array)
    area_sn, peak_max_indices = integration(d2signal, res)
    return res, area_sn, peak_max_indices


def compute_peak(batch_peakgroup: BatchRoi, maximum_peak_width, maximumzerofiles, minimum_intensity, ref_file, lamda,
                 file_list, cur_tar_cpd=None):
    # import warnings
    # warnings.simplefilter('error', RuntimeWarning)

    peak_res = []

    # 1. check if there have targeted compound
    # if tar_cpd:
    #     cur_tar_cpd = []
    #     # avg_mz = batch_peakgroup.mz_array.mean()
    #     avg_mz = batch_peakgroup.mz_array.mean()
    #     for cpd in tar_cpd:
    #         cpd_mz_range = (cpd['m/z'] - 0.01, cpd['m/z'] + 0.01)  # fixed m/z threshold
    #         mz_match = (avg_mz >= cpd_mz_range[0]) and (avg_mz <= cpd_mz_range[1])
    #         rt_match = (batch_peakgroup.rt_array[0] < cpd['RT(s)']) and (batch_peakgroup.rt_array[-1] > cpd['RT(s)'])
    #         if mz_match and rt_match:
    #             cur_tar_cpd.append(cpd)
    # else:
    #     cur_tar_cpd = []

    if not cur_tar_cpd:
        # 2. no targeted compound, untargeted detection
        res = UntargetedEdgeDetector(batch_peakgroup.i_array, lamda=lamda)
        batch_peakgroup.peak_pos, area_sn, peak_max_indices = res
        for k, peak in enumerate(batch_peakgroup.peak_pos):
            sn, max_i = (np.zeros(peak.shape[0], dtype=np.float32) for _ in range(2))
            rt, mz = None, None
            peak_res.append({})
            keep = True
            zero_file = 0
            for i in range(peak.shape[0]):
                cur_sample = file_list[i]
                peak_res[-1][cur_sample] = area_sn[k][i][0]
                if area_sn[k][i][0] == 0:
                    zero_file += 1
                max_i[i] = batch_peakgroup.i_array[i][peak_max_indices[k][i]]
                peak_start, peak_end = peak[i]
                sn[i] = area_sn[k][i][1]
                dur_time = batch_peakgroup.rt_array[peak_end] - batch_peakgroup.rt_array[peak_start]
                if cur_sample == ref_file:
                    rt = batch_peakgroup.rt_array[peak_max_indices[k][i]]
                    mz = batch_peakgroup.mz_array[peak_max_indices[k][i]]
                if dur_time > maximum_peak_width:
                    keep = False

            peak_res[-1]['rt'] = rt
            peak_res[-1]['mz'] = mz
            peak_res[-1]['group'] = (batch_peakgroup, k)
            peak_res[-1]['s/n'] = np.nanmean(
                sn)  # RuntimeWarning: Mean of empty slice peak_res[-1]['s/n'] = np.nanmean(sn)
            peak_res[-1]['max_intensity'] = np.max(max_i)
            peak_res[-1]['compound'] = None
            peak_res[-1]['tar_label'] = 1
            if keep is False or zero_file > maximumzerofiles or peak_res[-1]['max_intensity'] < minimum_intensity:
                del peak_res[-1]
    else:
        # 3. targeted compound, targeted extraction
        res = TargetedEdgeDetector(batch_peakgroup.i_array,
                                   batch_peakgroup.rt_array,
                                   cur_tar_cpd,
                                   lamda=lamda)
        # 处理res，前面的是targeted_compounds
        batch_peakgroup.peak_pos, area_sn, peak_max_indices = res
        for k, peak in enumerate(batch_peakgroup.peak_pos):
            sn, max_i = (np.zeros(peak.shape[0], dtype=np.float32) for _ in range(2))
            rt, mz = None, None
            peak_res.append({})
            for i in range(peak.shape[0]):
                cur_sample = file_list[i]
                peak_res[-1][cur_sample] = area_sn[k][i][0]
                max_i[i] = batch_peakgroup.i_array[i][peak_max_indices[k][i]]
                sn[i] = area_sn[k][i][1]
                if cur_sample == ref_file:
                    rt = batch_peakgroup.rt_array[peak_max_indices[k][i]]
                    mz = batch_peakgroup.mz_array[peak_max_indices[k][i]]
            peak_res[-1]['rt'] = rt
            peak_res[-1]['mz'] = mz
            peak_res[-1]['group'] = (batch_peakgroup, k)
            peak_res[-1]['s/n'] = np.nanmean(
                sn)  # RuntimeWarning: Mean of empty slice peak_res[-1]['s/n'] = np.nanmean(sn)
            peak_res[-1]['max_intensity'] = np.max(max_i)
            peak_res[-1]['compound'] = cur_tar_cpd['compound']
            peak_res[-1]['tar_label'] = 0       # Act as a sort key, make sure the target memeber in the former
    batch_peakgroup.mz_array = None  # release memory
    return peak_res


if __name__ == '__main__':
    pass
