import numpy as np
from scipy.signal import find_peaks


def non_maxs_1d(row):
    Z = np.zeros(row.shape)
    peaks = find_peaks(row)[0]
    Z[peaks] = row[peaks]
    return Z


def non_mins_1d(row):
    Z = np.zeros(row.shape)
    peaks = find_peaks(-row)[0]
    Z[peaks] = row[peaks]
    return Z


def find_path_dp_tar_ext(mI,
                         tar_index,
                         start_region_index,
                         end_region_index,
                         penalty=0.9):
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
    non_maxs_ll = non_maxs_1d(mS_p[-1, :])
    non_mins_ll = non_mins_1d(mS_n[-1, :])

    # for indx, start_indx, end_indx in zip(tar_index, start_region_index, end_region_index):
    # find peak start for one targeted compound, maximum seam in before 10 seconds is regarded as beginning seam
    # in non_maxs_ll[0, indx-1], find closest non-zero value from indx.
    if not start_region_index < tar_index < end_region_index - 1:
        return res, res_array
    valid_peak_start = start_region_index + np.argmax(non_maxs_ll[start_region_index:tar_index])
    if valid_peak_start > N - 3:
        valid_peak_start = N - 3
    colIdx = valid_peak_start
    seam = [np.clip(colIdx - 2, 0, N - 5)]
    for i in range(M - 2, -1, -1):
        colIdx = np.clip(mD_p[i + 1, colIdx] + colIdx, 0, N - 3)
        seam.append(np.clip(colIdx - 2, 0, N - 5))
    start_id = np.array(seam[::-1])
    # find peak end for one targeted compound, minimum sean in after 10 seconds is regarded as end seam
    valid_peak_end = end_region_index - np.argmin(non_mins_ll[tar_index + 1:end_region_index][::-1])
    if valid_peak_end > N - 1:
        valid_peak_end = N - 1
    colIdx = valid_peak_end
    seam = [np.clip(colIdx + 2, 0, N - 1)]
    for i in range(M - 2, -1, -1):
        colIdx = np.clip(mD_n[i + 1, colIdx] + colIdx, 0, N - 1)
        seam.append(np.clip(colIdx + 2, 0, N - 1))
    end_id = np.array(seam[::-1])

    if not (start_id < end_id).all():
        end_id = np.clip(start_id + (valid_peak_end - valid_peak_start), 0, N - 1)
    res_array[np.arange(M), start_id] = 1
    res_array[np.arange(M), end_id] = -1
    res.append(np.column_stack((start_id, end_id)))
    return res, res_array
