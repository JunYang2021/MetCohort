import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from batch_peakdetection import BatchRoi, compute_peak
from file_processing import MsFile
import multiprocessing as mp


def get_group_id(array: np.ndarray, delta_mz):
    tree = cKDTree(array.reshape(-1, 1))
    group_id = np.zeros_like(array, dtype=np.int32)
    group = 1
    for i in range(array.size):
        if group_id[i] == 0:
            ids = tree.query_ball_point(array[i], delta_mz)
            zero_mask = group_id[ids] == 0
            group_id[np.array(ids)[zero_mask]] = group
            group += 1
    return group_id


# 在代表中获得ROI
def grouping_roi_range(ranges_df: pd.DataFrame, qc_number, interfile_delta_mz=0.01):
    pad = 3
    # ranges_df (input): {'file': 'int16', 'mz': 'float32', 't0': 'float32', 't1': 'float32'}
    # 'interfile_delta_mz' can be a bit larger to ensure same roi are in one group
    ranges_df['group_id'] = get_group_id(ranges_df['mz'].values, interfile_delta_mz)
    # sort with 't0' in every same 'group_id'
    ranges_df.sort_values(['group_id', 't0'], inplace=True)
    # maximum value of previous rows of 't1' in this group
    ranges_df['t0_'] = ranges_df.groupby('group_id')['t1'].transform(lambda x: x.cummax().shift().fillna(-np.inf))
    ranges_df['t0_l_t0'] = ranges_df['t0_'] < ranges_df['t0']
    ranges_df['groupt_id'] = ranges_df['t0_l_t0'].cumsum()

    grouped = ranges_df.groupby('groupt_id')
    mzmean = grouped['mz'].mean()
    files = grouped['file'].apply(set)
    t0 = grouped['t0'].min()
    t1 = grouped['t1'].max()

    ranges_df = pd.DataFrame({'mzmean': mzmean,
                              'files': files,
                              't0': t0,
                              't1': t1}).reset_index(drop=True)
    # ranges_df (output): {'mzmean': 'float32', 't0': float32, 't1': float32}
    ranges_df = ranges_df[ranges_df['files'].map(len) >= int(qc_number * 0.5)].reset_index(drop=True)
    ranges_df['t0'] = ranges_df['t0'] - pad
    ranges_df['t1'] = ranges_df['t1'] + pad
    ranges_df.drop('files', axis=1, inplace=True)
    return ranges_df


def get_closest(mzarray, mz):
    pos = np.searchsorted(mzarray, mz)
    if pos == len(mzarray):
        res = pos - 1
    elif pos == 0:
        res = pos
    else:
        res = pos if (mzarray[pos] - mz) < (mz - mzarray[pos - 1]) else pos - 1
    return res


def extract_roi(msfile: MsFile, t0, t1, mzmean, delta_mz, ref_rt):
    # Not for reference file
    # t_array = []
    # i_array = []
    # for j, scan in enumerate(msfile.exp):
    #     cur_t = msfile.corrected_time[j]
    #     if t0 <= cur_t <= t1:
    #         t_array.append(cur_t)
    #         closest = get_closest(scan.mz, mzmean)
    #         if abs(scan.mz[closest] - mzmean) < delta_mz:
    #             i_array.append(scan.i[closest])
    #         else:
    #             i_array.append(0)
    #     elif cur_t > t1:
    #         break
    # i_array = np.interp(ref_rt, t_array, i_array, left=0, right=0)  # One-dimensional linear interpolation
    # return i_array
    t_mask = (msfile.corrected_time >= t0) & (msfile.corrected_time <= t1)
    t_array = msfile.corrected_time[t_mask]
    if t_array.size == 0:
        return np.zeros(ref_rt.shape, dtype=np.float32)
    i_array = np.zeros(t_array.shape, dtype=np.float32)
    indices = np.where(t_mask)[0]
    for h, j in enumerate(indices):
        cur_scan = msfile.exp[j]
        closest = get_closest(cur_scan.mz, mzmean)
        if abs(cur_scan.mz[closest] - mzmean) < delta_mz:
            i_array[h] = cur_scan.i[closest]
    i_array = np.interp(ref_rt, t_array, i_array, left=0, right=0)
    return i_array


# @profile
def peak_detection_roi_matrix(large_exps, range_row, file_list, ref_file, delta_mz,
                              maximum_peak_width, maximumzerofiles, min_intensity, lamda, cur_tar_cpd=None):
    mzmean, t0, t1 = range_row
    # print(mzmean, t0, t1)

    # Get matrix_rt and matrix_mz from reference file
    matrix_i = []
    matrix_rt_mask = (large_exps[ref_file].corrected_time >= t0) & (large_exps[ref_file].corrected_time <= t1)
    matrix_rt = large_exps[ref_file].corrected_time[matrix_rt_mask]
    ref_i_array = np.zeros(matrix_rt.shape, dtype=np.float32)
    matrix_mz = np.full(matrix_rt.shape, mzmean, dtype=np.float32)
    indices = np.where(matrix_rt_mask)[0]
    for h, j in enumerate(indices):
        scan = large_exps[ref_file].exp[j]
        closest = get_closest(scan.mz, mzmean)
        if abs(scan.mz[closest] - mzmean) < delta_mz:
            matrix_mz[h] = scan.mz[closest]
            ref_i_array[h] = scan.i[closest]

    # Get i_array from all the files to matrix_i
    for _file in file_list:
        if _file == ref_file:
            matrix_i.append(ref_i_array)
        else:
            matrix_i.append(extract_roi(msfile=large_exps[_file],
                                        t0=t0, t1=t1,
                                        mzmean=mzmean, delta_mz=delta_mz, ref_rt=matrix_rt))

    matrix_i = np.array(matrix_i, dtype=np.float32)

    if matrix_i.shape[1] <= 2:
        return []

    batch_roi = BatchRoi(matrix_rt,
                         matrix_i,
                         matrix_mz)
    # Feature detection
    peak_res = compute_peak(batch_roi, maximum_peak_width, maximumzerofiles, min_intensity,
                            ref_file, lamda, file_list, cur_tar_cpd)
    return peak_res


if __name__ == '__main__':
    from file_processing import MsFile

    exps = {}
    input_files_list = ['F:/Test_Peakmat/xgym data/mzml data/LC-001-1.mzML',
                        'F:/Test_Peakmat/xgym data/mzml data/LC-2-3.mzML',
                        'F:/Test_Peakmat/xgym data/mzml data/LC-001-2.mzML',
                        'F:/Test_Peakmat/xgym data/mzml data/LC-2-4.mzML',
                        'F:/Test_Peakmat/xgym data/mzml data/LC-001-3.mzML',
                        'F:/Test_Peakmat/xgym data/mzml data/LC-3-1.mzML',
                        'F:/Test_Peakmat/xgym data/mzml data/LC-3-2.mzML',
                        'F:/Test_Peakmat/xgym data/mzml data/LC-3-3.mzML',
                        'F:/Test_Peakmat/xgym data/mzml data/LC-3-4.mzML',
                        'F:/Test_Peakmat/xgym data/mzml data/LC-3-5.mzML'
                        ]
    for path in input_files_list:
        exps[path] = MsFile(path)

    range_row = (582.8608, 540.19836, 678.6957)
    import time
    t0 = time.time()
    peak_res = peak_detection_roi_matrix(exps, range_row, input_files_list,
                                         input_files_list[0], 0.01, 15, 2, 1000, 0.8, None)
    t1 = time.time()
    print(t1 - t0)
