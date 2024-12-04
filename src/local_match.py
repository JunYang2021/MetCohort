import numpy as np
from scipy import signal
from scipy.stats import t
import statsmodels.api as sm
from scipy.interpolate import interp1d
import pandas as pd
from scipy.interpolate import PchipInterpolator


def preprocess_pairs(res1, res2):
    combined = np.vstack((res2, res1)).T
    combined = combined[combined[:, 0].argsort()]
    unique_res2, unique_indices = np.unique(combined[:, 0], return_index=True)
    means_res1 = np.array([combined[combined[:, 0] == val, 1].mean() for val in unique_res2])
    return means_res1, unique_res2


def local_match(file1_roas, file2_eics, loess_function, per=68.27):
    # return f_rt, f_mz (file2 -> file1)
    # loess_function: file1 -> file2
    global outlier_indices, residual_cutoff
    file1_rt, file1_mz, file2_rt, file2_mz = [], [], [], []
    file1_rt_u, file1_mz_u, file2_rt_u, file2_mz_u = [], [], [], []
    pn = 5  # Candidate peak numbers is 5
    min_rt, max_rt = file2_eics['rt'][0], file2_eics['rt'][-1]

    for i in range(len(file1_roas)):
        if np.max(file2_eics['intensity'][i]) < file1_roas[i].i_max / 10:
            continue
        corr = signal.correlate(file2_eics['intensity'][i], file1_roas[i].i, mode='same')
        peaks, _ = signal.find_peaks(corr)
        sorted_peaks = sorted(peaks, key=lambda h: corr[h], reverse=True)
        if len(sorted_peaks) >= 2 and corr[sorted_peaks[0]] > corr[sorted_peaks[1]] * 5:
            file1_rt.append(file1_roas[i].rt[file1_roas[i].i_max_id])
            file1_mz.append(file1_roas[i].mzmean)
            file2_rt.append(file2_eics['rt'][sorted_peaks[0]])
            file2_mz.append(file2_eics['mz'][i][sorted_peaks[0]])
        elif len(sorted_peaks) == 1:
            file1_rt.append(file1_roas[i].rt[file1_roas[i].i_max_id])
            file1_mz.append(file1_roas[i].mzmean)
            file2_rt.append(file2_eics['rt'][sorted_peaks[0]])
            file2_mz.append(file2_eics['mz'][i][sorted_peaks[0]])
        elif len(sorted_peaks) >= 2:
            file1_rt_u.append(file1_roas[i].rt[file1_roas[i].i_max_id])
            file1_mz_u.append(file1_roas[i].mzmean)
            file2_rt_u.append(
                [file2_eics['rt'][sorted_peaks[j]] if j < len(sorted_peaks) else -1000 for j in range(pn)])
            file2_mz_u.append(
                [file2_eics['mz'][i][sorted_peaks[j]] if j < len(sorted_peaks) else -1000 for j in range(pn)])
    file1_rt, file2_rt, file1_mz, file2_mz = np.array(file1_rt), np.array(file2_rt), np.array(file1_mz), np.array(
        file2_mz)
    file1_rt_u, file2_rt_u, file1_mz_u, file2_mz_u = np.array(file1_rt_u), np.array(file2_rt_u), np.array(
        file1_mz_u), np.array(file2_mz_u)

    loess_residuals = np.abs(loess_function(file1_rt) - file2_rt)
    q = 0.01
    n = len(file1_rt)
    alpha = q * (n - np.arange(n)) / n
    p68 = np.percentile(loess_residuals, per, interpolation='linear')
    rsdr = p68 * n / (n - 2)
    rank_indices = loess_residuals.argsort()
    sorted_residuals = loess_residuals[rank_indices]
    # for i in range(int(0.7 * n) - 1, n):
    for i in range(int(0.8 * n) - 1, n):
        t_value = sorted_residuals[i] / rsdr
        p_value = 2 * (1 - t.cdf(t_value, n - 2))
        if p_value < alpha[i]:
            outlier_indices = rank_indices[i:]
            residual_cutoff = sorted_residuals[i - 1]
            break
        elif i == n - 1:
            outlier_indices = []
            residual_cutoff = sorted_residuals[i]
    mask = np.isin(np.arange(n), outlier_indices, invert=True)
    rt1 = file1_rt[mask]
    rt2 = file2_rt[mask]
    mz1 = file1_mz[mask]
    mz2 = file2_mz[mask]
    if file1_rt_u.size != 0:
        pred_file1_rt_u = loess_function(file1_rt_u)
        act_minus_pred = np.abs(file2_rt_u - pred_file1_rt_u.reshape(-1, 1))
        min_id_amp = np.argmin(act_minus_pred, axis=1)
        min_value_test = np.amin(act_minus_pred, axis=1) <= residual_cutoff
        a_rt1 = file1_rt_u[np.where(min_value_test)[0]]
        a_rt2 = file2_rt_u[np.where(min_value_test)[0], min_id_amp[np.where(min_value_test)[0]]]
        a_mz1 = file1_mz_u[np.where(min_value_test)[0]]
        a_mz2 = file2_mz_u[np.where(min_value_test)[0], min_id_amp[np.where(min_value_test)[0]]]
    else:
        a_rt1, a_rt2, a_mz1, a_mz2 = [], [], [], []

    res_rt1, res_rt2 = np.append(rt1, a_rt1), np.append(rt2, a_rt2)
    res_mz1, res_mz2 = np.append(mz1, a_mz1), np.append(mz2, a_mz2)

    # 'frac' can be adjusted to change the degree of non-linearity
    uni_res_rt1, uni_res_rt2 = preprocess_pairs(res_rt1, res_rt2)
    uni_res_rt1 = np.append(uni_res_rt1, [min_rt, max_rt])
    uni_res_rt2 = np.append(uni_res_rt2, [min_rt, max_rt])
    lowess = sm.nonparametric.lowess(uni_res_rt1, uni_res_rt2, frac=1 / 5)  # note, default frac=2/3. Mapping on file2.
    df = pd.DataFrame(lowess, columns=['rt2', 'rt1'])

    df = df.groupby('rt2', as_index=False).mean()
    f_rt = interp1d(df['rt2'], df['rt1'], bounds_error=False, fill_value='extrapolate')

    # LOESS regression on m/z
    uni_res_mz1, uni_res_mz2 = preprocess_pairs(res_mz1, res_mz2)
    lowess_mz = sm.nonparametric.lowess(uni_res_mz1, uni_res_mz2)
    df = pd.DataFrame(lowess_mz, columns=['mz2', 'mz1'])

    df = df.groupby('mz2', as_index=False).mean()
    f_mz = interp1d(df['mz2'], df['mz1'], bounds_error=False, fill_value='extrapolate')
    return f_rt, f_mz
    """uni_res_rt1, uni_res_rt2 = preprocess_pairs(res_rt1, res_rt2)
    uni_res_rt1 = np.append(uni_res_rt1, [min_rt, max_rt])
    uni_res_rt2 = np.append(uni_res_rt2, [min_rt, max_rt])

    lowess = sm.nonparametric.lowess(uni_res_rt1, uni_res_rt2, frac=1 / 5)  # note, default frac=2/3. Mapping on file2.
    df = pd.DataFrame(lowess, columns=['rt2', 'rt1'])
    df = df.groupby('rt2', as_index=False).mean()

    sorted_indices_rt = np.argsort(df['rt2'])
    uni_res_rt1 = df['rt1'][sorted_indices_rt]
    uni_res_rt2 = df['rt2'][sorted_indices_rt]

    pchip_rt = PchipInterpolator(uni_res_rt2, uni_res_rt1)

    # MZ regression
    uni_res_mz1, uni_res_mz2 = preprocess_pairs(res_mz1, res_mz2)
    lowess_mz = sm.nonparametric.lowess(uni_res_mz1, uni_res_mz2)
    df = pd.DataFrame(lowess_mz, columns=['mz2', 'mz1'])
    df = df.groupby('mz2', as_index=False).mean()

    sorted_indices_mz = np.argsort(df['mz2'])
    uni_res_mz1 = df['mz1'][sorted_indices_mz]
    uni_res_mz2 = df['mz2'][sorted_indices_mz]
    pchip_mz = PchipInterpolator(uni_res_mz2, uni_res_mz1)

    return pchip_rt, pchip_mz"""


def residuals(x, t1, t2):
    return x[0] + x[1] * t1 - t2  # This is a linear function.


def local_match_no_dtw(file1_roas, file2_eics, per=68.27):
    from scipy.optimize import least_squares
    # return f_rt, f_mz (file2 -> file1)
    # loess_function: file1 -> file2
    global outlier_indices, residual_cutoff
    file1_rt, file1_mz, file2_rt, file2_mz = [], [], [], []
    file1_rt_u, file1_mz_u, file2_rt_u, file2_mz_u = [], [], [], []
    pn = 5  # Candidate peak numbers is 5
    min_rt, max_rt = file2_eics['rt'][0], file2_eics['rt'][-1]
    between_file_fold = 10  # Control the number of high quality ROAs

    while len(file1_rt) == 0:
        for i in range(len(file1_roas)):
            if np.max(file2_eics['intensity'][i]) < file1_roas[i].i_max / between_file_fold:
                continue
            corr = signal.correlate(file2_eics['intensity'][i], file1_roas[i].i, mode='same')
            peaks, _ = signal.find_peaks(corr)
            sorted_peaks = sorted(peaks, key=lambda h: corr[h], reverse=True)
            if len(sorted_peaks) >= 2 and corr[sorted_peaks[0]] > corr[sorted_peaks[1]] * 5:
                file1_rt.append(file1_roas[i].rt[file1_roas[i].i_max_id])
                file1_mz.append(file1_roas[i].mzmean)
                file2_rt.append(file2_eics['rt'][sorted_peaks[0]])
                file2_mz.append(file2_eics['mz'][i][sorted_peaks[0]])
            elif len(sorted_peaks) == 1:
                file1_rt.append(file1_roas[i].rt[file1_roas[i].i_max_id])
                file1_mz.append(file1_roas[i].mzmean)
                file2_rt.append(file2_eics['rt'][sorted_peaks[0]])
                file2_mz.append(file2_eics['mz'][i][sorted_peaks[0]])
            elif len(sorted_peaks) >= 2:
                file1_rt_u.append(file1_roas[i].rt[file1_roas[i].i_max_id])
                file1_mz_u.append(file1_roas[i].mzmean)
                file2_rt_u.append(
                    [file2_eics['rt'][sorted_peaks[j]] if j < len(sorted_peaks) else -1000 for j in range(pn)])
                file2_mz_u.append(
                    [file2_eics['mz'][i][sorted_peaks[j]] if j < len(sorted_peaks) else -1000 for j in range(pn)])
        between_file_fold *= 2
    file1_rt, file2_rt, file1_mz, file2_mz = np.array(file1_rt), np.array(file2_rt), np.array(file1_mz), np.array(
        file2_mz)
    file1_rt_u, file2_rt_u, file1_mz_u, file2_mz_u = np.array(file1_rt_u), np.array(file2_rt_u), np.array(
        file1_mz_u), np.array(file2_mz_u)

    res_lsq = least_squares(residuals, np.ones(2), loss='soft_l1', f_scale=1.0, args=(file1_rt, file2_rt))
    loess_residuals = np.abs(residuals(res_lsq.x, file1_rt, file2_rt))
    q = 0.01
    n = len(file1_rt)
    alpha = q * (n - np.arange(n)) / n
    p68 = np.percentile(loess_residuals, per, interpolation='linear')
    rsdr = p68 * n / (n - 2)
    rank_indices = loess_residuals.argsort()
    sorted_residuals = loess_residuals[rank_indices]
    # for i in range(int(0.7 * n) - 1, n):
    for i in range(int(0.8 * n) - 1, n):
        t_value = sorted_residuals[i] / rsdr
        p_value = 2 * (1 - t.cdf(t_value, n - 2))
        if p_value < alpha[i]:
            outlier_indices = rank_indices[i:]
            residual_cutoff = sorted_residuals[i - 1]
            break
        elif i == n - 1:
            outlier_indices = []
            residual_cutoff = sorted_residuals[i]
    mask = np.isin(np.arange(n), outlier_indices, invert=True)
    rt1 = file1_rt[mask]
    rt2 = file2_rt[mask]
    mz1 = file1_mz[mask]
    mz2 = file2_mz[mask]
    if file1_rt_u.size != 0:
        pred_file1_rt_u = residuals(res_lsq.x, file1_rt_u, 0)
        act_minus_pred = np.abs(file2_rt_u - pred_file1_rt_u.reshape(-1, 1))
        min_id_amp = np.argmin(act_minus_pred, axis=1)
        min_value_test = np.amin(act_minus_pred, axis=1) <= residual_cutoff
        a_rt1 = file1_rt_u[np.where(min_value_test)[0]]
        a_rt2 = file2_rt_u[np.where(min_value_test)[0], min_id_amp[np.where(min_value_test)[0]]]
        a_mz1 = file1_mz_u[np.where(min_value_test)[0]]
        a_mz2 = file2_mz_u[np.where(min_value_test)[0], min_id_amp[np.where(min_value_test)[0]]]
    else:
        a_rt1, a_rt2, a_mz1, a_mz2 = [], [], [], []

    res_rt1, res_rt2 = np.append(rt1, a_rt1), np.append(rt2, a_rt2)
    res_mz1, res_mz2 = np.append(mz1, a_mz1), np.append(mz2, a_mz2)

    # 'frac' can be adjusted to change the degree of non-linearity
    uni_res_rt1, uni_res_rt2 = preprocess_pairs(res_rt1, res_rt2)
    uni_res_rt1 = np.append(uni_res_rt1, [min_rt, max_rt])
    uni_res_rt2 = np.append(uni_res_rt2, [min_rt, max_rt])
    lowess = sm.nonparametric.lowess(uni_res_rt1, uni_res_rt2, frac=1 / 5)  # note, default frac=2/3. Mapping on file2.
    df = pd.DataFrame(lowess, columns=['rt2', 'rt1'])

    df = df.groupby('rt2', as_index=False).mean()
    f_rt = interp1d(df['rt2'], df['rt1'], bounds_error=False, fill_value='extrapolate')

    # LOESS regression on m/z
    uni_res_mz1, uni_res_mz2 = preprocess_pairs(res_mz1, res_mz2)
    lowess_mz = sm.nonparametric.lowess(uni_res_mz1, uni_res_mz2)
    df = pd.DataFrame(lowess_mz, columns=['mz2', 'mz1'])

    df = df.groupby('mz2', as_index=False).mean()
    f_mz = interp1d(df['mz2'], df['mz1'], bounds_error=False, fill_value='extrapolate')
    return f_rt, f_mz
    """uni_res_rt1, uni_res_rt2 = preprocess_pairs(res_rt1, res_rt2)
    uni_res_rt1 = np.append(uni_res_rt1, [min_rt, max_rt])
    uni_res_rt2 = np.append(uni_res_rt2, [min_rt, max_rt])

    lowess = sm.nonparametric.lowess(uni_res_rt1, uni_res_rt2, frac=1 / 5)  # note, default frac=2/3. Mapping on file2.
    df = pd.DataFrame(lowess, columns=['rt2', 'rt1'])
    df = df.groupby('rt2', as_index=False).mean()

    sorted_indices_rt = np.argsort(df['rt2'])
    uni_res_rt1 = df['rt1'][sorted_indices_rt]
    uni_res_rt2 = df['rt2'][sorted_indices_rt]

    pchip_rt = PchipInterpolator(uni_res_rt2, uni_res_rt1)

    # MZ regression
    uni_res_mz1, uni_res_mz2 = preprocess_pairs(res_mz1, res_mz2)
    lowess_mz = sm.nonparametric.lowess(uni_res_mz1, uni_res_mz2)
    df = pd.DataFrame(lowess_mz, columns=['mz2', 'mz1'])
    df = df.groupby('mz2', as_index=False).mean()

    sorted_indices_mz = np.argsort(df['mz2'])
    uni_res_mz1 = df['mz1'][sorted_indices_mz]
    uni_res_mz2 = df['mz2'][sorted_indices_mz]
    pchip_mz = PchipInterpolator(uni_res_mz2, uni_res_mz1)

    return pchip_rt, pchip_mz"""


if __name__ == '__main__':
    pass
