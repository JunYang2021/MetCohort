from BTrees.OOBTree import OOBTree
from file_processing import MsFile
import numpy as np


class ROA:
    def __init__(self, scan, rt, i, mz, mzmean):
        self.scan = scan
        self.rt = rt
        self.i = i
        self.mz = mz
        self.mzmean = mzmean
        self.i_max = i[0]  # max(i)
        self.i_max_id = 0  # i.index(max(i))
        self.points = 1

    def __repr__(self):
        return 'mz = {:.4f}, rt = {:.2f} - {:.2f}, rt_max = {:.2f}'.format(self.mzmean,
                                                                           self.rt[0],
                                                                           self.rt[-1],
                                                                           self.rt[self.i_max_id])

    def plot(self):
        import matplotlib.pyplot as plt

        font = {'family': 'Arial', 'size': 16}
        plt.rc('font', **font)
        fig, ax = plt.subplots(figsize=(7, 2.5))
        ax.plot(self.rt, self.i, linewidth=2)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel('Retention time (s)')
        ax.set_ylabel('Intensity')
        # ax.set_title(self.__repr__())
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', which='both', left=True, right=False, labelleft=False)
        ax.locator_params(axis='x', nbins=4)
        # ax.locator_params(axis='y', nbins=3)
        plt.tight_layout()
        # plt.savefig('roa{0}.tiff'.format(self.mzmean), dpi=300)
        plt.show()


def create_new_roa(rois_tree, scan_number, time, i, mz, scans, window_width, delta_mz, time_array):
    rois_tree[mz] = ROA([scan_number, scan_number],
                        [time],
                        [i],
                        [mz],
                        mz)
    isbreak = False
    for id_scan in range(scan_number - 2, -1, -1):
        time_now = time_array[id_scan]
        if time - time_now < window_width / 2:
            mask = (scans[id_scan].mz > (rois_tree[mz].mzmean - delta_mz)) & (scans[id_scan].mz < (rois_tree[mz].mzmean + delta_mz))
            notzero_mz = scans[id_scan].mz[mask]
            notzero_i = scans[id_scan].i[mask]

            if np.any(notzero_i > i):
                isbreak = True
                break
            rois_tree[mz].scan[0] = id_scan + 1
            rois_tree[mz].rt.insert(0, time_now)
            if notzero_mz.size == 0:
                rois_tree[mz].i.insert(0, 0)
                rois_tree[mz].mz.insert(0, mz)
            else:
                rois_tree[mz].i.insert(0, np.sum(notzero_i))
                rois_tree[mz].mz.insert(0, np.sum(notzero_mz * notzero_i) / np.sum(notzero_i))
                rois_tree[mz].mzmean = (rois_tree[mz].mzmean * rois_tree[mz].points + np.sum(notzero_mz)) / (
                        rois_tree[mz].points + notzero_mz.size)
                rois_tree[mz].points += notzero_mz.size
            rois_tree[mz].i_max_id = rois_tree[mz].i.index(i)
        else:
            break
    if max(rois_tree[mz].i) > i or isbreak:
        # rois_tree.remove_items([mz])
        del rois_tree[mz]


def check_peak(i_array):
    from scipy import signal
    peaks, _ = signal.find_peaks(i_array)
    if not peaks.size:
        return False
    sorted_peaks = sorted(peaks, key=lambda h: i_array[h], reverse=True)
    if len(sorted_peaks) == 1 or (len(sorted_peaks) >= 2 and i_array[sorted_peaks[0]] > i_array[sorted_peaks[1]] * 5):
        return True
    return False


def roa_construction(msfile: MsFile, delta_mz=0.001, int_cutoff=0.01, window_width=30, progress_callback=None, alpha=0.1):
    """

    :param msfile: MSFile
    :param delta_mz:
    :param int_cutoff: Max_intensity cutoff
    :param window_width: unit is seconds
    :return:
    """
    # The only decision a user of the EWMA must make is the parameter alpha.
    # The parameter decides how important the current observation is in the calculation of the EWMA. The higher the
    # value of alpha, the more closely the EWMA tracks the original time series.
    time_max = msfile.original_time[-1]
    first = True
    ewma = []
    for scan in msfile.exp:
        if first:
            ewma.append(scan.i.sum())
            first = False
        elif not first:
            ewma.append(alpha * scan.i.sum() + (1 - alpha) * ewma[-1])
    roa_list = []
    process_rois = OOBTree()

    number = 1
    for scan in msfile.exp:
        if number == 1:
            number += 1
            continue
        ctime = msfile.original_time[number - 1]
        for mz, i in zip(*(scan.mz, scan.i)):
            if i > ewma[number - 1] * int_cutoff:  # add points or create new roi
                closest_mz, closest_item = None, None
                for candidate_mz, candidate_item in process_rois.items(min=mz - delta_mz, max=mz + delta_mz):
                    if closest_mz is None or abs(candidate_mz - mz) < abs(closest_mz - mz):
                        closest_mz, closest_item = candidate_mz, candidate_item

                if closest_item is None:
                    create_new_roa(process_rois, number, ctime, i, mz, msfile.exp, window_width, delta_mz, msfile.original_time)
                elif abs(closest_item.mzmean - mz) < delta_mz:
                    # add point to closest_item
                    if closest_item.scan[1] == number:
                        closest_item.mz[-1] = (closest_item.i[-1] * closest_item.mz[-1] + i * mz) / (
                                closest_item.i[-1] + i)
                        closest_item.i[-1] = closest_item.i[-1] + i
                        closest_item.mzmean = (closest_item.mzmean * closest_item.points + mz) / (
                                closest_item.points + 1)
                        closest_item.points += 1
                        closest_item.i_max = max(closest_item.i)
                        closest_item.i_max_id = closest_item.i.index(closest_item.i_max)
                    else:
                        closest_item.scan[1] = number
                        closest_item.rt.append(ctime)
                        closest_item.i.append(i)
                        closest_item.mz.append(mz)
                        closest_item.mzmean = (closest_item.mzmean * closest_item.points + mz) / (
                                closest_item.points + 1)
                        closest_item.points += 1
                        closest_item.i_max = max(closest_item.i)
                        closest_item.i_max_id = closest_item.i.index(closest_item.i_max)
                else:
                    # create new roi
                    create_new_roa(process_rois, number, ctime, i, mz, msfile.exp, window_width, delta_mz, msfile.original_time)
        # has finished one ctime step, now should ensure it is a good roi
        # crop the head; if tail is long enough, save it.
        to_delete = []
        for mz, roi in process_rois.items():  # 没有增加点的roi需要加一个点
            if roi.scan[1] != number:
                roi.scan[1] = number
                roi.rt.append(ctime)
                mask = (scan.mz > (mz - delta_mz)) & (scan.mz < (mz + delta_mz))
                notzero_mz = scan.mz[mask]
                notzero_i = scan.i[mask]

                if notzero_mz.size == 0:
                    roi.i.append(0)
                    roi.mz.append(roi.mzmean)
                else:
                    roi.i.append(np.sum(notzero_i))
                    roi.mz.append(np.sum(notzero_mz * notzero_i) / np.sum(notzero_i))
                    roi.mzmean = (roi.mzmean * roi.points + np.sum(notzero_mz)) / (
                            roi.points + notzero_mz.size)
                    roi.points += notzero_mz.size
            if roi.scan[0] != number:  # If is added in this ctime step, no need to crop or save
                while roi.rt[roi.i_max_id] - roi.rt[0] >= window_width / 2:
                    roi.scan[0] += 1
                    if roi.i[0] != 0:
                        roi.mzmean = (roi.mzmean * roi.points - roi.mz[0]) / (roi.points - 1)
                        roi.points -= 1
                    del roi.rt[0]
                    del roi.i[0]
                    del roi.mz[0]
                    roi.i_max_id -= 1
                if roi.rt[-1] - roi.rt[roi.i_max_id] >= window_width / 2 or roi.rt[-1] == time_max:
                    roa_list.append(roi)
                    to_delete.append(mz)
        for item in to_delete:
            del process_rois[item]
        number += 1
        if progress_callback is not None:
            progress_callback.emit(int(number * 100 / len(msfile.exp)))
    roa_list = [roa for roa in roa_list if roa.i_max_id != 0 and roa.i_max_id != len(roa.i) - 1 and check_peak(roa.i)]
    return roa_list


def get_eic_array(msfile: MsFile, mz_list, delta_mz=0.01):
    mz_list = np.array(mz_list)

    rt_array = msfile.original_time
    mz_array = np.zeros((mz_list.shape[0], rt_array.shape[0]), dtype=np.float32)
    i_array = np.zeros((mz_list.shape[0], rt_array.shape[0]), dtype=np.float32)

    d = 0
    for scan in msfile.exp:
        scan_mz, scan_i = np.array(scan.mz), np.array(scan.i)
        L = len(scan_mz)
        pos = np.searchsorted(scan_mz, mz_list)  # shape: (len(mz_list),)
        pos = np.clip(np.vstack((pos, pos - 1)), 0, L - 1)  # shape: (2, len(mz_list))
        closest = np.vstack((scan_mz[pos[0]], scan_mz[pos[1]]))  # shape: (2, len(mz_list))
        distance = np.abs(closest - mz_list)  # shape: (2, len(mz_list))
        min_distance = np.amin(distance, axis=0)  # shape: (len(mz_list), )
        min_distance_ = np.argmin(distance, axis=0)  # shape: (len(mz_list), )
        mask = min_distance <= delta_mz
        pos = pos[min_distance_, range(mz_list.shape[0])]
        mz_array[:, d] = np.where(mask, scan_mz[pos], mz_list)
        i_array[:, d] = np.where(mask, scan_i[pos], 0)
        d += 1

    return {'rt': rt_array, 'mz': mz_array, 'intensity': i_array}


def find_roi_range(msfile, file_list, delta_mz=0.005, continuous_points=3,
                   dropped_points=10):
    """
        Get all roi ranges from a mzml file
        :param msfile:
        :param data: A list include dictionaries to store ranges [{'file': id, 'mz': mzmean, 'to':start time, 't1': end time}]
        :return:
    """
    data = []
    id = file_list.index(msfile.file_path)
    process_ranges = OOBTree()

    for number, scan in enumerate(msfile.exp):
        ctime = msfile.corrected_time[number]
        for mz, i in zip(scan.mz, scan.i):
            if i != 0:
                closest_mz, closest_range = None, None
                for candidate_mz, candidate_range in process_ranges.items(min=mz - delta_mz, max=mz + delta_mz):
                    if closest_mz is None or abs(candidate_mz - mz) < abs(closest_mz - mz):
                        # One point can only be added to one roi
                        closest_mz, closest_range = candidate_mz, candidate_range
                if closest_range is None:
                    # rp: real points; si: last scan index of real points; cp: continuous points
                    process_ranges[mz] = {'mz': mz, 't0': ctime, 't1': ctime, 'rp': 1, 'si': number, 'cp': [1]}
                else:
                    del process_ranges[closest_mz]
                    closest_range['mz'] = (closest_range['mz'] * closest_range['rp'] + mz) / (closest_range['rp'] + 1)
                    closest_range['t1'] = ctime
                    closest_range['rp'] += 1
                    if number - closest_range['si'] >= 2:
                        closest_range['cp'].append(1)
                    elif number - closest_range['si'] == 1:
                        closest_range['cp'][-1] += 1
                    closest_range['si'] = number
                    process_ranges[closest_range['mz']] = closest_range
        to_delete = []
        for mz, t_range in process_ranges.items():
            if number - t_range['si'] == dropped_points:
                to_delete.append(mz)
                if max(t_range['cp']) >= continuous_points:
                    data.append({'file': id, 'mz': mz, 't0': t_range['t0'], 't1': t_range['t1']})
        for i in to_delete:
            del process_ranges[i]

    for mz, t_range in process_ranges.items():
        if max(t_range['cp']) >= continuous_points:
            data.append({'file': id, 'mz': mz, 't0': t_range['t0'], 't1': t_range['t1']})

    return data


if __name__ == '__main__':
    import time
    time0 = time.time()
    msfile = MsFile('E:\\batch2-pos data\\S0002.mzML')
    time1 = time.time()
    roas = roa_construction(msfile)
    # for scan in msfile.exp:
    #     for mz, i in zip(*(scan.mz, scan.i)):
    #         pass
    time2 = time.time()
    print(len(roas))
    print(time1 - time0, time2 - time1)
