from pyopenms import *
import numpy as np


def spectrum_to_vector(mz, intensity, min_mz=100, max_mz=1000, bin_width=0.5):
    """
        Convert a mass spectrum to a fixed-length vector using binning with np.histogram,
        including filtering out low-intensity values and normalizing the vector.

        Parameters:
        - mz: array of m/z values
        - intensity: array of intensity values
        - min_mz, max_mz: range of m/z values to consider
        - bin_width: width of each bin in m/z units

        Returns:
        - A normalized, fixed-length vector of aggregated intensities.
        """
    # 1. Filter out low-intensity m/z values
    max_intensity = np.max(intensity)
    threshold = 0.001 * max_intensity  # 0.1% of the maximum intensity
    mask = intensity > threshold
    mz_filtered = mz[mask]
    intensity_filtered = intensity[mask]

    # 2. Define bins
    num_bins = int(np.ceil((max_mz - min_mz) / bin_width))
    bins = np.linspace(min_mz, max_mz, num_bins + 1)

    # 3. Use np.histogram to aggregate intensities into bins
    vector, _ = np.histogram(mz_filtered, bins=bins, weights=intensity_filtered, density=False)

    # 4. Normalize the vector to make the sum of intensities equal to 1
    vector_sum = np.sum(vector)
    if vector_sum > 0:  # Prevent division by zero
        vector = vector / vector_sum

    return vector


class ms1scan:
    # MS1 scan
    def __init__(self, mz, i):
        self.mz = mz.astype(np.float32)
        self.i = i.astype(np.float32)


class MsFile:
    def __init__(self, filepath, mzrange=(100, 1000), rt_filter=None):
        self.file_path = filepath
        reader = MSExperiment()
        if filepath.endswith('.mzML'):
            MzMLFile().load(filepath, reader)
        elif filepath.endswith('.mzXML'):
            MzXMLFile().load(filepath, reader)
        self.exp = []
        self.original_time = []
        self.bin_spectrum = []
        if rt_filter is None:
            for scan in reader:
                if scan.getMSLevel() == 1:
                    self.original_time.append(scan.getRT())
                    _scan = ms1scan(scan.get_peaks()[0], scan.get_peaks()[1])
                    self.bin_spectrum.append(spectrum_to_vector(_scan.mz, _scan.i, min_mz=mzrange[0], max_mz=mzrange[1]))
                    self.exp.append(_scan)
        else:
            for scan in reader:
                if scan.getMSLevel() == 1:
                    ret_t = scan.getRT()
                    if rt_filter[0] <= ret_t <= rt_filter[1]:
                        _scan = ms1scan(scan.get_peaks()[0], scan.get_peaks()[1])
                        self.original_time.append(ret_t)
                        self.bin_spectrum.append(spectrum_to_vector(_scan.mz, _scan.i))
                        self.exp.append(_scan)
                    elif ret_t > rt_filter[1]:
                        break
        self.original_time = np.array(self.original_time, dtype=np.float32)
        self.corrected_time = self.original_time.copy()
        print(filepath, 'load.')

    def map_file(self, f_rt, f_mz):
        new_exp = []
        for j, scan in enumerate(self.exp):
            old_time = self.original_time[j]
            new_time = f_rt(old_time).item()
            self.corrected_time[j] = new_time
            new_mz = f_mz(scan.mz)
            sorted_indices = np.argsort(new_mz)
            new_exp.append(ms1scan(new_mz[sorted_indices], scan.i[sorted_indices]))
        self.exp = new_exp
        print(self.file_path, 'aligned.')
        # Release memory
        # self.bin_spectrum = None

    def save_file(self, new_path):
        tmp_exp = MSExperiment()
        for j, scan in enumerate(self.exp):
            tmp_spectrum = MSSpectrum()
            tmp_spectrum.setRT(self.corrected_time[j])
            tmp_spectrum.setMSLevel(1)
            tmp_spectrum.set_peaks((scan.mz, scan.i))
            tmp_spectrum.sortByPosition()
            tmp_exp.addSpectrum(tmp_spectrum)
        tmp_exp.sortSpectra(True)

        if new_path.endswith('.mzML'):
            MzMLFile().store(new_path, tmp_exp)
        elif new_path.endswith('.mzXML'):
            MzXMLFile().load(new_path, tmp_exp)
        print(new_path, 'stored.')


if __name__ == '__main__':
    pass
