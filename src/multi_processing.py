from local_match import local_match, local_match_no_dtw
from file_processing import MsFile
from dtw import dtw, loess_dtw
from get_chromatograms import get_eic_array


def load_data_parallel(file_path, min_mz, max_mz, rt_filter=None):
    return MsFile(file_path, mzrange=(min_mz, max_mz), rt_filter=rt_filter)


def align_data_parallel(msfile, msfile_ref, roas, roa_mz_list, delta_mz_eic, no_dtw=True):
    if msfile.file_path != msfile_ref.file_path:
        if no_dtw:
            file2_eic_array = get_eic_array(msfile, roa_mz_list, delta_mz_eic)
            f_rt, f_mz = local_match_no_dtw(roas, file2_eic_array)
            msfile.map_file(f_rt, f_mz)
            return msfile
        else:
            rt1, rt2 = dtw(msfile_ref.original_time,
                           msfile_ref.bin_spectrum,
                           msfile.original_time,
                           msfile.bin_spectrum)
            loess_func_1t2 = loess_dtw(rt2, rt1)

            file2_eic_array = get_eic_array(msfile, roa_mz_list, delta_mz_eic)
            f_rt, f_mz = local_match(roas, file2_eic_array, loess_func_1t2)
            msfile.map_file(f_rt, f_mz)
            return msfile
    else:
        return msfile