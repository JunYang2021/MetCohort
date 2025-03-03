import time
from pyopenms import *
import multiprocessing as mp
from multi_processing import load_data_parallel, align_data_parallel
from get_chromatograms import roa_construction, find_roi_range, get_eic_array
from dtw import dtw, loess_dtw
from file_processing import MsFile
from local_match import local_match, local_match_no_dtw


def load_files(file_list, rt_filter=None):
    time0 = time.time()
    length = len(file_list)
    # Read a temp file and get the m/z range
    temp_p = MSExperiment()
    if file_list[0].endswith('.mzML'):
        MzMLFile().load(file_list[0], temp_p)
    elif file_list[0].endswith('.mzXML'):
        MzXMLFile().load(file_list[0], temp_p)
    min_mz, max_mz = np.inf, -1
    for scan in temp_p:
        if scan.get_peaks()[0].size:
            _min, _max = scan.get_peaks()[0].min(), scan.get_peaks()[0].max()
            if min_mz > _min:
                min_mz = _min
            if max_mz < _max:
                max_mz = _max

    print("Loading {0} files...".format(length))
    exps = {}
    for file_path in file_list:
        exps[file_path] = MsFile(file_path, mzrange=(min_mz, max_mz), rt_filter=rt_filter)
    print("Data imported successfully.")
    time1 = time.time()
    print('Time spending: {0} seconds.'.format(time1 - time0))
    return exps


def correct_files(ref_file, roa_window_width, delta_mz_roa, intensity_threshold_roa,
                 delta_mz_eic, plot, write_new_files, new_file_directory, exps, file_list,
                 blank_file_list, baseline='rl'):
    # baseline = 'rl'(Robust linear fitting) or 'DTW'
    time0 = time.time()
    print('Start alignment...')
    print('ROA construction from reference file {0}...'.format(ref_file))
    roas = roa_construction(exps[ref_file], delta_mz_roa,
                            intensity_threshold_roa, roa_window_width)
    print('{0} ROAs found.'.format(len(roas)))
    if len(roas) < 100:
        print('Warning: Less than 100 ROAs were found. Correction may be inaccurate.')
    mz_list = [i.mzmean for i in roas]

    print(
        'Retention time alignment to Reference file...')

    for file_path in file_list:
        if file_path != ref_file and file_path not in blank_file_list:
            if baseline == 'rl':
                file2_eic_array = get_eic_array(exps[file_path], mz_list, delta_mz_eic)
                f_rt = local_match_no_dtw(roas, file2_eic_array)
                exps[file_path].map_file(f_rt)
            elif baseline == 'DTW':
                rt1, rt2 = dtw(exps[ref_file].original_time,
                               exps[ref_file].bin_spectrum,
                               exps[file_path].original_time,
                               exps[file_path].bin_spectrum)
                loess_func_1t2 = loess_dtw(rt2, rt1)

                file2_eic_array = get_eic_array(exps[file_path], mz_list, delta_mz_eic)
                f_rt = local_match(roas, file2_eic_array, loess_func_1t2)
                exps[file_path].map_file(f_rt)

    print('Alignment successfully.')
    time1 = time.time()
    print('Time spending: {0} seconds.'.format(time1 - time0))

    if plot:
        from bokeh.plotting import figure, output_file, save
        from bokeh.models import HoverTool, ColumnDataSource
        from bokeh.palettes import viridis
        import os, datetime
        import random
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
        pic_path = os.path.join(new_file_directory, f'Alignment results{timestamp}.html')
        output_file(pic_path)

        p = figure(width=1200, height=675, title="Alignment Results")
        colors = viridis(256)
        print('Making plot...')
        for i, path in enumerate(file_list):
            if path != ref_file:
                y_values = exps[path].original_time - exps[path].corrected_time
                source = ColumnDataSource(data={'x': exps[path].corrected_time,
                                                'y': y_values,
                                                'name': [path] * len(y_values)})
                p.line('x', 'y', color=random.choice(colors), source=source)

        hover = HoverTool(tooltips=[('File', '@name')])
        p.add_tools(hover)
        p.xaxis.axis_label = 'Retention Time (s) of Reference sample ({0})'.format(
            os.path.basename(ref_file))
        p.yaxis.axis_label = 'Deviation of Retention Time (s)'
        save(p)
        print('Plot saved successfully.')
    if write_new_files:
        import os
        input_filename = [os.path.basename(i) for i in file_list]
        new_path = [os.path.join(new_file_directory, i) for i in input_filename]
        print('Writing new files to {0}...'.format(new_file_directory))
        for i, path in enumerate(file_list):
            exps[path].save_file(new_path[i])
        print('Writing new files successfully.')
    return exps
    

# Only for test
def correct_files_only_dtw(ref_file, plot, write_new_files, new_file_directory, exps, file_list):
    for file_path in file_list:
        if file_path != ref_file:
            rt1, rt2 = dtw(exps[ref_file].original_time,
                           exps[ref_file].bin_spectrum,
                           exps[file_path].original_time,
                           exps[file_path].bin_spectrum)
            loess_func_2t1 = loess_dtw(rt1, rt2)
            exps[file_path].map_file(loess_func_2t1)
    if plot:
        from bokeh.plotting import figure, output_file, save
        from bokeh.models import HoverTool, ColumnDataSource
        from bokeh.palettes import viridis
        import os, datetime
        import random
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
        pic_path = os.path.join(new_file_directory, f'Alignment results{timestamp}.html')
        output_file(pic_path)

        p = figure(width=1200, height=675, title="Alignment Results")
        colors = viridis(256)
        print('Making plot...')
        for i, path in enumerate(file_list):
            if path != ref_file:
                y_values = exps[path].original_time - exps[path].corrected_time
                source = ColumnDataSource(data={'x': exps[path].corrected_time,
                                                'y': y_values,
                                                'name': [path] * len(y_values)})
                p.line('x', 'y', color=random.choice(colors), source=source)

        hover = HoverTool(tooltips=[('File', '@name')])
        p.add_tools(hover)
        p.xaxis.axis_label = 'Retention Time (s) of Reference sample ({0})'.format(
            os.path.basename(ref_file))
        p.yaxis.axis_label = 'Deviation of Retention Time (s)'
        save(p)
        print('Plot saved successfully.')
    if write_new_files:
        import os
        input_filename = [os.path.basename(i) for i in file_list]
        new_path = [os.path.join(new_file_directory, i) for i in input_filename]
        print('Writing new files to {0}...'.format(new_file_directory))
        for i, path in enumerate(file_list):
            exps[path].save_file(new_path[i])
        print('Writing new files successfully.')
    return exps
    


