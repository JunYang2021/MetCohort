# This version use DTW for profile alignment background.
# This version use pyopenms for LC-MS raw data reading and writing
import itertools
import multiprocessing as mp
import os
from pyopenms import *
import pandas as pd
from PyQt5 import QtWidgets, uic, QtCore
from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QTableWidgetItem, QCheckBox
from PyQt5.QtGui import QIcon
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
# from memory_profiler import profile
from get_chromatograms import roa_construction, find_roi_range
from batch_roi import grouping_roi_range, peak_detection_roi_matrix
from multi_processing import load_data_parallel, align_data_parallel

mp.freeze_support()


def sort_key(x):
    return x['mz'], x['rt'], x['tar_label']


def unique_key(x):
    return x['mz'], x['rt']


def unique_and_sorted(peaks):  # remove replicate peaks in results
    seen = set()
    unique_peaks = []

    for p in sorted(peaks, key=sort_key):
        key = unique_key(p)
        if key not in seen:
            seen.add(key)
            unique_peaks.append(p)

    return unique_peaks


class FileLoader(QThread):
    update_signal = pyqtSignal(str)
    finished_signal = pyqtSignal()

    def __init__(self, file_list, exps_dict, rt_filter):
        super(FileLoader, self).__init__()
        self.file_list = file_list
        self.exps = exps_dict
        self.rt_filter = rt_filter

    def run(self):
        try:
            import time
            time0 = time.time()
            length = len(self.file_list)
            # Read a temp file and get the m/z range
            temp_p = MSExperiment()
            if self.file_list[0].endswith('.mzML'):
                MzMLFile().load(self.file_list[0], temp_p)
            elif self.file_list[0].endswith('.mzXML'):
                MzXMLFile().load(self.file_list[0], temp_p)
            min_mz, max_mz = np.inf, -1
            for scan in temp_p:
                _min, _max = scan.get_peaks()[0].min(), scan.get_peaks()[0].max()
                if min_mz > _min:
                    min_mz = _min
                if max_mz < _max:
                    max_mz = _max

            self.update_signal.emit("Loading {0} files (Parallel)... (This step may be time consuming)".format(length))

            with mp.Pool(max(1, mp.cpu_count() // 4), maxtasksperchild=10) as pool:
                results_list = pool.starmap(load_data_parallel, [(file_path, min_mz, max_mz, self.rt_filter)
                                                                 for file_path in self.file_list])

            for file_path, result in zip(self.file_list, results_list):
                self.exps[file_path] = result

            self.update_signal.emit("Data imported successfully.")
            time1 = time.time()
            self.update_signal.emit('Time spending: {0} seconds.'.format(time1 - time0))
            self.finished_signal.emit()
        except Exception as e:
            import traceback
            print(f"Caught an exception: {e}")
            traceback.print_exc()


class CorrectionThread(QThread):
    update_signal = pyqtSignal(str)
    progress_signal = pyqtSignal(int)
    finished_signal = pyqtSignal()
    hide_bar_signal = pyqtSignal()
    activate_progress_signal = pyqtSignal()

    def __init__(self, ref_file, roa_window_width, delta_mz_roa, intensity_threshold_roa,
                 delta_mz_eic, plot, write_new_files, new_file_directory, exps, file_list):
        super(CorrectionThread, self).__init__()
        self.ref_file = ref_file
        self.roa_window_width = roa_window_width
        self.delta_mz_roa = delta_mz_roa
        self.intensity_threshold_roa = intensity_threshold_roa
        self.delta_mz_eic = delta_mz_eic
        self.plot = plot
        self.write_new_files = write_new_files
        self.new_file_directory = new_file_directory
        self.exps = exps
        self.file_list = file_list

    def run(self):
        global roas
        try:
            import time, os
            time0 = time.time()
            L = len(self.file_list)
            self.update_signal.emit('Start alignment...')
            self.update_signal.emit('ROA construction from reference file {0}...'.format(self.ref_file))
            roas = roa_construction(self.exps[self.ref_file], self.delta_mz_roa,
                                    self.intensity_threshold_roa, self.roa_window_width, self.progress_signal)
            self.progress_signal.emit(100)
            self.update_signal.emit('{0} ROAs found.'.format(len(roas)))
            if len(roas) < 100:
                self.update_signal.emit('Warning: Less than 100 ROAs were found. Correction may be inaccurate.')
            mz_list = [i.mzmean for i in roas]
            self.hide_bar_signal.emit()

            self.update_signal.emit('Retention time alignment to Reference file (Parallel)... (This step may be time consuming)')

            with mp.Pool(max(1, mp.cpu_count() // 4), maxtasksperchild=10) as pool:
                aligned_res = pool.starmap(align_data_parallel, [(self.exps[file_path], self.exps[self.ref_file],
                                                                  roas, mz_list, self.delta_mz_eic)
                                                                 for file_path in self.file_list])
            for file_path, result in zip(self.file_list, aligned_res):
                self.exps[file_path] = result

            self.update_signal.emit('Alignment successfully.')
            time1 = time.time()
            self.update_signal.emit('Time spending: {0} seconds.'.format(time1 - time0))

            if self.plot:
                from bokeh.plotting import figure, output_file, save
                from bokeh.models import HoverTool, ColumnDataSource
                from bokeh.palettes import viridis
                import os, datetime
                import random
                timestamp = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M")
                pic_path = os.path.join(self.new_file_directory, f'Alignment results{timestamp}.html')
                output_file(pic_path)

                p = figure(width=1200, height=675, title="Alignment Results")
                colors = viridis(256)
                self.update_signal.emit('Making plot...')
                self.activate_progress_signal.emit()
                for i, path in enumerate(self.file_list):
                    if path != self.ref_file:
                        y_values = self.exps[path].original_time - self.exps[path].corrected_time
                        source = ColumnDataSource(data={'x': self.exps[path].corrected_time,
                                                        'y': y_values,
                                                        'name': [path] * len(y_values)})
                        p.line('x', 'y', color=random.choice(colors), source=source)
                    self.progress_signal.emit(int((i + 1) * 100 / L))
                self.progress_signal.emit(100)

                hover = HoverTool(tooltips=[('File', '@name')])
                p.add_tools(hover)
                p.xaxis.axis_label = 'Retention Time (s) of Reference sample ({0})'.format(
                    os.path.basename(self.ref_file))
                p.yaxis.axis_label = 'Deviation of Retention Time (s)'
                save(p)
                self.update_signal.emit('Plot saved successfully.')
            if self.write_new_files:
                import os
                input_filename = [os.path.basename(i) for i in self.file_list]
                new_path = [os.path.join(self.new_file_directory, i) for i in input_filename]
                self.update_signal.emit('Writing new files to {0}...'.format(self.new_file_directory))
                self.activate_progress_signal.emit()
                for i, path in enumerate(self.file_list):
                    self.exps[path].save_file(new_path[i])
                    self.progress_signal.emit(int((i + 1) * 100 / L))
                self.progress_signal.emit(100)
                self.update_signal.emit('Writing new files successfully.')

            self.finished_signal.emit()
        except Exception as e:
            import traceback
            print(f"Caught an exception: {e}")
            traceback.print_exc()
            self.update_signal.emit(f"Caught an exception: {e}")
            if len(roas) == 0:
                self.update_signal.emit("Reason: Intensity threshold is too high to find ROA. "
                                        "Please reduce the value in 'Intensity threshold of ROA'.")
            if isinstance(e, IndexError) and str(e) == 'index -1 is out of bounds for axis 0 with size 0':
                print("Possible reason: There is one empty scan in current file. Please check the file conversion "
                      "process, or remove this file and try again.")


class PeakDetectionThread(QThread):
    update_signal = pyqtSignal(str)
    progress_signal = pyqtSignal(int)
    result_signal = pyqtSignal(object)
    finished_signal = pyqtSignal()
    activate_progress_signal = pyqtSignal()

    def __init__(self, min_files_zero, delta_mz,
                 continuous_points, dropped_points, exps,
                 file_list, maximum_peak_width, minimum_intensity, ref_file, qc_file_list, lamda, tar_cpd):
        # required_points改为continuous_time, dropped_points改为dropped_time
        super(PeakDetectionThread, self).__init__()
        self.delta_mz = delta_mz
        self.continuous_points = continuous_points
        self.dropped_points = dropped_points
        self.exps = exps
        self.file_list = file_list
        self.maximum_peak_width = maximum_peak_width
        self.maximumzerofiles = int(len(file_list) * (1 - min_files_zero))
        self.min_intensity = minimum_intensity
        self.ref_file = ref_file
        self.qc_file_list = qc_file_list
        self.lamda = lamda
        self.tar_cpd = tar_cpd  # A list of dictionaries

    def run(self):
        import time
        time0 = time.time()
        try:
            # Release memory
            for file in self.file_list:
                self.exps[file].bin_spectrum = None

            self.update_signal.emit('ROI detecting in {0} QC files (Parallel)...(This step may be time consuming)'.format(len(self.qc_file_list)))

            with mp.Pool(max(1, mp.cpu_count() // 4), maxtasksperchild=10) as pool:
                range_data = pool.starmap(find_roi_range, [(self.exps[file_path], self.file_list,
                                                            self.delta_mz, self.continuous_points, self.dropped_points)
                                                           for file_path in self.qc_file_list])
            # self.progress_signal.emit(100)
            range_data = list(itertools.chain.from_iterable(range_data))

            self.update_signal.emit('ROI grouping...(This step may be time consuming)')
            range_data = pd.DataFrame(range_data).astype({'file': 'int16', 'mz': 'float32',
                                                          't0': 'float32', 't1': 'float32'})
            range_data = grouping_roi_range(range_data,
                                            len(self.qc_file_list),
                                            self.delta_mz)
            self.update_signal.emit('QC ROI detection and grouping finished. '
                                    '{0} batchROIs found.'.format(len(range_data)))
            self.update_signal.emit('Untargeted feature detection on ROI matrix...')
            self.activate_progress_signal.emit()
            range_data = [tuple(x) for x in range_data.to_numpy()]
            L= len(range_data)
            peaks_res = []
            # self.progress_signal.emit(0)
            for _t, range_row in enumerate(range_data):
                peaks_res.append(peak_detection_roi_matrix(self.exps, range_row, self.file_list,
                                                           self.ref_file, self.delta_mz,
                                                           self.maximum_peak_width,
                                                           self.maximumzerofiles, self.min_intensity,
                                                           self.lamda))
                self.progress_signal.emit(int((_t + 1) * 100 / L))
            self.progress_signal.emit(100)
            if self.tar_cpd:
                self.update_signal.emit('Targeted feature detection of {0} compounds...'.format(len(self.tar_cpd)))
                self.progress_signal.emit(0)
                for _t, cur_tar_cpd in enumerate(self.tar_cpd):
                    peaks_res.append(peak_detection_roi_matrix(self.exps, (cur_tar_cpd['m/z'], cur_tar_cpd['RT(s)'] - 20, cur_tar_cpd['RT(s)'] + 20), self.file_list,
                                                               self.ref_file, self.delta_mz,
                                                               self.maximum_peak_width,
                                                               self.maximumzerofiles, self.min_intensity,
                                                               self.lamda, cur_tar_cpd))
                    self.progress_signal.emit(int((_t + 1) * 100 / L))
                self.progress_signal.emit(100)

            peaks_res = [j for peak_res in peaks_res for j in peak_res]
            print('Organize peak table.')
            peaks_res = unique_and_sorted(peaks_res)
            print('Unique and sort peak table.')
            self.update_signal.emit('Peak detection finished. {0} peaks found.'.format(len(peaks_res)))
            self.result_signal.emit(peaks_res)
            print('Emit peak table.')
            self.finished_signal.emit()
        except Exception as e:
            import traceback
            print(f"Caught an exception: {e}")
            traceback.print_exc()
        time1 = time.time()
        self.update_signal.emit('Time spending: {0} seconds.'.format(time1 - time0))


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        uic.loadUi("main_window.ui", self)
        self.setWindowIcon(QIcon('favicon.ico'))

        self.setMinimumSize(1600, 1200)
        self.progressBar.hide()

        self.button_list = [self.pushButton, self.pushButton_2, self.pushButton_3, self.pushButton_4]
        for i, button in enumerate(self.button_list):
            button.clicked.connect(lambda checked, i=i: self.change_page(i))

        self.pushButton.setCheckable(True)
        self.pushButton_2.setCheckable(True)
        self.pushButton_3.setCheckable(True)
        self.pushButton_4.setCheckable(True)

        # List to store file paths
        self.file_list = []
        self.exps = {}
        self.peaks_res = None
        self.file_loader = None
        self.correction_thread = None
        self.peak_detection_thread = None

        self.qc_file_list = []

        # Connect Browse and Start Importing buttons
        self.pushButton_5.clicked.connect(self.browse_files)
        self.pushButton_6.clicked.connect(self.start_importing)

        self.comboBox.addItems(self.file_list)
        self.comboBox_4.addItems(self.file_list)

        self.lineEdit_6.setDisabled(True)
        self.pushButton_8.setDisabled(True)
        self.comboBox_2.currentIndexChanged.connect(self.update_file_directory_setting)

        self.pushButton_8.clicked.connect(self.browse_directory)
        self.pushButton_7.clicked.connect(self.start_correcting)
        self.pushButton_9.clicked.connect(self.peak_detecting)
        self.pushButton_tar.clicked.connect(self.select_target_excel)
        self.targeted_extraction_file = None
        self.peaks_table.itemSelectionChanged.connect(self.update_plot)

        self.figure = Figure(figsize=(5, 3), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.layout_frame_2 = QtWidgets.QVBoxLayout(self.frame_2)
        self.layout_frame_2.addWidget(self.canvas)

        self.pushButton_10.clicked.connect(self.export_to_excel)
        self.pushButton_11.clicked.connect(self.export_pickle)
        self.pushButton_12.clicked.connect(self.import_pickle)
        self.pushButton_13.clicked.connect(self.one_click_processing)
        self.one_step_processing_active = False

        # Crop retention time
        self.lineEdit_7.setEnabled(False)
        self.lineEdit_8.setEnabled(False)
        self.checkBox.setChecked(True)
        self.checkBox.stateChanged.connect(self.toggle_crop_rt)

        # Set Columns of files
        self.tableWidget.setColumnWidth(0, 900)
        self.tableWidget.setColumnWidth(1, 100)

    def change_page(self, i):
        self.button_list[i].setChecked(True)
        self.stackedWidget.setCurrentIndex(i)

    def browse_files(self):
        self.file_list.clear()
        self.qc_file_list.clear()
        self.comboBox.clear()
        self.comboBox_4.clear()
        self.tableWidget.setRowCount(0)

        file_dialog = QFileDialog()
        files, _ = file_dialog.getOpenFileNames(self,
                                                "Select files",
                                                "",
                                                "MzML Files (*.mzML);;MzXML Files (*.mzXML)",
                                                options=QFileDialog.DontUseNativeDialog)
        for file in files:
            self.file_list.append(file)

            row_position = self.tableWidget.rowCount()
            self.tableWidget.insertRow(row_position)
            self.tableWidget.setItem(row_position, 0, QTableWidgetItem(file))

            qc_checkbox = QCheckBox()
            qc_checkbox.setChecked(False)  # initially not a QC file
            qc_checkbox.stateChanged.connect(
                lambda state, current_file=file: self.update_qc_file_list(state, current_file))
            self.tableWidget.setCellWidget(row_position, 1, qc_checkbox)

    def update_qc_file_list(self, state, file):
        if state == QtCore.Qt.Checked:
            if file not in self.qc_file_list:
                self.qc_file_list.append(file)
                self.comboBox.addItem(file)
                self.comboBox_4.addItem(file)
        elif state == QtCore.Qt.Unchecked:
            if file in self.qc_file_list:
                self.qc_file_list.remove(file)
                index = self.comboBox.findText(file)
                self.comboBox.removeItem(index)
                index_4 = self.comboBox_4.findText(file)
                self.comboBox_4.removeItem(index_4)

    def toggle_crop_rt(self, state):
        if state == 0:  # Unchecked state
            self.lineEdit_7.setEnabled(True)
            self.lineEdit_8.setEnabled(True)
        else:
            self.lineEdit_7.setEnabled(False)
            self.lineEdit_8.setEnabled(False)

    def start_importing(self):
        self.exps.clear()
        if len(self.file_list) == 0:
            QMessageBox.warning(self, "No files selected", "Please select files to process first!")
            return

        if len(self.qc_file_list) == 0:
            QMessageBox.warning(self, "No QC files selected", "At least one QC or representative file is needed!")
            return

        self.textEdit.clear()  # Clear the textEdit before importing

        if not self.checkBox.isChecked():
            start_rt = float(self.lineEdit_7.text())
            end_rt = float(self.lineEdit_8.text())
            rt_filter = (start_rt, end_rt)
        else:
            rt_filter = None

        self.file_loader = FileLoader(self.file_list, self.exps, rt_filter)
        self.file_loader.update_signal.connect(self.textEdit.append)
        if self.one_step_processing_active:
            self.file_loader.finished_signal.connect(self.start_correcting)
        self.file_loader.start()

    def show_progress_bar(self):
        self.progressBar.setValue(0)
        self.progressBar.show()

    def hide_progress_bar(self):
        self.progressBar.hide()

    def update_file_directory_setting(self, index):
        if self.comboBox_2.itemText(index) == "True":
            self.lineEdit_6.setDisabled(False)
            self.pushButton_8.setDisabled(False)
        else:
            self.lineEdit_6.setDisabled(True)
            self.pushButton_8.setDisabled(True)

    def browse_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Directory")
        self.lineEdit_6.setText(directory)

    def start_correcting(self):
        if not self.exps:
            QMessageBox.warning(self, "No data imported",
                                "Please import correct LC-MS data file first.")
            return
        try:
            ref_file = self.comboBox.currentText()
            roa_window_width = float(self.lineEdit.text())
            if not 0 <= roa_window_width <= 70:
                QMessageBox.warning(self, "Invalid ROA Window Width",
                                    "Width of ROA Window should be between 0 and 70.")
                return
            delta_mz_roa = float(self.lineEdit_2.text())
            if not 0 <= delta_mz_roa <= 0.1:
                QMessageBox.warning(self, "Invalid ROA Delta m/z",
                                    "Allowed delta m/z of ROA should be between 0 and 0.1.")
                return
            intensity_threshold_roa = self.doubleSpinBox.value() / 100.0
            delta_mz_eic = float(self.lineEdit_4.text())
            if not 0 <= delta_mz_eic <= 0.2:
                QMessageBox.warning(self, "Invalid EIC Delta m/z",
                                    "Allowed delta m/z of EIC should be between 0 and 0.2.")
                return
            plot = self.comboBox_2.currentText() == 'True'
            write_new_files = self.comboBox_3.currentText() == 'True'
            new_file_directory = self.lineEdit_6.text()
            if plot or write_new_files:
                if not new_file_directory:
                    QMessageBox.warning(self, "No Directory Selected",
                                        "Please choose a directory for the new files/picture.")
                    return
            # Start CorrectionThread
            self.correction_thread = CorrectionThread(
                ref_file,
                roa_window_width,
                delta_mz_roa,
                intensity_threshold_roa,
                delta_mz_eic,
                plot,
                write_new_files,
                new_file_directory,
                self.exps,
                self.file_list
            )
            self.progressBar.setValue(0)
            self.correction_thread.update_signal.connect(self.textEdit.append)
            self.progressBar.show()
            if self.one_step_processing_active:
                self.correction_thread.finished_signal.connect(self.peak_detecting)
            else:
                self.correction_thread.finished_signal.connect(self.hide_progress_bar)
            self.correction_thread.progress_signal.connect(self.progressBar.setValue)
            self.correction_thread.activate_progress_signal.connect(self.show_progress_bar)
            self.correction_thread.hide_bar_signal.connect(self.hide_progress_bar)
            self.correction_thread.start()
        except Exception as e:
            import traceback
            print(f"Caught an exception: {e}")
            traceback.print_exc()

    def select_target_excel(self):
        self.targeted_extraction_file, _ = QFileDialog.getOpenFileName(self, "Select Excel file", "",
                                                                       "Excel Files (*.xlsx *.xls)")
        if self.targeted_extraction_file:
            self.targeted_ext_file.setText(os.path.basename(self.targeted_extraction_file))
        else:
            self.targeted_ext_file.setText("No file selected")

    def peak_detecting(self):
        if not self.exps:
            QMessageBox.warning(self, "No data imported",
                                "Please import correct LC-MS data file first.")
            return
        min_files_zero = self.spinBox_4.value() / 100.0
        delta_mz = self.lineEdit_5.text()
        continuous_points = self.spinBox_2.value()
        dropped_points = self.spinBox_3.value()
        maximum_peak_width = self.spinBox_5.value()
        min_i = self.lineEdit_3.text()
        ref_file = self.comboBox_4.currentText()
        lamda = self.doubleSpinBox_2.value()

        if any(text.strip() == '' for text in [delta_mz, min_i]):
            QMessageBox.warning(self, "Invalid Input", "All fields must be filled.")
            return

        delta_mz = float(delta_mz)
        min_i = float(min_i)

        if not 0 < delta_mz < 0.1:
            QMessageBox.warning(self, "Invalid Delta m/z",
                                "Allowed delta m/z should be between 0 and 0.1.")
            return
        if not 0 <= min_i:
            QMessageBox.warning(self, "Invalid minimum intensity",
                                "Minimum intensity should be greater than or equal to 0.")
            return
        if self.targeted_extraction_file:
            from pandas.api.types import is_numeric_dtype
            try:
                tar_df = pd.read_excel(self.targeted_extraction_file)
            except Exception as e:
                QMessageBox.warning(self, "Failed to read Excel file", str(e))
                return
            tar_cpd_columns = ['compound', 'm/z', 'RT(s)']
            if not all(column in tar_df.columns for column in tar_cpd_columns):
                QMessageBox.warning(self, "Incorrect Excel file format",
                                    "Please ensure it contains the required columns.")
                return
            numeric_columns = ['m/z', 'RT(s)']
            if not all(is_numeric_dtype(tar_df[column]) for column in numeric_columns):
                QMessageBox.warning(self, "Incorrect data format",
                                    "m/z and RT should be numeric.")
                return
            tar_cpd = tar_df.to_dict('records')
        else:
            tar_cpd = None
        self.peak_detection_thread = PeakDetectionThread(min_files_zero,
                                                         delta_mz,
                                                         continuous_points,
                                                         dropped_points,
                                                         self.exps,
                                                         self.file_list,
                                                         maximum_peak_width,
                                                         min_i,
                                                         ref_file,
                                                         self.qc_file_list,
                                                         lamda,
                                                         tar_cpd)
        # self.progressBar.setValue(0)
        # self.progressBar.show()
        self.peak_detection_thread.update_signal.connect(self.textEdit.append)
        self.peak_detection_thread.result_signal.connect(self.populate_peaks_table)
        self.peak_detection_thread.progress_signal.connect(self.progressBar.setValue)
        self.peak_detection_thread.finished_signal.connect(self.hide_progress_bar)
        self.peak_detection_thread.activate_progress_signal.connect(self.show_progress_bar)
        self.peak_detection_thread.start()

    def populate_peaks_table(self, peaks_res):
        self.peaks_res = peaks_res
        self.peaks_table.clear()

        # display only the first 50 files
        displayed_file_list = self.file_list[:50]

        self.peaks_table.setRowCount(len(peaks_res))
        self.peaks_table.setColumnCount(5 + len(displayed_file_list))

        headers = ['mz', 'rt', 's/n', 'max_i', 'compound'] + [os.path.splitext(os.path.basename(file))[0] for file in
                                                              displayed_file_list]
        self.peaks_table.setHorizontalHeaderLabels(headers)

        for i, peak in enumerate(peaks_res):
            row_data = [peak['mz'], peak['rt'], peak['s/n'], peak['max_intensity'], peak['compound']]
            # row_data = [peak['mz'], peak['rt'], peak['s/n'], peak['max_intensity'], '']
            for file in displayed_file_list:
                row_data.append(peak[file])

            for j, data in enumerate(row_data):
                self.peaks_table.setItem(i, j, QtWidgets.QTableWidgetItem(str(data)))

    def update_plot(self):
        try:
            self.figure.clear()
            selected_rows = self.peaks_table.selectedItems()
            if selected_rows:
                row = selected_rows[0].row()
                peak = self.peaks_res[row]
                ax = self.figure.add_subplot(111)
                peak['group'][0].plot(peak['group'][1], ax=ax)
            self.canvas.draw()
            self.figure.tight_layout()
        except Exception as e:
            import traceback
            print(f"Caught an exception: {e}")
            traceback.print_exc()

    def table_to_df(self):
        # Create a DataFrame from the complete dataset
        headers = ['mz', 'rt', 's/n', 'max_i', 'compound'] + [os.path.splitext(os.path.basename(file))[0] for file in
                                                              self.file_list]
        data = []

        for peak in self.peaks_res:
            row_data = [peak['mz'], peak['rt'], peak['s/n'], peak['max_intensity'], peak['compound']]
            row_data += [peak[file] for file in self.file_list]
            data.append(row_data)

        df = pd.DataFrame(data, columns=headers)
        return df

    def export_to_excel(self):
        df = self.table_to_df()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, "Save File", "", "Excel Files (*.xlsx)",
                                                  options=options)

        if fileName:
            if not fileName.endswith('.xlsx'):  # Add .pkd extension if it's not there
                fileName += '.xlsx'
            self.textEdit.append("Start exporting. This process may be time-consuming. Don't close the window!")
            df.to_excel(fileName, index=False)
            self.textEdit.append('{0} exported successfully.'.format(fileName))

    def export_pickle(self):
        if self.peaks_res is None:
            QMessageBox.warning(self, "Export failed",
                                "No peak data. Please process peak data first.")
            return
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, "Save File", "", "Peak Data Files (*.pkd)",
                                                  options=options)

        if fileName:
            import pickle
            if not fileName.endswith('.pkd'):  # Add .pkd extension if it's not there
                fileName += '.pkd'
            self.textEdit.append("Start exporting. This process may be time-consuming. Don't close the window!")
            with open(fileName, 'wb') as f:
                pickle.dump(self.peaks_res, f)
                pickle.dump(self.file_list, f)
            self.textEdit.append('{0} exported successfully.'.format(fileName))

    def import_pickle(self):
        file_dialog = QFileDialog()
        file, _ = file_dialog.getOpenFileName(self, "Select files", "", "Peak Data Files (*.pkd)")
        if file:
            import pickle
            try:
                self.textEdit.append("Start importing. This process may be time-consuming. Don't close the window!")
                with open(file, 'rb') as f:
                    peaks_res = pickle.load(f)
                    self.file_list = pickle.load(f)
                self.populate_peaks_table(peaks_res)
                self.textEdit.append('{0} imported successfully.'.format(file))
            except (OSError, IOError) as e:
                QMessageBox.warning(self, "File Error", str(e))
            except pickle.UnpicklingError:
                QMessageBox.warning(self, "File Error", "The selected file is not a valid peak data (.pkd) file.")

    def one_click_processing(self):
        reply = QtWidgets.QMessageBox.question(self, 'Confirmation',
                                               "Please confirm that all parameters in the file import, "
                                               "data correction, and peak detection steps have been set and are "
                                               "complete.",
                                               QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                               QtWidgets.QMessageBox.No)
        if reply == QtWidgets.QMessageBox.Yes:
            self.one_step_processing_active = True
            self.start_importing()


# @profile
# def main_function():
#     app = QtWidgets.QApplication([])
#     window = MainWindow()
#     window.show()
#     app.exec_()


if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    window = MainWindow()
    window.show()
    app.exec_()
    # main_function()
