# MetCohort
![](./images/logo.png "MetCohort")

MetCohort is an untargeted liquid chromatography-mass spectrometry 
(LC-MS) data processing tool for large-scale metabolomics and exposomics. MetCohort can 
realize automatic correction of retention time and m/z between samples and precise 
feature detection. With innovative and robust data correction and feature detection algorithm,
MetCohort have a low false positive and false negative rate simultaneously. Feature table of high quality is generated after whole procedures, which 
significantly improves subsequent feature annotation and statistical analysis.

MetCohort now supports mzML or mzXML file format of LC-MS raw data. At
least one file need to be specified as quality control (QC) file. 
Then data correction and peak detection should be performed in order. To optimize the processing results, adjusting some parameters 
is necessary. Finally, feature detection results are visualized in the software and can be saved for 
subsequent analysis.

Detailed information about the algorithm and parameters are available in our article.

## Software interface

Data import:
![](./images/1.png "Software interface")
Data correction:
![](./images/2.png "Data correction")
Feature detection:
![](./images/3.png "Feature detection")
Results visualization:
![](./images/4.png "Results visualization")


## User guide
### Step 1: Download source code of MetCohort
Download source code to local system from [Github](https://github.com/JunYang2021/MetCohort) page.

### Step 2: Install Anaconda
If you don't already have Anaconda installed on your system, follow these steps:

1. Download the Anaconda installer for your operating system from the [Anaconda download page](https://www.anaconda.com/products/individual).
2. Run the installer and follow the on-screen instructions.

### Step 3: Set up the environment
After installing Anaconda, set up the Python environment for running MetCohort.
1. Open Anaconda Prompt (or your terminal).
2. Create a new environment using the following command:
```shell
cd path_to_MetCohort   # the directory you downloaded MetCohort repository, windows system should drive to specific hard drive first
conda env create -f environment.yml
conda activate metcohort_env
```

### Step 4: Run MetCohort
With configured environment in Anaconda, MetCohort GUI can be activated with following command in Anaconda prompt:
```shell
cd path_to_MetCohort   # the directory you downloaded MetCohort repository, windows system should drive to specific hard drive first
cd ./src
python main_panel.py
```


## Parameters introduction
### Data import
In Data Import stage, all the pending files (in mzML or mzXML format) need to be uploaded 
to the window. The QC files can be labeled. ROI detection is only performed on the labelled
QC files, which determine the construction of ROI matrix and ranges of features. Users can
select all the files as QC or only the representative files. For some experiments having 
blank or undesired signal at the beginning or ending of chromatographic gradient, users 
can uncheck the box of Auto range and set the real retention time in Crop retention time 
in seconds.


### Data alignment
In Data Alignment stage, users should select one reference file from labelled QC files 
for ROA detection. Related parameters are shown below:


_**Width of ROA window**_: The time width of the ROA window. Default is 30 seconds.


_**Allowed delta m/z of ROA**_: Allowed m/z deviation in the process of ROA detection. 
Because the ROA width is short, a relatively narrow value should be set. The default 
value is 0.001.


_**Intensity threshold of ROA**_: An intensity coefficient to control the numbers of detected
ROAs. Specifically, the intensity at the center of a detected ROA should exceed a dynamic
specified value, which is calculated as the current EWMA (Exponentially weighted moving 
average) of TIC (total ion chromatograms) multiplied by the coefficient. The default 
coefficient is 0.5%.


_**Allowed delta m/z of XIC**_: Allowed m/z deviation in the process of XIC extraction. 
It should be larger than allowed delta m/z of ROA to reach a good matching between files. 
The default value is 0.02.


_**Plot**_: Choose True or False to export the retention time deviation results of processed 
files in a html file. Abnormal files can be found from the resulting plot.


**_Write to new files_**: Choose True or False to export the raw data after data alignment
in mzML format. The exported files can be processed with other tools of LC-MS data 
processing.


_**New file directory**_: Select the exporting directory of new plots or files.


### Feature detection
In Feature Detection stage, related parameters are shown below:


_**Minimum number of non-zero peaks in a feature**_: Non-zero peaks should not exceed the 
specified proportion of all the samples. Default value is 80%.


_**Delta m/z**_: Allowed m/z deviation in the process of ROI detection and ROI matrix 
construction. Default value is 0.01. Actually, it can be adjusted to be narrower
following data alignment to enhance the peak resolution in feature detection.


_**Minimum number of continuous non-zero points**_: Allowed minimum number of continuous 
non-zero points in the process of ROI detection with centWave algorithm. Default
value is 3. In small-scale sample processing with a few QC files, the value should be
increased to decrease noise.


_**Maximum number of continuous zero points**_: Allowed maximum number of continuous zero 
points in the process of ROI detection with centWave algorithm. Default value is 10. 
It can control the number of ROI matrix through controlling the time length of detected 
ROIs.


_**Maximum width of peak**_: Allowed maximum chromatographic width of detected features. 
Default value is 15 seconds.


_**Minimum intensity of a feature**_: Allowed minimum intensity of a feature. 
If the maximum height of a feature in all the samples is lower than the value, 
the feature is filtered. Default value is 50.


_**Select a reference file**_: Select one reference file from the labelled QC 
files for feature representation. The retention time and m/z of one feature is 
represented by the retention time and m/z of the apex of the peak in the reference
file. It can be different from the reference file in raw data alignment. However, 
selecting the same file with reference file in raw data alignment is recommended,
because the feature representation can keep the same with any processing parameters 
with a fixed reference file. It can ease the process of feature matching between 
feature tables.


_**Entropy coefficient**_: A value ranging from 0 to 1 that controls the score 
threshold for feature determination in the ROI matrix. A larger entropy coefficient
corresponds to a higher threshold. The default value is set to 0.8. For large-scale
sample processing, users are advised to adjust this value lower.


_**Targeted extraction (optional)**_: A table of targeted extraction compounds in 
xlsx format. Users can optionally upload a table here for targeted extraction. 
The detailed style of the table is shown below:
![](./images/tar_emp.png "Targeted Extraction")


### View Results
After the feature detection stage, users can view the results on the View results
page of the window. When a feature is selected in the feature table, its chromatograms
are displayed. To conserve memory usage, chromatograms for a random selection of 
up to 50 files are shown. The shaded region in the chromatograms represents the 
integration region of the feature. Users can export the feature table in xlsx format 
using the button Export feature data (*.xlsx) in the window. Additionally, users can 
export and import the feature data within the window to save features and view them
in future sessions with Export feature data (*.pkd) and Import feature data (*.pkd)
buttons.

