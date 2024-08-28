# Peakmat v1.0
Peakmat is an untargeted liquid chromatography-mass spectrometry 
(LC-MS) data processing tool for metabolomics and exposomics. Peakmat can 
realize automatic correction of retention time and m/z between samples and precise 
feature detection. With innovative and robust data correction and feature detection algorithm,
Peakmat have a low false positive and false negative rate simultaneously. Feature table of high quality is generated after whole procedures, which 
significantly improves subsequent feature annotation and statistical analysis.

Peakmat now supports mzML or mzXML file format of LC-MS raw data. At
least one file need to be specified as quality control (QC) file. 
Then data correction and peak detection should be performed in order. To optimize the processing results, adjusting some parameters 
is necessary. Finally, feature detection results are visualized in the software and can be saved for 
subsequent analysis.

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
### Step 1: Download source code of Peakmat
Download source code to local system from [Github](https://github.com/JunYang2021/Peakmat) page.

### Step 2: Install Anaconda
If you don't already have Anaconda installed on your system, follow these steps:

1. Download the Anaconda installer for your operating system from the [Anaconda download page](https://www.anaconda.com/products/individual).
2. Run the installer and follow the on-screen instructions.

### Step 3: Set up the environment
After installing Anaconda, set up the Python environment for running Peakmat.
1. Open Anaconda Prompt (or your terminal).
2. Create a new environment using the following command:
```shell
cd path_to_peakmat   # the directory you downloaded Peakmat repository, windows system should drive to specific hard drive first
conda env create -f environment.yml
conda activate peakmat_env
```

### Step 4: Run Peakmat
With configured environment in Anaconda, Peakmat GUI can be activated with following command in Anaconda prompt:
```shell
cd path_to_peakmat   # the directory you downloaded Peakmat repository, windows system should drive to specific hard drive first
cd ./src
python main_panel.py
```
