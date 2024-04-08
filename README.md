# BrainAGE

![](BrainAgeScheme.png)

The current BrainAGE approach using Gaussian Process Regression (GPR) is described in the paper [BrainAGE: Revisited and reframed machine learning workflow](https://doi.org/10.1002/hbm.26632).

We were the first to introduce this idea using machine learning approaches to estimate BrainAGE in [2010](https://doi.org/10.1016/j.neuroimage.2010.01.005).

### Necessary steps to use BrainAGE

Although it is of course possible to use any grey and white matter segmentation, it is recommended to use CAT12 for pre-processing, as some CAT12 functions are necessary for BrainAGE and many of the tools are adapted for use with CAT12.

## 1. Preprocessing your MRI data

Download [CAT12](https://neuro-jena.github.io/cat) and get more detailed information [here](https://github.com/ChristianGaser/cat12) and use the [CAT12 manual](http://141.35.69.218/cat12/CAT12-Manual.pdf) and the following chapters to start with CAT12:
[Getting Started](https://neuro-jena.github.io/cat12-help/#get_started)
[Quick Start Guide](https://neuro-jena.github.io/cat12-html/cat_starting.html)

For BrainAGE, we use the affine registered segmentations for grey and white matter. Therefore, you can use the CAT12 defaults for segmentation and change the following options to speed up processing, as we can skip non-linear registration and surface extraction (unless you also want to analyse the voxel- and surface-based morphometry data):
- enable 'Grey matter -> DARTEL export -> Affine'
- enable 'White matter -> DARTEL export -> Affine'

If you don't want to analyze voxel-based morphometry data:
- disable 'Process Volume ROIs'
- disable 'Grey matter -> Modulated normalized'
- disable 'White matter -> Modulated normalized'
- disable 'Bias, noise and global intensity corrected T1 image -> Normalized'
- disable 'Deformation Fields -> Image->Template (forward)'

If you don't want to analyze surface-based morphometry data:
- disable 'Surface and thickness estimation'

You can also use the CAT12 shell scripts:

        cat_batch_cat.sh -ns -nm -rp your_T1_data.nii

The used flags are:

        -ns    skip surface and thickness estimation
        -nm    skip estimating modulated and warped segmentations and ROI measures
        -rp    additionally estimate affine registered segmentations

Finally, carefully check the quality of the pre-processed data using the sample-homogeneity tool in CAT12. More information can be found [here](https://neuro-jena.github.io/cat12-help/#module4).

## 2. Organize Pre-processed Data

The abbreviations 'rp1' and 'rp2' are used for grey and white matter segmentations that have been affinely registered. In order to prepare the data for the next step, the rp1 and rp2 files should be organised in separate folders that also indicate the CAT12 version used. Copy (or move) all rp1- and rp2 files to the respective folder:

        YourDataFolder/rp1_CAT12.9
        YourDataFolder/rp2_CAT12.9

## 3. Resample and Smooth Pre-processed Data

The function `BA_data2mat` can be used to prepare pre-processed data for machine learning analysis by saving spatially registered volumes as Matlab data matrices. The function applies a mask to the volume data to remove non-brain areas, ensuring that only relevant brain information is included in the output and resamples and smoothes the data with different sizes (i.e. 4/8mm resampling, 4/8mm smoothing).

        D.age   = load(YourAgeTextFile);  % Load your age information
        D.male  = load(YourMaleTextFile); % Load information about male/female (1/0)
        
        D.release='_CAT12.9';      % Release or version information of the data
        D.name = 'YourDataName';   % Base name for the output .mat file.
        D.data = {YourDataFolder}; % Cell array of strings; paths to data folders to be concatenated.
        
        BA_data2mat(D); % call BA_data2mat to save mat-files of resampled and smoothed data

The function will save segmented, smoothed, and resampled volumes as .mat files in the specified or default directories. For example, using the above parameters, the output files will be:

        s8_8mm_rp1_YourDataName_CAT12.9.mat
        s8_8mm_rp2_YourDataName_CAT12.9.mat
        s4_8mm_rp1_YourDataName_CAT12.9.mat
        s4_8mm_rp2_YourDataName_CAT12.9.mat
        s8_4mm_rp1_YourDataName_CAT12.9.mat
        s8_4mm_rp2_YourDataName_CAT12.9.mat
        s4_4mm_rp1_YourDataName_CAT12.9.mat
        s4_4mm_rp2_YourDataName_CAT12.9.mat

If you have defined multiple folders for `D.data` (e.g. D.data = {YourDataFolder1,YourDataFolder2}), the data from the specified folders will be concatenated (first 'YourDataFolder1', then 'YourDataFolder2') and saved in a single .mat file format for each segmentation. This can be useful if your data consists of different groups (i.e. controls and patients). Please take care that the order of all pre-processed files always corresponds to the age and male information.

## 4. BrainAGE Estimation


<https://www.markdownguide.org>
<fake@example.com>

