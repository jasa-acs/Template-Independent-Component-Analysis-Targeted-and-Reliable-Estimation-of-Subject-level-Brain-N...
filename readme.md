# Template Independent Component Analysis: Targeted and Reliable Estimation of Subject-level Brain Networks Using Big Data Population Priors


# Author Contributions Checklist Form

## Data

### Abstract

Study subjects included healthy adult participants from the Human Connectome
Project (HCP) 500-subject (HCP-500) and 900-subject (HCP-900) releases. For all subjects, the resting-state functional magnetic resonance (fMRI) data was analyzed, which was minimally preprocessed and projected to the cortical surface using the HCP fMRISurface pipeline. The HCP data releases include this preprocessed data. In addition, group-level independent component analysis (ICA) maps from the HCP 500 Subjects Parcellation, Timeseries, Netmats (HCP 500 PTN) release were also utilized. More information about the HCP data releases is available at https://www.humanconnectome.org/study/hcp-young-adult/data-releases.



### Availability 

The HCP data is publicly available. To download it, first register with ConnectomeDB at https://db.humanconnectome.org/ then accept the HCP data use terms at http://connectome. uardev.com/study/hcp-young-adult/data-use-terms. There are several
options for obtaining the data, including through Amazon S3, Aspera Connect, and Connectome in a Box (now retired).

The HCP-500 and HCP-900 are available as part of the WU-Minn HCP Data - 1200 Subjects dataset, linked at the top of the main ConnectomeDB page (https://db.humanconnectome.org/). Using Aspera Connect, imaging data can be downloaded for specific subsets of subjects (e.g., S500 Release Subjects, S900 Release Subjects) through the “Download Image Data” menu from the main page. It is also possible to download data for individual subjects by entering the “Explore Subsets” menu on the main page, selecting a subject ID, and clicking “Download Images” at the top of the page. The files used in this analysis are available in the Resting State fMRI FIX-Denoised (Compact) package.

The HCP data is organized in a nested directory structure. For each subject [SUB], the
necessary files are located [SUB]/MNINonLinear/Results. In our analysis, we accessed the HCP data through Connectome in a Box, in which subjects are organized into hard disks, so there is also a disk directory in the file path (i.e., [DISK]/[SUB]/MNINonLinear/Results). To perform the analysis using data obtained through Amazon S3 or Aspera Connect, one would simply exclude the disk directories in the file paths. Within the Results directory for each subject, for this analysis we utilized four resting-state fMRI sessions:
• rfMRI_REST1_LR/ rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii (419MB)
• rfMRI_REST1_RL/ rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii (419MB)
• rfMRI_REST2_LR/ rfMRI_REST2_LR_Atlas_hp2000_clean.dtseries.nii (419MB)
• rfMRI_REST2_RL/ rfMRI_REST2_RL_Atlas_hp2000_clean.dtseries.nii (419MB)


### Description
Separate subsets of HCP subjects were used for (1) template estimation and (2)
testing the performance of template ICA and other methods.

Template estimation. For template estimation, the HCP-500 was used, specifically 461 subjects for whom all 4 resting-state sessions of data were collected. Note that the HCP-500 is a subset of the HCP-900.

In addition, a set of group ICA maps at the Q=25 resolution available as part of the HCP500-PTN were used in template estimation. This data release is described at
https://www.humanconnectome.org/study/hcp-young-adult/article/beta-release-group-ica-mapsnode-timeseries-network-matrices but appears to be no longer be available for download. We therefore share the file containing the group ICA maps used in this analysis, named melodic_IC.dscalar.nii (9MB) and originally located in HCP500_Parcellation_Timeseries_Netmats /groupICA/groupICA_3T_Q1-
Q6related468_MSMsulc_d25.ica in the HCP500-PTN release.

Testing. For testing, twenty subjects were randomly sampled from the HCP-900, excluding the 461 subjects used for template estimation. These subject IDs are 645450, 917558, 121416, 664757, 390645, 656253, 894067, 926862, 287248, 134425, 825048, 449753, 168745, 188448, 820745, 129331, 130821, 894774, 156536, and 178647.

## Code

### Abstract

Template ICA is implemented in the templateICA MATLAB toolbox, available for download at https://github.com/mandymejia/templateICA. The code to perform the analysis presented in this paper is divided into three different code files. First, in the MATLAB script Templates.m, templates are estimated based on a training set of subjects. Second, in the MATLAB script TemplateICA_main.m, we apply the template ICA algorithm to a separate set of test subjects. Both MATLAB scripts rely on the templateICA toolbox. Finally, the R script visualize.R creates plots for the paper.

### Description

The templateICA toolbox contains functions to perform template estimation and apply the
template ICA algorithm. It is available at https://github.com/mandymejia/templateICA. A
README file (https://github.com/mandymejia/templateICA/blob/master/README.md) describes the functions and directory structure of the toolbox.
The code to perform the analysis presented in the paper consists of three files, to be executed in the following order. Arrows (“à”) indicate where intermediate helper files are produced. Data files are indicated in italics; functions or code files are indicated in bold.

1. Templates.m

Estimate templates (mean and variance maps) using resting-state fMRI data from a
large group of training subjects and make plots for Figure 1.
a. Identify the HCP-500 subjects with 4 resting-state fMRI sessions. à subjects.mat
b. Read in group ICA maps from melodic_IC.dscalar.nii (9MB).
c. Read in resting-state fMRI data for each subject and session. Perform dual
regression to obtain subject-level ICA map estimates at both visits.
d. Calculate mean and between-subject variance maps (the “templates”). à
maps_mean.dtseries.nii (6.5MB), maps_var.dtseries.nii (6.5MB)
e. Fits mixture of gaussians (MoG) distribution to the values in each mean map (for
Figure 1).

Relies on the following toolboxes and data files:
o templateICA – https://github.com/mandymejia/templateICA
o cifti-matlab – https://github.com/Washington-University/cifti-matlab
o Group ICA maps (melodic_IC.dscalar.nii, 9MB)
o HCP resting-state fMRI sessions from HCP-500 (see Data Availability)
Computation time: The main computational task is, for each of the 461 subjects, reading
in the data and performing dual regression to get IC estimates. This takes ~2.3 minutes
per subject, for a total of 18 hours if the subjects were run in sequence. Using a parfor
loop with 12 parallel workers, the computation time is reduced to approximately 2 hours.
The remaining tasks (computing the mean and variance maps) are very fast.


2. TemplateICA_main.m
Estimate subject-level ICs using template ICA and dual regression and quantify reliability
of the estimates.
a. Read in templates (mean and variance maps)
b. Randomly sample 20 test subjects from HCP-900 (excluding training subjects) à
subjects_samp.mat
c. For each subject, estimate ICs based on different scan durations
(T=400,800,1200) using the following methods:
o Dual regression
o Template ICA (fast EM)
à S_DR_*00.mat (250MB/each), S_tICA_*00.mat (250MB/each), *=4,8,12
d. Compute reliability measures of IC estimates produced from each method à
ICC_DR_*00.dscalar.nii (6MB/each), ICC_tICA_*00.dscalar.nii (6MB/each),
I2C2_*00.csv (<1KB/each), *=4,8,12

Relies on the following toolboxes/functions and data files:
o templateICA – https://github.com/mandymejia/templateICA
o cifti – Contains ciftiopen.m, ciftisavereset.m from the HCP. These functions
rely on the Connectome Workbench and the gifti toolbox. See instructions at
https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ#HCPUs
ersFAQ-2.HowdoyougetCIFTIfilesintoMATLAB?
o Group ICA maps (melodic_IC.dscalar.nii, 9MB)
o HCP resting-state fMRI sessions from 20 test subjects (see Data Availability)
Computation time: Running template ICA for each subject and visit requires, on average,
approximately 5 minutes for T=400, 10 minutes for T=800 and 20 minutes for T=1200.
Subjects and visits can be parallelized to reduce the total computation time. Reading in
the data, performing dual regression, and computing the reliability of the IC estimates
are very fast.

3.visualize.R, in which plots for the figures are created.
Computation time: negligible

We also used the wb_view application of the Connectome Workbench for visualization of the CIFTI files produced through the analysis (ending in .dscalar.nii and .dtseries.nii). This software is publicly available at http://www.humanconnectome.org/software/get-connectome-workbench.

Note that visualization requires the user to supply a pair (for the left and right hemispheres) of GIFTI surface files (ending in surf.gii) to map the values in the CIFTI files to the locations on the cortical surface. We used the files Q1-Q6_R440.*.inflated.32K_fs_LR.surf.gii, *=L,R. These are
available through the Workbench v1.0 Tutorial Dataset, access to which is described in the
tutorial document available at https://www.humanconnectome.org/storage/app/media/
documentation/tutorials/Connectome_WB_Tutorial_v1.0.pdf. For convenience, we also include both files here.

## Instructions for Use

### Reproducibility

All Figures related to the data analysis described in the paper can be
reproduced by accessing the HCP data for the subjects listed below, downloading Connectome Workbench, and executing the scripts above. To download HCP data, you must register with ConnectomeDB at https://db.humanconnectome.org. To visualize CIFTI files (Figures 20 and 21), Connectome Workbench is required – see instructions at the end of Code Description.