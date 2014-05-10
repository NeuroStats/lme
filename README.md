lme
===

Linear mixed effects (lme) Matlab toolbox. Implements both univariate and mass-univariate analyses. Jorge Luis Bernal Rusiel, 2012.

If you use these tools in your analyses please cite:

Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu M.R., 2013. Spatiotemporal Linear Mixed Effects Modeling for the Mass-univariate Analysis of Longitudinal Neuroimage Data, Neuroimage 81C, 358-370.

Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear Mixed Effects Models, Neuroimage 66C, 249-260.

in addition to (for the longitudinal image processing):

Reuter M., Schmansky N.J., Rosas H.D., Fischl B, 2012.Within-Subject Template Estimation for Unbiased Longitudinal Image Analysis, NeuroImage 61(4), pp. 1402-1418, 2012, http://dx.doi.org/10.1016/j.neuroimage.2012.02.084

These Matlab tools are freely distributed and intended to help neuroimaging researchers when analyzing longitudinal neuroimaging (LNI) data. The statistical analysis of such type of data is arguable more challenging than the cross-sectional or time series data traditionally encountered in the neuroimaging field. This is because the timing associated with the measurement occasions and the underlying biological process under study are not usually under full experimental control.

There are two aspects of longitudinal data that require correct modeling: The mean response over time and the covariance among repeated measures on the same individual. I hope these tools can serve for such modeling purpose as they provide functionality for exploratory data visualization, model specification, model selection, parameter estimation, inference and power analysis including sample size estimation. They are specially targeted to be used with Freesurfer's data but can be used with any other data as long as they are loaded into Matlab and put into the appropriate format. Here are some recommendations about how to use these tools.

The lme Matlab toolbox is also distributed within Freesurfer since Freesurfer's v5.2 release.  An optional sample dataset which can be used to become familiar with the LME Matlab tools can be found at ftp://surfer.nmr.mgh.harvard.edu/pub/data/long_mixed_effects_tools-data.tar.gz

OUTLINE
1-Preparing your data
2-Model specification 
3-Parameter estimation 
4-Model selection
5-Inference
6-Power analysis
7-Example data analyses



1-Preparing your data

There are two types of analyses that can be done: univariate and mass-univariate. The first step is to load your data into Matlab. If you are working with Freesurfer then univariate data (eg. Hippocampal volume) can be loaded using Qdec tables. There are, under the Qdec directory, some simple example scripts for reading, modifying and writing Freesurfer's Qdec tables. 

In order to read mass-univariate data you should use the following scripts:

fs_read_label.m
fs_read_surf.m
fs_read_Y.m

The last two depend on Freesurfer's scripts so you need to have installed Freesurfer software package and included the Freesurfer's matlab subdirectory in your Matlab's search path.

Previously, the mass-univariate data can be generated in Freesurfer by running variants of the following commands: 

mris_preproc --qdec qdec.table.dat --target study_average --hemi lh --meas thickness --out lh.thickness.mgh (assembles your thickness data into a single lh.thickness.mgh file)

mri_surf2surf --hemi lh --s study_average --sval lh.thickness.mgh --tval lh.thickness_sm10.mgh --fwhm-trg 10 --cortex  --noreshape (smooths the cortical thickness maps with FWHM=10 mm. Note the --cortex and --noreshape options)

Then you can load the cortical thickness data lh.thickness_sm10.mgh into Matlab using fs_read_Y.m

eg.

[Y,mri] = fs_read_Y('lh.thickness_sm10.mgh');

You should also read the spherical surface (lh.sphere) and cortex label (lh.cortex.label) of study_average.

eg.

lhsphere = fs_read_surf('$FsDir/freesurfer/subjects/fsaverage/surf/lh.sphere');
lhcortex = fs_read_label('$FsDir/freesurfer/subjects/fsaverage/label/lh.cortex.label');

Once you have your data in Matlab you need to build your design matrix. For computational efficiency reasons, these tools require the data ordered according to time for each individual (that is, your design matrix needs to have all the repeated assessments for the first subject, then all for the second and so on). You can use the script:

sortData 

For example, if you have a longitudinal Qdec table containing four columns: "fsid" (Freesurfer ID), "fsid-base" (subject-specific template name that can be used as a subject ID), "time" (eg. time from baseline) and "group"(eg. a binary variable) then you can use the following code:

Qdec = fReadQdec('qdec.table.dat');
Qdec = rmQdecCol(Qdec,1); (removes first column)
sID = Qdec(2:end,1); (grabs the subjects' IDs)
Qdec = rmQdecCol(Qdec,1); (removes the subjects'ID column)
M = Qdec2num(Qdec); (converts to a numeric matrix)
[M,Y,ni] = sortData(M,1,Y,sID); (sorts the data)

where M,Y,ni are now: the ordered assessment variables, the ordered data and a vector with the number of repeated measures for each subject respectively. A design matrix X can then be built from matrix M to represent a desired longitudinal model. For example, a simple linear model containing a group by time interaction can be obtained with the following design matrix:

X = [ones(length(M),1) M M(:,1).*M(:,2)];

and the contrast [0 0 0 1] can be used to test that interaction. 



2-Model specification 

Analysis of longitudinal data should always start by simple graphical displays of the data. A natural way is through the use of time plots. A time plot is a scatter-plot with the responses on the vertical axis and the measurement times on the horizontal axis. In general it is usually more informative to display a time plot of the smoothed mean response over time. This graphical display can be quite enlightening and can provide the basis for choosing an appropriate model for the analysis of the longitudinal changes. The toolbox provides the following two functions:

lme_timePlot

lme_lowessPlot

In the mass-univariate setting it is more difficult to determine the trend in the mean response over time since it may vary across different regions of the cortex. In the absence of an a priori hypothesis you can assume a simple linear trend or, alternatively, you can plot the longitudinal trajectories in different regions (eg.from the Freesurfer's parcelation) and select the most complex trend for your model (eg. quadratic in time over linear in time, etc...). We recommend to order the columns of your design matrix in the following way: First, the intercept term (which is a column of ones); second, the time covariate; third, any time-varying covariates (eg. time^2); fourth, the group covariates of interest (eg. a binary variable indicating whether the subject is an Alzheimer's patient or not) and their interactions with the time-varying covariates; finally any other nuisance time-invariant covariate (eg. gender). See below in section 7 for an example of a mass-univariate analysis and the corresponding model for the design matrix. 



3-Parameter estimation 

A very flexible and versatile approach for analyzing longitudinal continuous data is the linear mixed effects (LME) regression paradigm. This paradigm can provide parsimonious models for both the trend in the mean response over time and the covariance among repeated measures on the same individual. The toolbox provides three types of tools for fitting LME models: 

a)Univariate tools
lme_fit_FS (eg. lhstats = lme_fit_FS(X,[1 2],Y,ni);)
lme_fit_EM
lme_fit_NR 

b)Mass-univariate tools
lme_mass_fit_vw (eg. lhstats = lme_mass_fit_vw(X,[1 2],Y,ni,lhcortex);)
lme_mass_fit 

c)Novel mass-univariate tools (spatiotemporal models)
Spatiotemporal models are more powerful to detect effects in your data than traditional vertex-wise models when two or more random effects are included in the longitudinal statistical model. Fitting these models usually requires less computation time than the above vertex-wise mass-univariate tools. Here you should first compute initial temporal covariance component estimates using:

lme_mass_fit_init or lme_mass_fit_EMinit

these estimates can then be used to segment the brain into homogeneous regions of vertices with similar covariance parameters by:

lme_mass_RgGrow

The spatiotemporal model can then be fitted by using:

lme_mass_fit_Rgw

together with the previous segmentation and initial covariance estimates. 



4-Model selection 

Here a likelihood ratio test can be used to compare a model with q random effects against a model with q+1 random effects using

Univariate
lme_LR

Mass-univariate
lme_mass_LR (a correction for multiple comparisons is then required)
 


5-Inference

Univariate
lme_F

Mass-univariate
lme_mass_F (a correction for multiple comparisons is then required)

At this time the only tools provided for multiple comparisons correction are

lme_mass_FDR (Original Benjamini and Hochberg FDR procedure)
lme_mass_FDR2 (A powerful two-stage FDR procedure)



6-Power analysis
Only univariate tools are provided:

lme_plannedPower
lme_plannedSampleSize
lme_realizedPower



7-Example data analyses 

Here are two example analyses using ADNI data processed with Freesurfer


a)Univariate

In this case the interest was in determining the differences in mean hippocampal volume change over time among four groups of individuals: healthy controls (HC), stable mild cognitive impairment patients (sMCI), MCI patients who converted to Alzheimer's disease (cMCI) during the follow-up period and patients with Alzheimer's disease (AD) at baseline. A Qdec table file, qdec.table.dat, was used to process in Freesurfer all time point scans from all subjects included in the study. It contains the following columns:

1) Freesurfer's fsid (fsid)
2) subject ID (fsid-base)
3) time from baseline (time)
4) group (HC=1, sMCI=2, cMCI=3, AD=4)
5) Apoe4 (E4:carriers=1, non-carriers=0)
6) Gender (female=1, male=0)
7) Age at baseline (BslAge)
8) Education (in years).

The qdec.table.dat file was then used to collect volumetric data for all subjects within the Freesurfer's Qdec interface. This resulted in a new Qdec table which had the previous columns plus three new columns:

9) left Hippocampal volume
10) right Hippocampal volume
11) total intracranial volume (ICV).

The data was loaded into Matlab using:

Qdec = fReadQdec('qdec.table.dat');
Qdec = rmQdecCol(Qdec,1);
sID = Qdec(2:end,1);
Qdec = rmQdecCol(Qdec,1);
M = Qdec2num(Qdec);
Y = M(:,7:9);
M = M(:,1:6);

All the data was ordered according to time for each subject using

[M,Y,ni] = sortData(M,1,Y,sID);

where M,Y,ni are now: the ordered assessment variables, the ordered data (Hippocampal and intracranial volumes) and a vector with the number of repeated measures for each subject respectively. The mean response trends over time were then plotted

lme_lowessPlot(M(:,1),Y(:,1)+Y(:,2),0.70,M(:,2));

These plots revealed a linear trajectory over time for the four groups. So, the following LME model was proposed:

Yij = ß1 + ß2*tij +ß3*sMCIi + ß4*sMCIi *tij + ß5*cMCIi + ß6*cMCIi *tij + ß7*ADi + ß8*ADi *tij + ß9*E4i + ß10*E4i *tij + ß11*Genderi + ß12*BslAgei + ß13*Educationi + ß14*ICVi + b1i + b2i + eij

Thus, our design matrix X was built from matrix M to represent the previous model. The design matrix X then had 14 columns:

1) intercept (all ones)
2) time (tij)
3) one for sMCI and zero otherwise (sMCIi)
4) colum 3) .* time
5) one for cMCI and zero otherwise (cMCIi)
6) colum 5) .* time
7) one for AD and zero otherwise (ADi)
8) colum 7) .* time
9) E4
10) E4 .* time
11) Gender
12) BslAge
13) Education
14) ICV (converted to liters).

The parameters of this model were estimated using

total_hipp_vol_stats = lme_fit_FS(X,[1 2],Y(:,1)+Y(:,2),ni);

Note that this model contains two random effects corresponding to the intercept and slope terms (b1i, b2i). This model can be compared with the model with only a single random effect corresponding to the intercept (b1i):

total_hipp_vol_stats_1RF = lme_fit_FS(X,[1],Y(:,1)+Y(:,2),ni);

using the likelihood ratio test:

lr = lme_LR(total_hipp_vol_stats.lreml,total_hipp_vol_stats_1RF.lreml,1);

The previous test was significant (pval<0.0001) indicating that the model with two random effects is significantly better than the model with a single random effect. To test if there was any difference in the rate of change over time among the four groups the following contrast was used

C = [zeros(3,3) [1 0 0 0 0; -1 0 1 0 0; 0 0 -1 0 1] zeros(3,6)];

and the F-test

F_C = lme_F(total_hipp_vol_stats,C);

The p-value of the test was in F_C.pval.



b)Mass-univariate

In this case the interest was in determining regionally specific differences in cortical thickness atrophy rate over time among the previous four groups of individuals: HC, sMCI, cMCI and AD. Thus we had the same Qdec table of the previous example. The following Freesurfer's commands were used to generate the spatial cortical thickness data: 

mris_preproc --qdec qdec.table.dat --target fsaverage --hemi lh --meas thickness --out lh.thickness.mgh 

mri_surf2surf --hemi lh --s fsaverage --sval lh.thickness.mgh --tval lh.thickness_sm10.mgh --fwhm-trg 10 --cortex  --noreshape 

Then the spatially smoothed data lh.thickness_sm10.mgh were loaded into Matlab using

[Y,mri] = fs_read_Y('lh.thickness_sm10.mgh');

As in the previous example, the data were ordered according to time for each individual:

Qdec = fReadQdec('qdec.table.dat');
Qdec = rmQdecCol(Qdec,1);
sID = Qdec(2:end,1);
Qdec = rmQdecCol(Qdec,1);
M = Qdec2num(Qdec);
[M,Y,ni] = sortData(M,1,Y,sID);

The fsaverage's spherical surface (lh.sphere) and cortex label (lh.cortex.label) were read with:

lhsphere = fs_read_surf('$FsDir/freesurfer/subjects/fsaverage/surf/lh.sphere');
lhcortex = fs_read_label('$FsDir/freesurfer/subjects/fsaverage/label/lh.cortex.label');

A priori, we expected a linear model of cortical thickness atrophy over time. However, in order to account for any possible consistent non-linearity across vertices we carried out a model selection procedure starting with a quadratic model over time. So, the following LME model was initially proposed: 

Yij = ß1 + ß2*tij + ß3*t²ij + ß4*sMCIi + ß5*sMCIi*tij + ß6*sMCIi*t²ij + ß7*cMCIi + ß8*cMCIi*tij + ß9*cMCIi*t²ij + ß10*ADi + ß11*ADi*tij + ß12*ADi*t²ij + ß13*E4i + ß14*E4i*tij + ß15*Genderi + ß16*BslAgei + ß17*Educationi + b1i + b2i + eij

Thus, our design matrix X was built from matrix M to represent the previous model.Then initial vertex-wise temporal covariance estimates were computed with:

[lhTh0,lhRe] = lme_mass_fit_EMinit(X,[1 2],Y,ni,lhcortex,3);

These covariance estimates were segmented into homogeneous regions using:

[lhRgs,lhRgMeans] = lme_mass_RgGrow(lhsphere,lhRe,lhTh0,lhcortex,2,95);

Here both lhTh0 and lhRgMeans maps were overlaid onto lhsphere and visually compared each other to ensure that they were similar enough (the essential spatial organization of the initial covariance estimates was not lost after the segmentation procedure, otherwise the above input value 2 must be reduced to 1.8 and so on):

surf.faces = lhsphere.tri;
surf.vertices = lhsphere.coord';
figure; p1 = patch(surf);
set(p1,'facecolor','interp','edgecolor','none','facevertexcdata',Th0(1,:)');
figure; p2 = patch(surf); set(p2,'facecolor','interp','edgecolor','none','facevertexcdata',lhRgMeans(1,:)');

The spatiotemporal LME model was fitted using:

lhstats = lme_mass_fit_Rgw(X,[1 2],Y,ni,lhTh0,lhRgs,lhsphere);

Note that there is no need to pass the lhcortex label into the previous estimation function since non-cortex vertices (codified as 0 in lhRgs segmentation vector) are not automatically considered for estimation. If it is desired to compare the previous LME model with the one with a single random effect (or the one with three random effects) the same segmentation must be used to estimate the parameters for both models. Here, in order to ensure validity of both models, the segmentation obtained for the most complex model must be used. This is because the same permisibility values for the segmentation procedure result in less and larger regions for the model with lesser number of random effects and using different segmentations for different models will invalid the subsequent likelihood ratio test. For instance, to compare the previous LME model with the one with a single random effect for the intercept term, the following steps can be followed: 

i) Fit the model with one random effect using the segmentation obtained from the previous model:
lhTh0_1RF = lme_mass_fit_EMinit(X,[1],Y,ni,lhcortex,3);
lhstats_1RF = lme_mass_fit_Rgw(X,[1],Y,ni,lhTh0_1RF,lhRgs,lhsphere);

ii) Apply the likelihood ratio test:
LR_pval = lme_mass_LR(lh_stats,lhstats_1RF,1);

iii) Correct for multiple comparisons:
dvtx = lme_mass_FDR2(LR_pval,ones(1,length(LR_pval)),lhcortex,0.05,1);

In our case most cortex vertices survived the above correction (length(dvtx) > 80% of length(lhcortex)) and therefore we assumed that the model with two random effects was significantly better than the model with a single random effect.

We then tested the null hypothesis of no group differences in the quadratic term:

CM.C = [zeros(3,5) [1 0 0 0 0 0 0;-1 0 0 1 0 0 0;0 0 0 -1 0 0 1] zeros(3,5)];
F_lhstats = lme_mass_F(lhstats,CM);

and multiple comparisons were corrected with:

dvtx = lme_mass_FDR2(F_lhstats.pval,F_lhstats.sgn,lhcortex,0.05,0);
    
As dvtx was empty, which means that no vertex survived the multiple comparison correction, then all quadratic terms were removed from the model. Thus a linear model over time with two random effects, as in the univariate case, was then fitted using the above functions (lme_mass_fit_EMinit, lme_mass_RgGrow and lme_mass_fit_Rgw) and the null hypothesis of no group differences in the rate of change over time among the four groups was tested with the following contrast:

CM.C = [zeros(3,3) [1 0 0 0 0; -1 0 1 0 0; 0 0 -1 0 1] zeros(3,5)];

F_lhstats = lme_mass_F(lhstats,CM);

The Freesurfer significance map was then written to disc for visualization and post-processing in Freesurfer (eg. using tksurfer tool):

fs_write_fstats(F_lhstats,mri,'sig.mgh','sig'); 

Any ß coefficient can be saved with the code (eg. the second coefficient ß2)

nv=length(lhstats);
Beta2 = zeros(1,nv);
for i=1:nv
  if ~isempty(lhstats(i).Bhat)
     Beta2(i) = lhstats(i).Bhat(2);
  end;
end;
mri1 = mri; mri1.volsz(4) = 1;
fs_write_Y(Beta2,mri1,'Beta2.mgh');

Alternatively, a correction for multiple comparisons can be carried out directly here with

[detvtx,sided_pval,pth] = lme_mass_FDR2(F_lhstats.pval,F_lhstats.sgn,lhcortex,0.05,0);

then the sided p-values can be saved with

fs_write_Y(sided_pval,mri1,'spval.mgh');

and finally the Freesurfer's mri_surfcluster cmd can be used with spval.mgh, lh.cortex.label and pth to generate the set of detected clusters and their anatomical coordinates. 
