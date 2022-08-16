# Meniscus_T1rho_and_T2star
Matlab code for calculating T1rho and T2* for knee meniscus segmented in OsiriX.

The analysis sequence is similar in structure to the sequence in MRI_T1rho_T2star_reliability repository.  Please see the ImageAnalysisPipeline6.pdf file.

There are two analyzes:  1. Using the entire segmented meniscus (mri_m_fit.m) and 2. Using the segmented meniscus with the top and bottom pixels removed (eroded) (mri_em_fit.m).

Please see the comments in the M-files.
