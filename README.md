# Final-year-project-AG
Improving a direct infusion mass spectrometry data processing pipeline
File Formatting

Combining_Spectra.rmd: R Markdown script for combining spectra to ensure compatibility with MALDIquant.
MALDIquant_Formatter.ipynb: Jupyter Notebook for standardizing MALDIquant output.
DIMSpy_Formatter.ipynb: Jupyter Notebook for standardizing DIMSpy output.
MZmine_Formatter.ipynb: Jupyter Notebook for standardizing MZmine3 output.
Pipeline Runs

DIMSpy_SNR2.0_ppm223_NC_Pipeline.ipynb: Jupyter Notebook for the initial DIMSpy run with default parameters.
DIMSpy_SNR3.5_ppm555_NC_Pipeline.ipynb: Jupyter Notebook for the DIMSpy run with optimized parameters.
MALDIquant_Tol5e-6_keep400_C_Pipeline.R: R script for the MALDIquant pipeline run.
Data Processing

Normalization.R: R script for evaluating the effects of normalization on the data.
Corrected_PCA-DIMSpy.R: R script for batch correction using DIMSpy data.
Corrected_PCA-MALDIquant.R: R script for batch correction using MALDIquant data.
Corrected_PCA-MZmine.R: R script for batch correction using MZmine data.
Statistics

Clustering.R: R script for clustering the data and extracting base peaks.
Univariate_Testing.ipynb: Jupyter Notebook for univariate testing and multiple corrections.
Volcano_Plots.R: R script for creating volcano plot figures.
Venn_Diagram.ipynb: Jupyter Notebook for creating Venn diagram figures.
