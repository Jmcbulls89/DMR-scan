# DMR-scan

The DMR-scan pipeline was for the purpose of identifying differentially methylatyed regions in "Parental experience modifies the Mimulus methylome" (Colicchio et al., BMC Genomics, 2018) to identify regions of differential methylation from WGBS data.  
This pipeline provides an extension to analyses after reads have been mapped (for example by BMap: http://itolab.med.kyushu-u.ac.jp/BMap/usage.html ) followed by nucleotide percent methylation smoothing using methimpute (BMC Genomics, Taudt et al., 2018 : https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4641-x , https://github.com/ataudt/methimpute/tree/master/R).  It utilizes a PELT changepoint detection approach followed by explicit glmer testing with a bobyqa equalizer to call consistent regions of differential methylation.

The code included here contains the pipeline used to analyze data as well as a simulated data set that can be used to compare different paramter options or could be used to compare this method with other existing methodologies.
