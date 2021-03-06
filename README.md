# Technologies: 
Project is created with:
* R: version 3.5.2 (2018-12-20)
# Introduction: BENIN
This repository stores the source code for BENIN: Biologically Enhanced Network INference
BENIN infers the gene regulatory network by combining time series gene expression data with prior knowledge data such as genome wide location data (ChiP-chip), transcription factor binding sites (TFBS), perturbation data or ortholog regulatory information from closely related model organisms

The BENIN source code is in the R folder

BENIN has been applied to:
* Whitfield Human Hela data (http://genome-www.stanford.edu/Human-CellCycle/Hela/data.shtml). The source code for applying BENIN to human Hela data is accessible in: realdatanetworkinference/humanhelanetwork
* DREAM4 challenge data (http://dreamchallenges.org/project/dream4-in-silico-network-challenge/). The source code of applying BENIN to the DREAM4 challenge data is accessible on dreamfourcode/ folder

# Inspiration
We reuse and modified part of the code the R  package: boot (https://github.com/cran/boot)
We modified the "tsboot" funtion at the folling line
".....statistic(ran.gen(ts.orig[inds, ], n.sim, ran.args),inds, ...)". Note that the boot function is used to resample the time series expression data and apply BENIN on the bootstrap samples.

We further use and modified the code from R package MRIaggr to compute evaluation metrics such as area under ROC  curve, when applying BENIN to human Hela data. We modified the function "modcalcAUPRC"

# Citing

If you used BENIN for your research, please cite:

* Wonkap, S. K., & Butler, G. (2020). BENIN: Biologically enhanced network inference. Journal of Bioinformatics and Computational Biology, 18(03), 2040007.
APA
	
* Kamgnia, S., & Butler, G. (2019, December). BENIN: combining knockout data with time series gene expression data for the gene regulatory network inference. In Proceedings of the Tenth International Conference on Computational Systems-Biology and Bioinformatics (pp. 1-9).  

