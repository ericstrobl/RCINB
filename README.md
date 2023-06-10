# Root Causal Inference with Negative Binomials (RCI-NB)

RCI-NB is an algorithm for discovering patient-specific root causes of disease from single cell RNA-sequencing (scRNA-seq) data. scRNA-seq datasets contain non-negative counts measured by an error prone sequencing process. RCI-NB accounts for count-based measurement error by separating negative binomial distributions into their gamma and Poisson components; the gamma distributions form a fully identifiable but latent post non-linear causal model representing the true RNA expression levels, which we only observe with Poisson corruption. RCI-NB identifies patient-specific root causal contributions from scRNA-seq datasets by integrating regression and goodness of fit testing procedures that bypass Poisson measurement error.

The ``Experiments`` folder contains code to replicate the experimental results in the paper. Please cite the article if you use any of the code in this repository.

# Installation

Press the green button up top and download the zip file. Then:

> library(devtools)

> install_local("Directory of RCINB-main.zip")

> library(RCINB)

# Run the Algorithm

Generate the grouth truth DAG with 10 variables and an expected neighborhood size of two:
> DAG = generate_NB_DAG_r(10,2)

Generate ten thousand samples from the DAG:
> dataAll = sample_NB_DAG_r(10000, DAG)

Run RCI-NB:
> out = RCI_NB(data$data,DAG$Y,data$C,data$Xp)

Print expected Shapley values:
> print(out$shaps)
