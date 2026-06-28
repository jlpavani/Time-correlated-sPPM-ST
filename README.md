# Time-correlated-sPPM-ST

Data and R code to support Pavani et al (2026).  In this paper, we adopt a probabilistic approach to time-dependent regionalization, conceptualizing it as a random partition of geographic space at each time point, with the sequence of spatial partitions exhibiting time dependency. This methodology facilitates inference regarding the temporal dynamics of clusters. We employ a product partition prior for the random partitions at each time point, introducing temporal correlation among partitions through the temporal structure associated with prior cohesions. To explore partition search space effectively and ensure spatially constrained clustering, we utilize random spanning trees. This research is motivated by a pertinent applied problem: the identification of spatial and temporal patterns associated with mosquito-borne diseases. Given the overdispersion inherent in this type of data, we propose a spatio-temporal Poisson mixture model in which both mean and dispersion parameters vary according to spatio-temporal covariates. We apply the proposed model to analyze weekly reported cases of dengue from 2018 to 2023 in the Southeast region of Brazil. Additionally, we assess modeling performance using simulated data. Results indicate that our model is competitive in analyzing the temporal evolution of spatial clustering.

## Citation

If you find this code helpful and use it in your work, please cite our paper:

> **Pavani, J.**; Loschi, R. H.; Quintana, F. A.: Modeling temporal dependence in a sequence of spatial random partitions driven by spanning trees: an application to mosquito-borne diseases. *Annals of Applied Statistics*, 20(2), 1388-1408, 2026. [[DOI](https://doi.org/10.1214/26-AOAS2172)]

```bibtex
@Article{Pavani2025,
  author  = {Pavani, Jessica and Loschi, Rosangela Helena and Quintana, Fernando Andr{\'e}s},
  journal = {Annals of Applied Statistics},
  title   = {Modeling temporal dependence in a sequence of spatial random partitions driven by spanning trees: an application to mosquito-borne diseases},
  year    = {2026},
  number  = {2},
  pages   = {1388--1408},
  volume  = {20},
  doi     = {10.1214/26-AOAS2172}
}
