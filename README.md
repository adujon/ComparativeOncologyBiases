# ComparativeOncologyBiases

This repository contains the codes and raw datasets used to perform the analyses in the following publication:

The cost of fame for popular species: strong risk of biases in comparative oncology of captive species

Dujon AM, Courtalon J, Asselin K, Thomas F. 2026 The Cost of Fame For Popular Species: Strong Risk of Bias in Comparative Oncology of Captive Species. Proceedings of the Royal Society B: Biological Sciences (doi:https://doi.org/10.1098/rspb)

Link to the publication: https://royalsocietypublishing.org/rspb/article/293/2072/20260504/482052/The-cost-of-fame-for-popular-species-strong-risk

Corresponding author: antoine.dujon@yahoo.fr

--------------------------------------------------------------------------------------------------------------------------------------------------

The rawdata folder contains the datasets used in the publication

The phylogenies folder contains the phylogenetic trees used in the Bayesian phylogenetic models

The Regressions folder contains the R scripts used to fit the model for each taxa and study (codes are fully annotated)

Raw cancer prevalence data was sourced from the following publications

Kapsetaki, S. E., Z. T. Compton, J. Dolan, V. Κ. Harris, W. Mellon, S. M. Rupp, E. G. Duke, T. M. Harrison, S. Aksoy, M. Giraudeau, O. Vincze, K. J. McGraw, A. Aktipis, M. Tollis, A. Μ. Boddy, and C. C. Maley. 2024. Life history traits and cancer prevalence in birds. Evol Med Public Health, doi: 10.1093/emph/eoae011

Kapsetaki, S.E., Compton, Z.T., Mellon, W., Vincze, O., Giraudeau, M., Harrison, T.M., et al. (2024b). Germline mutation rate predicts cancer mortality across 37 vertebrate species. Evol Med Public Health, 12, 122–128. https://doi.org/10.1093/emph/eoae016

Vincze, O., F. Colchero, J.-F. Lemaître, D. A. Conde, S. Pavard, M. Bieuville, A. O. Urrutia, B. Ujvari, A. M. Boddy, C. C. Maley, F. Thomas, and M. Giraudeau. 2022. Cancer risk across mammals. Nature 601:263–267 https://doi.org/10.1038/s41586-021-04224-5

Compton, Z. T., W. Mellon, V. K. Harris, S. Rupp, D. Mallo, S. E. Kapsetaki, M. Wilmot, R. Kennington, K. Noble, C. Baciu, L. N. Ramirez, A. Peraza, B. Martins, S. Sudhakar, S. Aksoy, G. Furukawa, O. Vincze, M. Giraudeau, E. G. Duke, S. Spiro, E. Flach, H. Davidson, C. I. Li, A. Zehnder, T. A. Graham, B. V. Troan, T. M. Harrison, M. Tollis, J. D. Schiffman, C. A. Aktipis, L. M. Abegglen, C. C. Maley, and A. M. Boddy. 2024. Cancer Prevalence across Vertebrates. Cancer Discovery, https://doi.org/10.1158/2159-8290.CD-24-0573

Dujon, A.M., Vincze, O., Lemaitre, J.-F., Alix-Panabières, C., Pujol, P., Giraudeau, M., et al. (2023b). The effect of placentation type, litter size, lactation and gestation length on cancer risk in mammals. Proceedings of the Royal Society B: Biological Sciences, 290, 20230940. https://doi.org/10.1098/rspb.2023.0940

IMPORTANT: Vincze et al. provide two cancer-mortality metrics: the incidence of cancer mortality, or ICM, and the cumulative mortality risk, or CMR. ICM is a survival/competing-risk estimate of the probability that an adult individual dies from cancer, while accounting for age structure and censored individuals. Vincze et al. designed this metric to reduce biases arising from left truncation and right censoring, including individuals that were still alive at the time of data extraction and whose eventual cause of death was therefore unknown. For this reason, we used ICM as our main cancer-mortality metric.
Because the ICM values were provided as point estimates, we converted them into derived event counts so that they could be analysed using binomial regression models. Specifically, we calculated:
derived cancer-mortality events = round(ICM × number of necropsied individuals)

--------------------------------------------------------------------------------------------------------------------------------------------------

Phylogenies were obtained from TimeTree websiste (https://timetree.org/)

Bibliometric data was obtained form the Scopus database (https://www.scopus.com) 

Wikipedia page views were obtained using the Wikimedia website (https://commons.wikimedia.org/) using the custom made script "Get Wikipedia page views.R". The script uses ISO language codes to access Wikipedia pages in different languages and download the number of page views for each species based on its scientific (binomial) name.

Economic data was obtained from the International Monetary Fund website (https://www.imf.org/)
