# ComparativeOncologyBiases

This repository contains the codes and raw datasets used to perform the analyses in the following publication:

The cost of fame for popular species: strong risk of biases in comparative oncology of captive species

Dujon AM, Courtalon J, Asselin K, Thomas F. 2026 The Cost of Fame For Popular Species: Strong Risk of Bias in Comparative Oncology of Captive Species. Proceedings of the Royal Society B: Biological Sciences (doi:https://doi.org/10.1098/rspb)

Link to the publication: https://royalsocietypublishing.org/rspb/article/293/2072/20260504/482052/The-cost-of-fame-for-popular-species-strong-risk

If you are accessing this project through the Zenodo repository, you can also find it on GitHub, where additional post-publication reanalyses may be posted: https://github.com/adujon/ComparativeOncologyBiases

Corresponding author: antoine.dujon@yahoo.fr

--------------------------------------------------------------------------------------------------------------------------------------------------

The rawdata folder contains the datasets used in the publication

The phylogenies folder contains the phylogenetic trees used in the Bayesian phylogenetic models

The Regressions folder contains the R scripts used to fit the model for each taxa and study (codes are fully annotated)

The post-publication reanalysis folder contain addtional analyses conducted after recieving feedback on our published manuscript (see bellow for details). Those reanalyses confirm our results.

--------------------------------------------------------------------------------------------------------------------------------------------------
Data sources:

Raw cancer prevalence data and other life history trait metrics were sourced from the following publications

Kapsetaki, S. E., Z. T. Compton, J. Dolan, V. Κ. Harris, W. Mellon, S. M. Rupp, E. G. Duke, T. M. Harrison, S. Aksoy, M. Giraudeau, O. Vincze, K. J. McGraw, A. Aktipis, M. Tollis, A. Μ. Boddy, and C. C. Maley. 2024. Life history traits and cancer prevalence in birds. Evol Med Public Health, doi: 10.1093/emph/eoae011

Kapsetaki, S.E., Compton, Z.T., Mellon, W., Vincze, O., Giraudeau, M., Harrison, T.M., et al. (2024b). Germline mutation rate predicts cancer mortality across 37 vertebrate species. Evol Med Public Health, 12, 122–128. https://doi.org/10.1093/emph/eoae016

Vincze, O., F. Colchero, J.-F. Lemaître, D. A. Conde, S. Pavard, M. Bieuville, A. O. Urrutia, B. Ujvari, A. M. Boddy, C. C. Maley, F. Thomas, and M. Giraudeau. 2022. Cancer risk across mammals. Nature 601:263–267 https://doi.org/10.1038/s41586-021-04224-5

Compton, Z. T., W. Mellon, V. K. Harris, S. Rupp, D. Mallo, S. E. Kapsetaki, M. Wilmot, R. Kennington, K. Noble, C. Baciu, L. N. Ramirez, A. Peraza, B. Martins, S. Sudhakar, S. Aksoy, G. Furukawa, O. Vincze, M. Giraudeau, E. G. Duke, S. Spiro, E. Flach, H. Davidson, C. I. Li, A. Zehnder, T. A. Graham, B. V. Troan, T. M. Harrison, M. Tollis, J. D. Schiffman, C. A. Aktipis, L. M. Abegglen, C. C. Maley, and A. M. Boddy. 2024. Cancer Prevalence across Vertebrates. Cancer Discovery, https://doi.org/10.1158/2159-8290.CD-24-0573

Dujon, A.M., Vincze, O., Lemaitre, J.-F., Alix-Panabières, C., Pujol, P., Giraudeau, M., et al. (2023b). The effect of placentation type, litter size, lactation and gestation length on cancer risk in mammals. Proceedings of the Royal Society B: Biological Sciences, 290, 20230940. https://doi.org/10.1098/rspb.2023.0940

IMPORTANT: Vincze et al. provide two cancer-mortality metrics: the incidence of cancer mortality, or ICM, and the cumulative mortality risk, or CMR. ICM is a survival/competing-risk estimate of the probability that an adult individual dies from cancer, while accounting for age structure and censored individuals. Vincze et al. designed this metric to reduce biases arising from left truncation and right censoring, including individuals that were still alive at the time of data extraction and whose eventual cause of death was therefore unknown. For this reason, we used ICM as our main cancer-mortality metric.
Because the ICM values were provided as point estimates, we converted them into derived event counts so that they could be analysed using binomial regression models. Specifically, we calculated:
derived cancer-mortality events = round(ICM × number of necropsied individuals)

--------------------------------------------------------------------------------------------------------------------------------------------------

Phylogenies were obtained from TimeTree website (https://timetree.org/)

Bibliometric data was obtained from the Scopus database (https://www.scopus.com) using the citation overview tool after searching for the species scientific names (aka their binomial name) between quotation marks: for example Ailurus fulgens is searched as "Ailurus fulgens".

Wikipedia page views were obtained using the Wikimedia website (https://commons.wikimedia.org/) using the custom made script "Get Wikipedia page views.R". The script uses ISO language codes (contained in the ISO Codes .csv file) to access Wikipedia pages in different languages and download the number of page views for each species based on its scientific (binomial) name. 

Economic data was obtained from the International Monetary Fund website (https://www.imf.org/)

--------------------------------------------------------------------------------------------------------------------------------------------------
POST-PUBLICATION REANALYSES

Following feedback received on our manuscript after publication, we are now providing additional analyses to assess the robustness of our study. These analyses, together with the corresponding code, figures, and raw datasets, have been uploaded to the “post-publication reanalysis” folder.

We reran the univariate models linking cancer-mortality risk to species popularity using CMR-based metrics, rather than ICM, following the same methodology described in the paper. This was done to ensure that the choice of ICM did not substantially affect the results. These analyses produced results similar to those obtained using ICM, which is consistent with Vincze et al., who reported good agreement between ICM and CMR. The same variables that were positively associated with mortality risk when using ICM were also positively associated with mortality risk when using CMR.

One of our popularity metrics was calculated using a custom R script that queried Wikimedia page-view data, and this script has now been uploaded to the GitHub repository. Our approach used scientific names, rather than common names, across language-specific Wikipedia pages. We chose this strategy to reduce issues associated with translating common names across languages and to provide a proxy for broad public interest in each species. The values used in our study were therefore based on scientific names and language-specific Wikipedia pages, rather than common-name searches. Searches based on common names can return different pages, redirects, and cross-language mappings, which can produce substantially different page-view totals. The website https://pageviews.wmcloud.org/langviews provides a tool called Langviews, which allows users to search common names across Wikipedia pages and obtain total page-view numbers. To evaluate whether our results were robust to this alternative way of estimating species popularity, we collected page-view estimates for each species in the Vincze et al. dataset using the Langviews tool from the Wikimedia Pageviews website. We used the date range 01/07/2015 to 31/12/2025, selected “en.wikipedia.org” as the source project, “all” as the platform, “user” as the agent, and searched using the English common name of each species. We then compared the page-view estimates obtained from our original proxy with those obtained using the Langviews-based proxy and found good agreement between the two approaches. We also reran the phylogenetic binomial regression models predicting cancer-mortality risk using both the original proxy and the Langviews-based proxy, and the results were highly consistent. These analyses indicate that both approaches capture the same broad signal: species with higher public popularity, estimated using Wikipedia page-view data, tend to have higher reported cancer-mortality risk.

To assess whether our conclusions depended our representation ICM which was transformed as a count for the binomial framework we used, we conducted an additional analysis using the ZIlogis procedure provided in the supplementary material of Vincze et al using the point estimate ICM values and considering this variable as continuous. This approach separates the analysis into two components. The first component models whether cancer was detected in a species in at least one individual as a binary response, with species classified as having no recorded cancer mortality (0) or recorded cancer mortality (1). Each species-level observation in this first component was weighted by the logarithm of the number of necropsies. The second fits a phylogenetic Gaussian regression to the logit-transformed ICM estimates for species in which cancer was detected, using the weighting scheme implemented by Vincze et al., with weights calculated as 1/log⁡("number of necropsies") and the same model of trait evolution that they used (Pagel’s model). Separate univariate models were fitted for each popularity metric. Coefficients from both components were exponentiated and are presented with their 95% confidence intervals. These analyses showed that species popularity was significantly associated with ICM when modelled as a continuous variable. Among species in which cancer mortality was detected, the number of Wikipedia page views, residual Wikipedia page views, and GDP per capita were positively associated with ICM. By contrast, the number of animals held in captivity and the number of zoological institutions in which a species was represented were negatively associated with ICM within this cancer-detected subset. Code is available on the GitHub repository of the project along with all the reanalyses. 

The effect sizes obtained from the two-part analysis are not directly comparable with those from the frequency-weighted binomial analysis we conducted in our publication because the models estimate different quantities, treat zero-valued ICM estimates differently, and use necropsy sample size in different ways. The purpose of this analysis is not to claim that the two approaches estimate identical quantities, but to determine whether associations with popularity were also recovered when ICM was analysed using the two-part procedure provided by Vincze et al.. The binomial models relate reconstructed cancer-mortality frequencies to popularity across the complete set of species while retaining information on both cancer occurrence and cancer mortality within a single model. By contrast, the ZIlogis approach decomposes the data into a binary cancer-detection process and a continuous mortality process conditional on cancer detection. This decomposition by the ZIlogis function results in some loss of information: in the first component all non-zero ICM values are collapsed into a simple cancer detected/not detected classification, while in the second component species with zero ICM values are excluded entirely. Consequently, species with zero ICM contribute to the detection component but not to the continuous component. This distinction is particularly relevant because cancer detection may itself be associated with species popularity, monitoring intensity and the amount of veterinary information available. Associations estimated only among species in which cancer was detected can therefore differ in magnitude or even direction from associations estimated across all species. In addition, the two approaches use necropsy sample size differently, altering the relative influence of data-rich and data-poor species.
