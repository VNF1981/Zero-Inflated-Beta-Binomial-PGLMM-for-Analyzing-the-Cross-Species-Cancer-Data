## LFSR - False Discovery Rate for Bayesian Models 
Since we run separate Bayesian models for individual predictors, there is no partial pooling across the predictor coefficients to account for multiple comparisons. Therefore, we should account for multiple testing before reporting significant associations. Traditional methods such as Bonferroni and Benjamini-Hochberg are designed for frequentist analyses based on p-values. Instead, we will use a Bayesian approach based on the local false sign rate (lfsr), following Stephens (2017).
The method works as follows:
1.	For each predictor, we calculate the posterior probabilities of its regression coefficient being greater or smaller than zero.
2.	We calculate the lfsr, which represents the posterior probability of incorrectly identifying the direction of the association.
3.	We rank the predictors based on their lfsr values in ascending order and calculate their cumulative mean.
4.	We identify the largest set of predictors for which the cumulative false sign rate remains at or below 5%.
5.	As an additional conservative criterion, we require that the 95% credible interval excludes zero before reporting an association as statistically supported.
This approach allows us to account for multiple comparisons while staying within the Bayesian framework. More specifically, it controls the posterior expected proportion of reported associations with an incorrect direction of effect. The script automatically identifies the completed models, extracts the posterior samples, performs the calculations, checks model diagnostics, and generates the results in CSV and Excel formats. It can handle any number of predictors and performs the correction separately for each response variable by default. Note that there is no need to rerun the ZiBBPGLMM models, as the script uses the existing posterior samples stored in the output files. I have checked the statistical calculations and compatibility with the current ZiBBPGLMM source code. However, the script still needs to be tested on your actual model outputs.

Reference for your papers:
Stephens, M. (2017). False discovery rates: a new deal. Biostatistics, 18(2), 275–294. https://doi.org/10.1093/biostatistics/kxw041

