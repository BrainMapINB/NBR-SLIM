# NBR-SLIM

Here are the scripts and preprocessed data to illustrate the R package NBR using the dataset SWU-SLIM.

The [Scripts](https://github.com/BrainMapINB/NBR-SLIM/tree/main/Scripts) directory contains the scripts to:
* [01-Preprocessing](https://github.com/BrainMapINB/NBR-SLIM/blob/main/Scripts/01-Preprocessing.R) preprocesses the data founded in <http://fcon_1000.projects.nitrc.org/indi/retro/southwestuni_qiu_index.html>.
* [02-Analyses](https://github.com/BrainMapINB/NBR-SLIM/blob/main/Scripts/02-Analyses.R) contains the commands to generate the inputs to run NBS ([Zalesky et al., 2010](https://doi.org/10.1016/j.neuroimage.2010.06.041)) and [NBR](https://CRAN.R-project.org/package=NBR) models in the preprocessed data.
* [03-Visualization](https://github.com/BrainMapINB/NBR-SLIM/blob/main/Scripts/03-Visualization.R) plots the results for each model tested.

The [Data](https://github.com/BrainMapINB/NBR-SLIM/tree/main/Data) directory contains the preprocessed phenotypic and Fisher's Z functional connectivity grouped by lobe anatomy.

The [Analyses](https://github.com/BrainMapINB/NBR-SLIM/tree/main/Analyses) directory contains the inputs to run NBS, as well as the results after being run.
