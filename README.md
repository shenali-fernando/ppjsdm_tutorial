# ppjsdm_tutorial
This repo contains a four-part in-depth tutorial for using the ppjsdm package that contains the saturated pairwise interaction Gibbs points process model with an example dataset.

We recommend starting with the vignette 'Part0_simpledata'. This vignettes gives an example of applying the model to simple datasets, where one can easily fit the model, in order to get a feel for how the model works.

We also recommend reading the [paper](https://doi.org/10.1111/rssc.12596) that introduces this model. The PPJSDM_model_overview also gives a plain-language summary of the workings of the model, a description of the model equation, and some notes on model behaviour.

Parts 1-4 of this tutorial shows how to apply this model to more complex datasets. These parts of the tutorial use the Smithsonian Institute's ForestGEO plot data. ForestGEO is a global network of fully-mapped forest plots, where all censuses have been conducted in the same way so that the data for all plots have a similar structure. As this model works best on fully-mapped, static datasets, the ForestGEO datasets are particularly useful. 

Parts 1-4 of the tutorial are set out in this way: 

* Part 1 gives package installation help and data loading and cleaning. It then gives an in-depth explanation of the various parameters and how to specify these
* Part 2 shows how to fit the model and some common errors that might occur at this stage
* Part 3 gives examples of visualisation of the coefficients including coefficient plots, heat maps, chord diagrams, and also how to extract coefficients into a dataframe to customise other plots
* Part 4 shows how to predict using the model, and goes through computing and understanding the log-Papangelou conditional intensity, computing AUC scores and comparing to the Poisson model

More examples of fitting the ppjsdm model to other datasets are given in the folder marked 'other_examples'. These examples tend to be simpler systems and use data freely available from R. 
