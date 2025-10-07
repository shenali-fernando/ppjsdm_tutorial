# ppjsdm_tutorial
This repo contains a four-part in-depth tutorial for using the ppjsdm package that contains the saturated pairwise interaction Gibbs points process model with an example dataset. I recommend starting with a simple example to get a feel for how the model works with a dataset, this can be found in the vignette 'Part0_simpledata'. This vignette runs through two relatively simple dataset where one can easily fit the model. 

Parts 1-4 of this tutorial show how to apply this model more complex datasets with multiple species. These parts of the tutorial use ForestGEO data 

Part 1 gives package installation help as well as code and explanation of the various parameters and how to assign values to these. Part 2 shows how to fit the model and some common errors that might occur at this stage. Part 3 gives examples of visualisation of the coefficients including coefficient plots, heat maps, chord diagrams, and also how to extract coefficients into a dataframe to customise other plots. Lastly Part 4 shows how to predict using the model, and goes through computing and understanding the log-Papangelou conditional intensity, computing AUC scores and comparing to the Poisson model. 

More examples of fitting the ppjsdm model to other datasets are given in the folder marked 'other_examples'. These examples tend to be simpler systems and use data freely available from R. Lastly, the PPJSDM_model_overview.doc gives an explanation of the model and its parameters and how everything comes together. 
