# dataconla2023
Resources for my talk at Data Con LA 2023: "Predicting Purchases, Rare Diseases, and More: Using Ordinal Regression to Estimate Rare Event Probabilities"

* The slides are available in the file **Data Con LA 2023.pdf**.
* The file **pres_figs.R** contains the R code that generates the figures in slides 7 - 13. (It will automatically load needed functions from the file **model_functions.R**. You do not need to open the file **model_functions.R** manually, just make sure both of these files are in the same folder on your computer when you run **pres_figs.R**.)
* The file **proportional_odds_code.R** contains the R code from slides 21 - 27 that train a proprtional odds model on the soup data set. It also contains the example code from slide 24 that demonstrates how to create an ordered categorical variable in R.

## Other resources
* For more resources related to the method [PRESTO](https://proceedings.mlr.press/v202/faletto23a.html) developed by me and Prof. Jacob Bien at USC (including code implementing PRESTO), please see [this GitHub repo](https://github.com/gregfaletto/presto).
* [This link](https://stats.oarc.ucla.edu/r/dae/ordinal-logistic-regression/) has more helpful information about implementing the proportional odds model in R, and [this link](https://www.statsmodels.org/stable/examples/notebooks/generated/ordinal_regression.html#Logit-ordinal-regression:) has more information about implementing it in Python.
* If you have any other questions, feel free to reach out to me on [Twitter](https://twitter.com/GregoryFaletto), [LinkedIn](https://www.linkedin.com/in/gregfaletto/), or by email at [gregory.faletto@marshall.usc.edu](mailto:gregory.faletto@marshall.usc.edu)!
