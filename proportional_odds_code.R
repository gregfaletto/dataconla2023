# Example of creating an ordinal categorical variable from scratch
# Categorical vector of age groupings
ages <- c(">35", "<25", ">35", "25-29", "30-35")

ages

# Specify that they are ordinal variables with the given levels
factor_ages <- factor(ages, order = TRUE, 
                             levels = c("<25", "25-29", "30-35",">35"))

factor_ages

# Training proportional odds model on soup data set

# Run this line of code if you haven't already installed the ordinal package
# install.packages(“ordinal”)
library(ordinal)

# Load the soup data set
data(soup)

# Take a look at the data
str(soup)

# Get the proportion of observations that lie in each class
summary(soup$SURENESS)/length(soup$SURENESS)

# Run this line of code if you haven't already installed the MASS package
# install.packages(“MASS”)

# Train proportional odds model
library(MASS)
model_soup <- polr(SURENESS ~ PROD + DAY + SOUPTYPE + SOUPFREQ +
                   COLD + EASY + GENDER + AGEGROUP + LOCATION,
                   data=soup)

model_soup

# Get predicted probabilities of lying in each class for each
# observation in the data set
predict(model_soup, data=soup, type="probs")
