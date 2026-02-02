library(ggplot2)
library(e1071) # to calculate empirical skewness and kurtosis

mass_predictions <- read.csv('Data/mass_predictions.csv')
mass_predictions$hum.fem.circ..mm. <- as.double(mass_predictions$hum.fem.circ..mm.)

specimens <- read.csv('Data/specimens.csv', skip = 17, sep='')

head(mass_predictions)

#################
### PLOT DATA ###
#################
# remove missing values from (log) circumferences
x <- mass_predictions$hum.fem.circ..mm.
x <- x[!is.na(x)]
x_log <- log10(x)

# Find the maximum of x and use its index to get it out of mass_predictions
max_x <- max(x)
max_index <- which(mass_predictions$hum.fem.circ..mm. == max_x)
max_species <- mass_predictions[max_index, 'genus.and.species']
cat("Maximum C_FH:", max_x, "\n")
cat("Species with max C_FH:", max_species, "\n")

# Find the second maximum of x and use its index to get it out of mass_predictions
second_max_x <- sort(x, decreasing = TRUE)[2]
second_max_index <- which(mass_predictions$hum.fem.circ..mm. == second_max_x)
second_max_species <- mass_predictions[second_max_index, 'genus.and.species']
cat("Second Maximum C_FH:", second_max_x, "\n")
cat("Species with second max C_FH:", second_max_species, "\n")

min_x <- min(x)
min_index <- which(mass_predictions$hum.fem.circ..mm. == min_x)
min_species <- mass_predictions[min_index, 'genus.and.species']
cat("Minimum C_FH:", min_x, "\n")
cat("Species with min C_FH:", min_species, "\n")

## Strong negative skewness
ggplot(mass_predictions, aes(x=log10(hum.fem.circ..mm.))) + geom_histogram(aes(y=after_stat(density)), bins=30) + geom_density(color='red', linewidth=1)

# Add secondary x-axis with exponentiated values
ggplot(mass_predictions, aes(x=log10(hum.fem.circ..mm.))) + 
  geom_histogram(aes(y=after_stat(density)), bins=30) + 
  geom_density(color='red', linewidth=1) +
  scale_x_continuous(
    name = "log10(hum.fem.circ..mm.)",
    sec.axis = sec_axis(~10^(.), name = "hum.fem.circ..mm.")
  )

# g = m_3/m_2^{3/2} (where m are then sample moments ...) 
skewness(x)
skewness(x_log) #-0.46

## Strong deviations from normal distribution -> need other parametric methods!
ggplot(mass_predictions, aes(sample=hum.fem.circ..mm.)) + geom_qq() + geom_qq_line()
# shorter left tail, shorter right tail it seems
ggplot(mass_predictions, aes(sample=log10(hum.fem.circ..mm.))) + geom_qq() + geom_qq_line()
# log transformed data gives the negative skewness. Doesn't solve the non normality problem
# kan ook met Kolmogorov-Smirnov test (meer algemeen, ik ken Shapiro Wilk eigenlijk niet)
shapiro.test(mass_predictions$hum.fem.circ..mm.) # p=0.04
shapiro.test(log10(mass_predictions$hum.fem.circ..mm.)) #p=0.007

saveRDS(x_log, file = "Data/logCFH")

