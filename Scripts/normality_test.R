# Read data(modify data file) 
data <- read.table("/media/d/9C50B69250B6731E/HHCLab/Malaria/Pf7/ngs.sanger.ac.uk/production/malaria/Resource/34/Pf7_vcf/Results/Scikit-allel/TRAP_Dxy.tsv", header = FALSE, sep = "\t")
head(data)
# Create a histogram
hist(data[, 3], main="Histogram")

# Create a Q-Q plot
qqnorm(data[,3])
qqline(data[,3])

# Test for normality 
shapiro.test(data[, 3])

Africa_Africa <- data[1:12,]
Africa_Africa
Asia_Asia <- data[13:32,]
Asia_Asia
Africa_Asia <- data[33:52,]
Africa_Asia

# Create group labels based on table names
group_labels <- rep(c("Africa_Africa", "Asia_Asia", "Africa_Asia"), times = c(nrow(Africa_Africa), nrow(Asia_Asia), nrow(Africa_Asia)))

# Combine the data from all tables into a single dataframe
all_data <- rbind(Africa_Africa, Asia_Asia, Africa_Asia)


#install.packages("car")
#library(car)             # Loads 'car'

# Perform Levene's test
#levene_test <- leveneTest(all_data[,3] ~ group_labels)

# View the test results
#print(levene_test)

library(tidyverse) 

anova_result <-aov(all_data[,3]~ group_labels, data = data)

Africa_t_test <-t.test(Africa_Africa[,3], Africa_Asia[,3])
Asia_t_test <-t.test(Asia_Asia[,3], Africa_Asia[,3])

anova_result
Africa_t_test
Asia_t_test

