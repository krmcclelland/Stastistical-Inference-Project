# Statistical Inference Project
#  Part 1 Simulation Exercise

## Overview
The Statistical Inference Project will analyze and compare exponential distribution and the Central Limit Theorem (CLT) in R.  Exponential distribution is a probability analyses of a continuous event over period of time measuring the lower and higher ends of something that predicts some type of occurrence being measure in the future. The properties are continuous and moves to the right.  There are a couple of means we will look at; Sample Means and Theoretical Means.  Sample Means by using facts and observations whereas, Theoretical Means deals with the what-if's or variables.  

In R we can simulate mean by $(\bar{x} = \frac{\sum x_i}{n})$ using the function means("data"), and theoretical using 'rexp(n, $\lambda$)' function, with lambda $\lambda$ representing the rate, mean ($\mu$) or the average of the data, is $\mu = \frac{1}{\lambda}$, and standard deviation $\sigma = \frac{1}{\lambda}$.     

The Central Limit Theorem states as long as there is a sample size of 30 or greater, the skewed data elements within the data-set will normalize and reshape into a bell curve over time. The formula is $\bar{x}_{n}$ ~ $N(\mu, \frac{\sigma }{\sqrt{n}})$

In comparing the exponential distribution and Central Limit Theorem, the sample means of $1000$ in groups of 40 $(n = 40)$, and standard deviation $({s} = \sqrt{s^2})$  of 0.2 $(\lambda = 0.2)$ is approximately normal using $N(\frac{1}{0.2}, \frac{1}{0.2}{\sqrt{40}})$.  

```{r Data Prep, echo = FALSE}
set.seed(200)
lambda <- 0.2
n <- 40
sim <- 1000
```

```{r Simulation Results, include=FALSE}
sim_exp <- replicate(sim, rexp(n, lambda))
head(sim_exp)
```

```{r Simulation Exponential Results, include=FALSE}
means_exp <- apply(sim_exp, 2, mean) 
head(means_exp)
```
## Simulation
Mean Comparison
```{r Mean Comparision, echo=FALSE}
samtheor_mean = NULL
for (i in 1:1000) {
  samtheor_mean = c(samtheor_mean, mean(rexp(n, lambda)))
}

head(samtheor_mean)
```
## Requirement 1: Compare Sample vs Theoretical Mean
Mean Actual
```{r Mean Actual Results, echo=FALSE}
act_mean <- mean(means_exp)
act_mean
```
Mean Theory
```{r Mean Theory Results, echo=FALSE}
theo_mean <- 1/lambda
theo_mean
```

``` {r "Sample vs Theoretical Means"}
library(ggplot2)
samtheor_mean_df <- as.data.frame(samtheor_mean)
g <- ggplot(samtheor_mean_df, aes(x = samtheor_mean))
g <- g + geom_histogram(binwidth = .2, color = "black", fill = "dark green") +
  geom_vline(xintercept = theo_mean, color = "yellow") +
  geom_vline(xintercept = act_mean, color = "red") +
  labs(x = "40 Simulated Exponential Samples", y = "Frequency", title = "Sample vs. Theoretical Mean")
g
```
![Sample-vs-theoretical Mean](https://user-images.githubusercontent.com/17502725/139456921-d56264bf-0c4e-4b82-9067-dff57ebd481e.png)


## Requirement 2: Compare Variance Distribution of Sample vs. Theoretical
```{r Standard Deviation Actual, include=FALSE}
act_sd <- sd(means_exp)
act_sd
```

```{r Standard Deviation Theory, include=FALSE}
theo_sd <- ((1/lambda) * (1/sqrt(n)))
theo_sd
```
Theory Variance
``` {r Theory Variance Results, echo = FALSE}
theo_var <- theo_sd^2
theo_var
```
Actual Variance
``` {r Actual Variance Results, echo = FALSE}
act_var <- var(means_exp)
act_var
```

```{r Distribution Approximation to Normal, echo=FALSE}
tab_den <-matrix (c(theo_mean, act_mean, 
                  theo_var, act_var), 
                  ncol = 2, byrow = TRUE)
colnames(tab_den) <- c("Theoretical", "Sample")
rownames(tab_den) <- c("Mean", "Variance")
tab_den <- as.table(tab_den)
tab_den
```
In comparison of the above table both the Theoretical and Sample Means are within close proximity to one another. 

## Requirement 3: Distribution Approximation to Normal
```{r Dist Approx to Normal, echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE}
samtheor_mean_df <- as.data.frame(samtheor_mean)
d <- ggplot(samtheor_mean_df, aes(x=samtheor_mean))
d <- d + geom_histogram(binwidth = .2, color = "black", fill = "orange", aes(y = ..density..)) +
    stat_function(fun = dnorm, args = list(mean = theo_mean, sd = sd(samtheor_mean)), 
    color = "red", size = 1) +
    stat_density(geom = "line", color = "blue", size = 1) +
    labs(x = "Mean of 40 Simulated Exponential Samples", y = "Density", title = "Density Simulation")
d
```
![Density Simulation](https://user-images.githubusercontent.com/17502725/139457673-2553f6a7-bbe0-42b6-838d-492d8e53e932.png)

The proximity of the density (blue) line and the normal (red) line demonstrates the distribution is near normal.  

## Verification of Confidence Levels
```{r Veri Conf Levels, echo=FALSE}
round((act_mean + c(-1, 1)* qnorm(.975)*sd(samtheor_mean)/sqrt(sim)),3)
```
The confidence levels above are within the 95% of the theoretical means of 5.0, and between the above range.

******
# Part 2 - Guinea Pig Odontoblasts Cell Tooth Growth - Basic Inferential Data Analysis

## Overview
Tooth Growth basic analysis is to evaluate whether Vitamin C effects the tooth growth in Guinea Pigs.  The data consist of 60 observation with three (3) variables; length, supplement, and dose.  Supplements (supp) were given in one of three levels (0.5 mg, 1.0 mg, and/or 2.0mg ) to the Guinea Pigs either given by orange juice (OJ) or ascorbic acid (VC).  The project - part two requirements are:

* Load the ToothGrowth data set and perform some basic exploratory data analyses
* Provide a basic summary of the data.
* Use confidence intervals and/or hypothesis tests to compare tooth growth by supp and dose.
* State your conclusions and the assumptions needed for your conclusions.

## Loading Libraries: library(ggplot2)| library(UsingR) | library(datasets) | library(dplyr)
```{r Load Libraries, include=FALSE}
library(ggplot2)
library(UsingR)
library(datasets)
library(dplyr)
```

## Requirement 1 - Loading Data
### Sampling of Data Within ToothGrowth Data Set
```{r TG Data Frame, echo=FALSE}
data(ToothGrowth)
TG <- ToothGrowth
str (TG)
```
### Review Sampling of Data Set Content
```{r TG Head, echo=FALSE}
head(TG)
```

```{r Dose Levels, include=FALSE}
unique(TG$dose)
```

## Requirement 2 - Summary of Data St Content
```{r Data Content, echo = FALSE}
summary(TG)
```
### Dose Mean
```{r Dose Mean, echo=FALSE}
md = split(TG$len, TG$dose)
sapply(md, mean)
```

```{r TG Dose Variable as a Factor, echo=FALSE}
TG$dose <- as.factor(TG$dose)
```
### tapply Length / Supplement
```{r tapply, echo=FALSE}
tapply(TG$len, TG$supp, var)
```

## Requirement 3 - Comparison of Tooth Length to Supplement and Dose 

### Supplement Response Rate
```{r Supplement Response Rate, echo=FALSE, fig.height=3, fig.width=7}
TGPlot_supp <- ggplot(data = TG, aes(x = supp, y = len)) + geom_boxplot(aes(fill = supp), color = "black")

TGPlot_supp
```
![ComparisonofToothLengthtoSupplementandDose](https://user-images.githubusercontent.com/17502725/139458568-9b66fc54-f5c5-4584-84cb-015b30835360.png)

### Supplement Dose by Type and Ondontoblasts Cell Growth
```{r Supp Dose by Type and Tooth Length Growth, echo=FALSE, fig.height=3.85, fig.width=7}
dat = ToothGrowth
coplot(len ~ dose | supp,dat = dat, panel = panel.smooth, xlab = "Supplement Dose and Type", 
       ylab = "Length of Tooth Growth") 
```
![SupplementDosebyTypeandOndontoblastsCellGrowth](https://user-images.githubusercontent.com/17502725/139460140-ff05481f-1371-4cef-ac7d-28b7d0c03344.png)

### Length of Odontoblast Cell vs. Dose
```{r TGPLOT, echo=FALSE, fig.height=2.75, fig.width=7}
library(ggplot2)
TGPlot <- ggplot(data = TG, aes(x = dose, y = len)) +
            geom_boxplot(aes(fill = dose)) + facet_grid(~supp) +
            xlab("Dose") + ylab("Length of Odontoblasts Cell") 
TGPlot
```
![DosebyLengthOndontoblastsCell](https://user-images.githubusercontent.com/17502725/139459398-f1697ff1-41d8-42f3-8213-2362749826a6.png)

```{r Dose Filter, include=FALSE}
library(dplyr)
d1 <- filter(TG, dose == 0.5)
d2 <- filter(TG, dose == 1.0)
d3 <- filter(TG, dose == 2.0)
```

## t.test Summary
* P-Value 
    * 0.5 dose = 0.006359
    * 1.0 dose = 0.001038
    * 2.0 dose = 0.9639
* 95% Confidence Levels
    * 0.5 dose = 1.719057 / 8.780943
    * 1.0 dose = 2.802148 / 9.057852
    * 2.0 dose = -3.79807 / 3.63807
* Sample Means in OJ Group / VC Group
    * 0.5 dose = 13.23 / 7.98 
    * 1.0 dose = 22.70 / 16.77
    * 2.0 dose = 26.06 / 26.14
    

```{r pvalue Comparision, include=FALSE}
t.test(len ~ supp, paired = F, var.equal = F, data = d1)
t.test(len ~ supp, paired = F, var.equal = F, data = d2) 
t.test(len ~ supp, paired = F, var.equal = F, data = d3)
```

## Conclusion
Assuming the ToothGrowth data-set representation of the population, and the Sample Means is consistent with the Central Limit Theorem (CTL), it can be concluded that the supplements have a positive correlation on the growth of the Guinea Pigs Odontoblasts cell length in the increasing dosage.  It appears that vitamin supplements do not grow the Odontoblasts Cell as fast as Orange Juice when given in 0.5 mg and / or 1.0 mg doses as it does in 2.0 mg dose.  The p-values and means are within normal expected values and a NULL hypothesis can be rejected.
