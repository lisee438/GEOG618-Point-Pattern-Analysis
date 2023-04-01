# GEOG618 - Point Pattern Analysis
GEOG618 Point Patter Analysis Term Project. Formulated by Lisa Tang, Arghavan Zarandi, and Tirtha Harshilkumar Gajjar.

# Point Pattern Analysis
Point Pattern Analysis is a topic under spatial analysis that can support explanations on various point patterns. 
![github-centrography](https://user-images.githubusercontent.com/118564598/228592282-5f54fa85-f859-4f84-9510-6b906ea7198a.svg)

The following workshop follows the Intro to GIS and Spatial Analysis text, specifically referring to chapter 11, which can be found here: https://mgimond.github.io/Spatial/chp11_0.html 

## Necessary Tools and Packages to be used:
This workshop will be conducted in R-studio and the dataset to be used is the "spatstat" package available in R. An R file is included as well to provide access to the entire code. To prepare for this workshop, please have R readily downloaded and the "spatstat" packaged installed so that it can be easily loaded. 

```Ruby
install.packages("spatstat")
require(spatstat)
data("paracou")
class("paracou")
str("paracou")
summary("paracou")
```

## 1. Plot Data Points
First we will plot our data points to visualize our dataset and understand what is happening spatially:
```Ruby
plot(paracou, cols=2:3, chars=c(16,3), main = "Locations of Kimboto Trees")
```
![Location Plot](https://user-images.githubusercontent.com/118564598/229259222-705b7ff8-1bba-41a0-bc5e-6f71c9e47e13.png)

## 2. Global Global Intensity
To understand the basic global pattern of our points, we can compute a simple function to determine the pattern's overall density. 
```Ruby
density <- summary(paracou)$intensity
```
Running the global intensity code will provide the estimated point density value of 0.004205303, which tells us that the average intensity is 0.004205303 points per square metre. 

## 3. Local Density
However, computing the global density does not provide us with much information on the local patterns. To understand the more complex point patterns that will be addressed in this dataset, we will go through an overview of the G- and F-functions to provide us with more answers.

**Null Hypothesis:**
If the Kimboto Trees are not normally distributed, then the data follows the Poisson distribution. 

**Alternative Hypothesis:**
If the Kimboto Trees are normally distributed, then there is a higher chance the data is spatially clustered. 

### 3.1 F(d) Function
We must use the F(d) and G(d) functions to determine whether our data follows the Poisson distribution. The distance at which at least 95% of the points have a neighbour is 19.57309 units (cornell). The proportion of points with a neighbour within 0.05 units is 0.007202148 (cornell).
Lastly, we use the plot function to demonstrate the theoretical value of F(r) and the reduced-sample estimate or “border” method of edge correction. 


```Ruby
#F-Function
F <- Fest(as(paracou, "ppp"))
class(F)
summary(F)
F$r[min(which(F$rs >= 0.95))]
F$rs[min(which(F$r >= 0.05))]
plot(F, cbind(rs,theo) ~ theo,
     main="F-function, Kimboto Trees")
```
![F-function](https://user-images.githubusercontent.com/118564598/229259917-f2247738-72fa-4ebe-863f-822de267083b.png)

### 3.2 G(d) Function
The G(d) function allows us to determine that the reduced sample estimate or “border” curve falls above the Poisson distribution curve and this implies that our dataset does not follow the Poisson distribution. Additionally, the distance at which at least 95% of the points have a neighbour is 16.14386 units (cornell). The proportion of points with a neighbour within 0.05 units is 0.002262443 (cornell). So far we have formed two subjective opinions about which patterns are consistent with CSR. However, our F-function and G-function tests are not consistent with each other and, therefore, we need to compute G- and F-function envelope tests. 


```Ruby
G <- Gest(as(paracou, "ppp"))
class(G)
summary(G)
G$r[min(which(G$rs >= 0.95))]
G$rs[min(which(G$r >= 0.05))]
plot(G, cbind(rs,theo) ~ theo,
     main="G-function, Kimboto Trees")
```
![G-function](https://user-images.githubusercontent.com/118564598/229259923-906d38f5-b9a6-42a9-9a1c-17c230a3f9d9.png)

### 3.3 F(d) Function Envelope
After completing the F and G function to formulate the visual confirmation whether our data follows the Poisson distribution, we can then use a test to assess our data against the F and G functions. This is done by utilizing an envelope test on the F and G functions to simulate whether our data fits into the F and G curves and compare it against our model. Once the envelopes are plotted, we can then see whether our data is contained inside, or if it is not, then where or at what point does it not. To do this, we calculate the following code structure for the F-test.

```Ruby
rmax.kimboto <- G$r[min(which(G$rs > 0.95))]
r <- seq(0, rmax.kimboto, by = 0.005)
envkimboto <- envelope(as(paracou, "ppp"), fun=Fest,
                       r=r, nrank=1, nsim=99)
plot(envkimboto, xlim=c(0, rmax.kimboto),
     main="Kimboto Trees, F-function Envelope")
```
![F-function Envelope](https://user-images.githubusercontent.com/118564598/229259982-2eba44cd-837b-4087-b14e-9da4e56ef9d2.png)

Since our observed F-function falls below the envelope, this implies that our data is mostly inhomogeneous and that we cannot fully assume CSR. However, we cannot confidently confirm this by running one test. So, to strengthen our case, we compute the same envelope test but with the G-function.

### 3.4 G(d) Function Envelope
The G-test envelope allows us to examine if there are any other types of pattern based on the G-function.

```Ruby
rmax.kimboto <- G$r[min(which(G$rs > 0.95))]
r <- seq(0, rmax.kimboto, by = 0.005)
envkimboto <- envelope(as(paracou, "ppp"), fun=Gest,
                       r=r, nrank=1, nsim=99)
plot(envkimboto, xlim=c(0, rmax.kimboto),
     main="Kimboto Trees, G-function Envelope")
```
![G-function Envelope](https://user-images.githubusercontent.com/118564598/229259989-2cc8ab5f-a803-4b05-afb2-bff52adf50f9.png)

 In this case, our observed G-function almost aligns with the envelope, and even fits within the envelope at around the 9 metre mark. This would be a correct observation based on the fact that fitting into the envelope would mean our data gets more and more dispersed as the distances increase. We can also see that the data has some form of clustering and a bit of dispersion as well. The G-test helps us by providing more concrete evidence that our data is not completely random, which means we can start formulating whether we would reject the null hypothesis or not. 

### 3.5 Likelihood Ratio Test
The anova.ppm function in spatstat package performs analysis of deviance for two or more fitted models with Poisson interaction terms. The Likelihood Ratio Test argument, provides a p value which helps interpreting the results for the most suitable hypothesis.

```Ruby
print(m.pois <- ppm(paracou, trend = ~1, interaction = NULL))
m.ts1 <- ppm(paracou, trend = ~polynom(x, y, 1),
             interaction = NULL)
m.ts2 <- ppm(paracou, trend = ~polynom(x, y, 2),
             interaction = NULL)

anova(m.ts2, m.ts1, m.pois, test="LRT")

```
The anova test gives us the following results in R: 
```
Analysis of Deviance Table

Model 1: ~x + y + I(x^2) + I(x * y) + I(y^2) 	 Poisson
Model 2: ~x + y 	 Poisson
Model 3: ~1 	 Poisson
  Npar Df Deviance  Pr(>Chi)    
1    6                          
2    3 -3  -44.159 1.396e-09 ***
3    1 -2   -2.739    0.2543    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

The value under Pr is a p-value, which is greater than 0.05 therefore it is statistically significant to say that the base model for this point pattern data favours the Null hypothesis. That is, the point pattern distribution is Inhomogeneous Poisson Distribution. 

To elaborate, an inhomogeneous Poisson process is a Poisson point process with a Poisson parameter set as some location-dependent function in the underlying space on which the Poisson process is defined. For the Alternate hypothesis, this was our presumption, we had considered the intensity as a function of the point location. And this can be easily understood by the Kernel Density Estimation, which accounts for the probability density distribution of the points in a given space.

### 3.6 Kernel Density Estimation
After completing a p-value test to confirm our suspicions, we can compute a final test to visually validate our assumptions. This can be done by producing a Kernel Density Estimation (KDE). The KDE considers the local intensity for our dataset. The following code was used to run the KDE on our dataset and also provides us with an output plot. 

```Ruby
K1 <- density(paracou) 
plot(K1, main=NULL, las=1)
contour(K1, add=TRUE)
```
![KernelPlot](https://user-images.githubusercontent.com/118564598/229258994-858dd8e2-ad39-4493-be4a-e5cef433414d.png)

The KDE plot also confirms our predictions that there is evidence of some spatial clustering, specifically centred on the right, where there are more observed points per squared metres. 

