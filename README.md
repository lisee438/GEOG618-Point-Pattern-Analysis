# GEOG618 - Point Pattern Analysis
GEOG618 Point Patter Analysis Term Project. Formulated by Lisa Tang, Arghavan Zarandi, and Tirtha Harshilkumar Gajjar.

# Point Pattern Analysis
Point Pattern Analysis is a topic under spatial analysis that can support explanations on various point patterns. 
![github-centrography](https://user-images.githubusercontent.com/118564598/228592282-5f54fa85-f859-4f84-9510-6b906ea7198a.svg)

The following workshop follows the Intro to GIS and Spatial Analysis text, specifically referring to chapter 11, which can be found here: https://mgimond.github.io/Spatial/chp11_0.html 

## Necessary Tools and Programs to be used:
This workshop will be conducted in R-studio and the dataset to be used is the "spatstat" package available in R. An R file is included as well to provide access to the entire code. To prepare for this workshop, please have R readily downloaded and the "spatstat" packaged installed so that it can be easily loaded in class. 

```Ruby
install.packages("spatstat")
library("spatstat")
require(spatstat)
data("paracou")
class("paracou")
str("paracou")
summary("paracou")
```

## 1. Plot Data Points
First we will plot our data points to visualize our dataset and understand what is happening spatially:
```Ruby
#Plot of Kimboto Trees
plot(paracou, cols=2:3, chars=c(16,3), main = "Locations of Kimboto Trees")
```

## 2. Global Density
To understand the basic global pattern of our points, we can compute a simple function to determine the pattern's overall density. 
<img width="264" alt="Global-density" src="https://user-images.githubusercontent.com/118564598/228600097-fbc0f92c-8c8b-42c8-896d-7b7d7976a1ef.png">
```Ruby
lamb <- summary(paracou)$intensity
```
Running the "lamb" in the console will provide the estimated point density value of 0.004205303, which tells us that the average intensity is 0.004205303 points per square metre. 

## 3. Local Density
However, computing the global density does not provide us with much information on the local patterns. To understand the more complex point patterns that will be addressed in this dataset, we will go through an overview of the G- and F-functions to provide us with more answers.

**Null Hypothesis:**
If the Kimboto Trees are not normally distributed, then the data follows the Poisson distribution. 

**Alternative Hypothesis:**
If the Kimboto Trees are normally distributed, then there is a higher chance the data is spatially clustered. 

### 3.1 F-Function


```Ruby
F <- Fest(as(paracou, "ppp"))
class(F)
summary(F)
plot(F, main = "F-function, Kimboto Trees")
plot(F, xlim = c(0, 0.16), main = "F-function, Kimboto Trees")
1 - exp(-paracou$n * pi * (0.16)^2)
G$r[min(which(F$rs >= 0.95))]
G$rs[min(which(F$r >= 0.05))]
plot(F, cbind(rs,theo) ~ theo,
     main="F-function, Kimboto Trees")
```

### 3.2 G-Function


```Ruby
G <- Gest(as(paracou, "ppp"))
class(G)
summary(G)
plot(G, main = "G-function, Kimboto Trees")
plot(G, xlim = c(0, 0.16), main = "G-function, Kimboto Trees")
1 - exp(-paracou$n * pi * (0.16)^2)
G$r[min(which(G$rs >= 0.95))]
G$rs[min(which(G$r >= 0.05))]
plot(G, cbind(rs,theo) ~ theo,
     main="G-function, Kimboto Trees")
```

### 3.3 F-function Envelope

```Ruby
rmax.kimboto <- G$r[min(which(G$rs > 0.95))]
r <- seq(0, rmax.kimboto, by = 0.005)
envkimboto <- envelope(as(paracou, "ppp"), fun=Fest,
                       r=r, nrank=1, nsim=99)
plot(envkimboto, xlim=c(0, rmax.kimboto),
     main="Kimboto Trees, F-function Envelope")
```

### 3.4 G-function Envelope

```Ruby
rmax.kimboto <- G$r[min(which(G$rs > 0.95))]
r <- seq(0, rmax.kimboto, by = 0.005)
envkimboto <- envelope(as(paracou, "ppp"), fun=Gest,
                       r=r, nrank=1, nsim=99)
plot(envkimboto, xlim=c(0, rmax.kimboto),
     main="Kimboto Trees, G-function Envelope")
```

### 3.5 Likelihood Ratio Test

```Ruby
print(m.pois <- ppm(paracou, trend = ~1, interaction = NULL))
m.ts1 <- ppm(paracou, trend = ~polynom(x, y, 1),
             interaction = NULL)
m.ts2 <- ppm(paracou, trend = ~polynom(x, y, 2),
             interaction = NULL)

anova(m.ts2, m.ts1, m.pois, test="LRT")
```

### 3.6 Kernel Density Estimation

```Ruby
K1 <- density(paracou) 
plot(K1, main=NULL, las=1)
contour(K1, add=TRUE)
```
