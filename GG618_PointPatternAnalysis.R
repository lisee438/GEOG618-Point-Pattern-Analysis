install.packages("spatstat")
require(spatstat)
data("paracou")
paracou

density <- summary(paracou)$intensity

plot(paracou, cols=2:3, chars=c(16,3), main = "Locations of Kimboto Trees")

#F-Function
F <- Fest(as(paracou, "ppp"))
class(F)
summary(F)
F$r[min(which(F$rs >= 0.95))]
F$rs[min(which(F$r >= 0.05))]
plot(F, cbind(rs,theo) ~ theo,
     main="F-function, Kimboto Trees")

G <- Gest(as(paracou, "ppp"))
class(G)
summary(G)
G$r[min(which(G$rs >= 0.95))]
G$rs[min(which(G$r >= 0.05))]
plot(G, cbind(rs,theo) ~ theo,
     main="G-function, Kimboto Trees")

rmax.kimboto <- G$r[min(which(G$rs > 0.95))]
r <- seq(0, rmax.kimboto, by = 0.005)
envkimboto <- envelope(as(paracou, "ppp"), fun=Fest,
                       r=r, nrank=1, nsim=99)
plot(envkimboto, xlim=c(0, rmax.kimboto),
     main="Kimboto Trees, F-function Envelope")

rmax.kimboto <- G$r[min(which(G$rs > 0.95))]
r <- seq(0, rmax.kimboto, by = 0.005)
envkimboto <- envelope(as(paracou, "ppp"), fun=Gest,
                       r=r, nrank=1, nsim=99)
plot(envkimboto, xlim=c(0, rmax.kimboto),
     main="Kimboto Trees, G-function Envelope")

print(m.pois <- ppm(paracou, trend = ~1, interaction = NULL))
m.ts1 <- ppm(paracou, trend = ~polynom(x, y, 1),
             interaction = NULL)
m.ts2 <- ppm(paracou, trend = ~polynom(x, y, 2),
             interaction = NULL)

anova(m.ts2, m.ts1, m.pois, test="LRT")

K1 <- density(paracou) 
plot(K1, main=NULL, las=1)
contour(K1, add=TRUE)

