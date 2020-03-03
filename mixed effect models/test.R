
library(lme4)

# library(lattice)
# xyplot(incidence/size ~ period|herd, cbpp, type=c('g','p','l'),
#        layout=c(3,5), index.cond = function(x,y)max(y))
# (gm1 <- glmer(incidence/size ~ period + (1 | herd),
#               data = cbpp, family = binomial))



library(brms)

library(HSAUR2)
data(toenail)


gm2 <- glmer(outcome~treatment+(1|patientID),
               data=toenail,
               family=binomial)


mod = brms::brm(outcome~treatment+(1|patientID),
                data=toenail,
                family=bernoulli)
