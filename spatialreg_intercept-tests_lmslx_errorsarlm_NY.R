library(spatialreg)
library(spdep)
library(sf)

library(multcomp)

NY8 <- st_read(system.file("shapes/NY8_utm18.shp", package="spData"))
NY_nb <- poly2nb(NY8)
NY_lw <- nb2listw(NY_nb, style="minmax")
NY_lw_w <- nb2listw(NY_nb, style="W")





#### --------- Tests for SDEM --------- ####

### SDEM with unstandardized W and intercept in lagX
m1b <- spatialreg::errorsarlm(Z ~ PEXPOSURE, 
                              data=NY8, listw=NY_lw, Durbin=TRUE)
# summary(m1b)
# summary(spatialreg::impacts(m1b))

# Impact equals coef
all.equal(spatialreg::impacts(m1b)[[1]]$indirect, m1b$coefficients["lag.PEXPOSURE"], 
          tolerance = 1e-5, check.attributes = FALSE)



### SDEM with unstandardized W, intercept in lagX, and multiple covariates
m1c <- spatialreg::errorsarlm(Z ~ PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                              data=NY8, listw=NY_lw, Durbin=TRUE)
# summary(m1c)
# summary(spatialreg::impacts(m1c))

# Impacts equal coefs
m1c.imps <- spatialreg::impacts(m1c)
all.equal(c(m1c.imps[[1]]$indirect["PEXPOSURE"],
            m1c.imps[[1]]$indirect["PCTAGE65P"],
            m1c.imps[[1]]$indirect["PCTOWNHOME"]), 
          c(m1c$coefficients["lag.PEXPOSURE"],
            m1c$coefficients["lag.PCTAGE65P"],
            m1c$coefficients["lag.PCTOWNHOME"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc <- summary(glht(m1c, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m1c.imps[[1]]$total["PCTAGE65P"],
            summary(m1c.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)



### SDEM with standardized W, without intercept in lagX, and multiple covariates
m1d <- spatialreg::errorsarlm(Z ~ PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                              data=NY8, listw=NY_lw_w, Durbin=TRUE)
# summary(m1d)
# summary(spatialreg::impacts(m1d))

# Impacts equal coefs
m1d.imps <- spatialreg::impacts(m1d)
all.equal(c(m1d.imps[[1]]$indirect["PEXPOSURE"],
            m1d.imps[[1]]$indirect["PCTAGE65P"],
            m1d.imps[[1]]$indirect["PCTOWNHOME"]), 
          c(m1d$coefficients["lag.PEXPOSURE"],
            m1d$coefficients["lag.PCTAGE65P"],
            m1d$coefficients["lag.PCTOWNHOME"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc <- summary(glht(m1d, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m1d.imps[[1]]$total["PCTAGE65P"],
            summary(m1d.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)



### Test Durbin is formula, with standardized W
m1e <- spatialreg::errorsarlm(Z ~ PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                              data=NY8, listw=NY_lw_w, 
                              Durbin = ~ PEXPOSURE + PCTAGE65P)
# summary(m1e)
# summary(spatialreg::impacts(m1e))

# Impacts equal coefs
m1e.imps <- spatialreg::impacts(m1e)
all.equal(c(m1e.imps[[1]]$indirect["PEXPOSURE"],
            m1e.imps[[1]]$indirect["PCTAGE65P"]), 
          c(m1e$coefficients["lag.PEXPOSURE"],
            m1e$coefficients["lag.PCTAGE65P"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc <- summary(glht(m1e, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m1e.imps[[1]]$total["PCTAGE65P"],
            summary(m1e.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)


### Test Durbin is formula, with unstandardized W and intercept
# Not allowed, switches to Durbin = TRUE
m1f <- spatialreg::errorsarlm(Z ~ PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                         data=NY8, listw=NY_lw, 
                         Durbin = ~ PEXPOSURE + PCTAGE65P)
# summary(m1f)
# summary(spatialreg::impacts(m1f))
# 
# # Impacts equal coefs
# m1f.imps <- spatialreg::impacts(m1f)
# all.equal(c(m1f.imps[[1]]$indirect["PEXPOSURE"],
#             m1f.imps[[1]]$indirect["PCTAGE65P"]), 
#           c(m1f$coefficients["lag.PEXPOSURE"],
#             m1f$coefficients["lag.PCTAGE65P"]), 
#           tolerance = 1e-5, check.attributes = FALSE)
# 
# # Total impact is linear combination
# lc <- summary(glht(m1f, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
# all.equal(c(m1f.imps[[1]]$total["PCTAGE65P"],
#             summary(m1f.imps)$semat["PCTAGE65P", "Total"]), 
#           c(lc$test$coefficients,
#             lc$test$sigma), 
#           tolerance = 1e-5, check.attributes = FALSE)






#### --------- Repeat tests for SLX--------- ####

### SLX with unstandardized W and intercept in lagX
m2b <- spatialreg::lmSLX(Z ~ PEXPOSURE, 
                              data=NY8, listw=NY_lw, Durbin=TRUE)
# summary(m2b)
# summary(spatialreg::impacts(m2b))

# Impact equals coef
all.equal(spatialreg::impacts(m2b)[[1]]$indirect, m2b$coefficients["lag.PEXPOSURE"], 
          tolerance = 1e-5, check.attributes = FALSE)



### SLX with unstandardized W, intercept in lagX, and multiple covariates
m2c <- spatialreg::lmSLX(Z ~ PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                              data=NY8, listw=NY_lw, Durbin=TRUE)
# summary(m2c)
# summary(spatialreg::impacts(m2c))

# Impacts equal coefs
m2c.imps <- spatialreg::impacts(m2c)
all.equal(c(m2c.imps[[1]]$indirect["PEXPOSURE"],
            m2c.imps[[1]]$indirect["PCTAGE65P"],
            m2c.imps[[1]]$indirect["PCTOWNHOME"]), 
          c(m2c$coefficients["lag.PEXPOSURE"],
            m2c$coefficients["lag.PCTAGE65P"],
            m2c$coefficients["lag.PCTOWNHOME"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc <- summary(glht(m2c, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m2c.imps[[1]]$total["PCTAGE65P"],
            summary(m2c.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)



### SLX with standardized W, without intercept in lagX, and multiple covariates
m2d <- spatialreg::lmSLX(Z ~ PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                              data=NY8, listw=NY_lw_w, Durbin=TRUE)
# summary(m2d)
# summary(spatialreg::impacts(m2d))

# Impacts equal coefs
m2d.imps <- spatialreg::impacts(m2d)
all.equal(c(m2d.imps[[1]]$indirect["PEXPOSURE"],
            m2d.imps[[1]]$indirect["PCTAGE65P"],
            m2d.imps[[1]]$indirect["PCTOWNHOME"]), 
          c(m2d$coefficients["lag.PEXPOSURE"],
            m2d$coefficients["lag.PCTAGE65P"],
            m2d$coefficients["lag.PCTOWNHOME"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc1 <- summary(glht(m2d, linfct = c("PEXPOSURE + lag.PEXPOSURE = 0")))
lc2 <- summary(glht(m2d, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m2d.imps[[1]]$total["PEXPOSURE"],
            summary(m2d.imps)$semat["PEXPOSURE", "Total"],
            m2d.imps[[1]]$total["PCTAGE65P"],
            summary(m2d.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)



### Test Durbin is formula, with standardized W
m2e <- spatialreg::lmSLX(Z ~ PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                              data=NY8, listw=NY_lw_w, 
                              Durbin = ~ PEXPOSURE + PCTAGE65P)
# summary(m2e)
# summary(spatialreg::impacts(m2e))

# Impacts equal coefs
m2e.imps <- spatialreg::impacts(m2e)
all.equal(c(m2e.imps[[1]]$indirect["PEXPOSURE"],
            m2e.imps[[1]]$indirect["PCTAGE65P"]), 
          c(m2e$coefficients["lag.PEXPOSURE"],
            m2e$coefficients["lag.PCTAGE65P"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc1 <- summary(glht(m2e, linfct = c("PEXPOSURE + lag.PEXPOSURE = 0")))
lc2 <- summary(glht(m2e, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m2e.imps[[1]]$total["PEXPOSURE"],
            summary(m2e.imps)$semat["PEXPOSURE", "Total"],
            m2e.imps[[1]]$total["PCTAGE65P"],
            summary(m2e.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)



### Test Durbin is formula, with unstandardized W and intercept
m2f <- spatialreg::lmSLX(Z ~ PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                         data=NY8, listw=NY_lw, 
                         Durbin = ~ PEXPOSURE + PCTAGE65P)
# summary(m2f)
# summary(spatialreg::impacts(m2f))

# Impacts equal coefs
m2f.imps <- spatialreg::impacts(m2f)
all.equal(c(m2f.imps[[1]]$indirect["PEXPOSURE"],
            m2f.imps[[1]]$indirect["PCTAGE65P"]), 
          c(m2f$coefficients["lag.PEXPOSURE"],
            m2f$coefficients["lag.PCTAGE65P"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc1 <- summary(glht(m2f, linfct = c("PEXPOSURE + lag.PEXPOSURE = 0")))
lc2 <- summary(glht(m2f, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m2f.imps[[1]]$total["PEXPOSURE"],
            summary(m2f.imps)$semat["PEXPOSURE", "Total"],
            m2f.imps[[1]]$total["PCTAGE65P"],
            summary(m2f.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)





#### --------- Tests intercept in SLX formula --------- ####


### SLX with unstandardized W and intercept in lagX
m3b <- spatialreg::lmSLX(Z ~ -1 + PEXPOSURE, 
                         data=NY8, listw=NY_lw_w, Durbin=TRUE)
# summary(m3b)
# summary(spatialreg::impacts(m3b))

# Impact equals coef
m3b.imps <- spatialreg::impacts(m3b)
all.equal(m3b.imps[[1]]$indirect, m3b$coefficients["lag.PEXPOSURE"], 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc <- summary(glht(m3b, linfct = c("PEXPOSURE + lag.PEXPOSURE = 0")))
all.equal(c(m3b.imps[[1]]$total,
            summary(m3b.imps)$semat["PEXPOSURE", "Total"]), 
          c(lc$test$coefficients,
            lc$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)



### SLX with unstandardized W, no intercept in X, and multiple covariates
m3c <- spatialreg::lmSLX(Z ~ - 1 + PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                         data=NY8, listw=NY_lw, Durbin=TRUE)
# summary(m3c)
# summary(spatialreg::impacts(m3c))

# Impacts equal coefs
m3c.imps <- spatialreg::impacts(m3c)
all.equal(c(m3c.imps[[1]]$indirect["PEXPOSURE"],
            m3c.imps[[1]]$indirect["PCTAGE65P"],
            m3c.imps[[1]]$indirect["PCTOWNHOME"]), 
          c(m3c$coefficients["lag.PEXPOSURE"],
            m3c$coefficients["lag.PCTAGE65P"],
            m3c$coefficients["lag.PCTOWNHOME"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc1 <- summary(glht(m3c, linfct = c("PEXPOSURE + lag.PEXPOSURE = 0")))
lc2 <- summary(glht(m3c, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m3c.imps[[1]]$total["PEXPOSURE"],
            summary(m3c.imps)$semat["PEXPOSURE", "Total"],
            m3c.imps[[1]]$total["PCTAGE65P"],
            summary(m3c.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)



### Test Durbin is formula, with standardized W
m3e <- spatialreg::lmSLX(Z ~ - 1 + PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                         data=NY8, listw=NY_lw_w, 
                         Durbin = ~ PEXPOSURE + PCTAGE65P)
# summary(m3e)
# summary(spatialreg::impacts(m3e))

# Impacts equal coefs
m3e.imps <- spatialreg::impacts(m3e)
all.equal(c(m3e.imps[[1]]$indirect["PEXPOSURE"],
            m3e.imps[[1]]$indirect["PCTAGE65P"]), 
          c(m3e$coefficients["lag.PEXPOSURE"],
            m3e$coefficients["lag.PCTAGE65P"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc1 <- summary(glht(m3e, linfct = c("PEXPOSURE + lag.PEXPOSURE = 0")))
lc2 <- summary(glht(m3e, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m3e.imps[[1]]$total["PEXPOSURE"],
            summary(m3e.imps)$semat["PEXPOSURE", "Total"],
            m3e.imps[[1]]$total["PCTAGE65P"],
            summary(m3e.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)

# In this setup, the total impact equals the direct impact if omitted from Durbin=formula
all.equal(c(m3e.imps[[1]]$total["PCTOWNHOME"]), 
          c(m3e.imps[[1]]$direct["PCTOWNHOME"]), 
          tolerance = 1e-5, check.attributes = FALSE)



### Test Durbin is formula, with unstandardized W and no intercept, and -1 in Durbin formula
m3f <- spatialreg::lmSLX(Z ~ -1 + PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                         data=NY8, listw=NY_lw, 
                         Durbin = ~ - 1 + PEXPOSURE + PCTAGE65P)
# summary(m3f)
# summary(spatialreg::impacts(m3f))

# Impacts equal coefs
m3f.imps <- spatialreg::impacts(m3f)
all.equal(c(m3f.imps[[1]]$indirect["PEXPOSURE"],
            m3f.imps[[1]]$indirect["PCTAGE65P"]), 
          c(m3f$coefficients["lag.PEXPOSURE"],
            m3f$coefficients["lag.PCTAGE65P"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc1 <- summary(glht(m3f, linfct = c("PEXPOSURE + lag.PEXPOSURE = 0")))
lc2 <- summary(glht(m3f, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m3f.imps[[1]]$total["PEXPOSURE"],
            summary(m3f.imps)$semat["PEXPOSURE", "Total"],
            m3f.imps[[1]]$total["PCTAGE65P"],
            summary(m3f.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)



### Test Durbin is formula, with unstandardized W and no intercept, no intercept declaration in Durbin formula
m3fb <- spatialreg::lmSLX(Z ~ -1 + PEXPOSURE + PCTAGE65P+ POP8 + PCTOWNHOME, 
                         data=NY8, listw=NY_lw, 
                         Durbin = ~ PEXPOSURE + PCTAGE65P)
# summary(m3fb)
# summary(spatialreg::impacts(m3fb))

# Impacts equal coefs
m3fb.imps <- spatialreg::impacts(m3fb)
all.equal(c(m3fb.imps[[1]]$indirect["PEXPOSURE"],
            m3fb.imps[[1]]$indirect["PCTAGE65P"]), 
          c(m3fb$coefficients["lag.PEXPOSURE"],
            m3fb$coefficients["lag.PCTAGE65P"]), 
          tolerance = 1e-5, check.attributes = FALSE)

# Total impact is linear combination
lc1 <- summary(glht(m3fb, linfct = c("PEXPOSURE + lag.PEXPOSURE = 0")))
lc2 <- summary(glht(m3fb, linfct = c("PCTAGE65P + lag.PCTAGE65P = 0")))
all.equal(c(m3fb.imps[[1]]$total["PEXPOSURE"],
            summary(m3fb.imps)$semat["PEXPOSURE", "Total"],
            m3fb.imps[[1]]$total["PCTAGE65P"],
            summary(m3fb.imps)$semat["PCTAGE65P", "Total"]), 
          c(lc1$test$coefficients,
            lc1$test$sigma,
            lc2$test$coefficients,
            lc2$test$sigma), 
          tolerance = 1e-5, check.attributes = FALSE)







#### --------- To be ignored --------- ####
# 
# 
# 
# ### Test SDM
# m3 <- spatialreg::lagsarlm(Z ~ PEXPOSURE  + PCTAGE65P, data=NY8, listw=NY_lw, Durbin=TRUE)
# summary(m3)
# spatialreg::impacts(m3, listw=NY_lw)
# 
# 
# # Manually
# W <- listw2mat(NY_lw)
# M <- solve(diag(1, nrow = length(m3$residuals)) - m3$rho * W)
# b <- m3$coefficients[c("PEXPOSURE", "PCTAGE65P")]
# lb <- m3$coefficients[c("lag.PEXPOSURE", "lag.PCTAGE65P")]
# for(i in 1:2){
#   Z <-  M %*% (diag(b[i], length(m3$residuals)) + W * lb[i])
#   #dir imps
#   print(mean(diag(Z)))
#   #indir imps
#   diag(Z) <- 0
#   print(mean(colSums(Z)))
# }




