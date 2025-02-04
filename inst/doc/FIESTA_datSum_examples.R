## ----setup--------------------------------------------------------------------
library(FIESTA)
options(scipen = 6)

## ----echo=-1------------------------------------------------------------------
data.table::setDTthreads(2)

## ----ex1----------------------------------------------------------------------
sumdat1 <- datSumTree(tree = WYtree,
                      tsumvarlst = c("BA", "VOLCFNET"),
                      tfilter = "STATUSCD == 1")

## Returned list items
names(sumdat1)

## The first six rows of the summarized data table.
head(sumdat1$treedat)

## The summarized variable names
sumdat1$sumvars

## The query used to get data (use message to output in pretty format)
message(sumdat1$treeqry)

## ----ex2----------------------------------------------------------------------
sumdat2 <- 
  datSumTree(tree = WYtree,
             seed = WYseed,
             tsumvarlst = c("BA", "VOLCFNET", "TPA_UNADJ"),
             tsumvarnmlst = c("BA_LIVE", "VOLNET_LIVE", "COUNT"),
             bydomainlst = "SPCD",
             tderive = list(SDI = '(POWER(DIA / 10, 1.605)) * TPA_UNADJ'),
             woodland = "N",
             seedlings = "Y",
             tfilter = "STATUSCD == 1")

## Returned list items
names(sumdat2)

## The first six rows of the summarized data table.
head(sumdat2$treedat)

## The summarized variable names
sumdat2$sumvars

## The query used to get data (use message to output in pretty format)
message(sumdat2$treeqry)

## ----ex3----------------------------------------------------------------------
## First, find unique species in WYtree
spcdlst <- sort(unique(WYtree$SPCD))
## specify new class values for each unique species in WYtree
spcdlut <- data.frame(SPCD = spcdlst,
                      SPCDCL = c("C","W","W","C","C","C","W","C","C","C","C","H","H","W","H","H","H","H","H"))

## Next, find unique diameters in WYtree
dialst <- sort(unique(WYtree$DIA))
## specify break values to define new diameter class
diabrks <- c(0,20,40,80)

sumdat3 <- 
  datSumTree(tree = WYtree,
             seed = WYseed,
             tsumvarlst = c("BA", "VOLCFNET", "TPA_UNADJ"),
             tsumvarnmlst = c("BA_LIVE", "VOLNET_LIVE", "COUNT"),
             bydomainlst = c("SPCD", "DIA"),
             tderive = list(SDI = '(POWER(DIA / 10, 1.605)) * TPA_UNADJ'),
             domclassify = list(SPCD = spcdlut, DIA = diabrks),
             woodland = "N",
             seedlings = "Y",
             tfilter = "STATUSCD == 1")

## Returned list items
names(sumdat3)

## The first six rows of the summarized data table.
head(sumdat3$treedat)

## The summarized variable names
sumdat3$sumvars

## The query used to get data (use message to output in pretty format)
message(sumdat3$treeqry)

## ----ex4----------------------------------------------------------------------

sumdat4 <- 
  datSumTree(tree = WYtree,
             tderive = list(LIVE_BA = "SUM(power(DIA, 2) * 0.005454 * TPA_UNADJ * (CASE WHEN STATUSCD = 1 THEN 1 ELSE 0 END))",
                            DEAD_BA = "SUM(power(DIA, 2) * 0.005454 * TPA_UNADJ * (CASE WHEN STATUSCD = 2 THEN 1 ELSE 0 END))",
                            SDI = "SUM((POWER(DIA / 10, 1.605)) * TPA_UNADJ)",
                            QMD = "sqrt(SUM(power(DIA, 2) * 0.005454 * TPA_UNADJ) / (SUM(TPA_UNADJ) * 0.005454))",
                            MEAN_DIA = "AVG(DIA)",
                            MEDIAN_DIA = "MEDIAN(DIA)",
                            LIVELESS20 = "SUM(TPA_UNADJ * (CASE WHEN DIA < 10 THEN 1 ELSE 0 END))",
                            LIVE10to30 = "SUM(TPA_UNADJ * (CASE WHEN DIA > 10 AND DIA <= 30 THEN 1 ELSE 0 END))"))
                          
## Returned list items
names(sumdat4)

## The first six rows of the summarized data table.
head(sumdat4$treedat)

## The summarized variable names
sumdat4$sumvars

## The query used to get data (use message to output in pretty format)
message(sumdat4$treeqry)

