## ----setup, include = F-------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(message = F, warning = F)

## ---- include=FALSE-----------------------------------------------------------
# Sets up output folding
hooks = knitr::knit_hooks$get()
hook_foldable = function(type) {
  force(type)
  function(x, options) {
    res = hooks[[type]](x, options)
    
    if (isFALSE(options[[paste0("fold.", type)]])) return(res)
    
    paste0(
      "<details><summary>", type, "</summary>\n\n",
      res,
      "\n\n</details>"
    )
  }
}
knitr::knit_hooks$set(
  output = hook_foldable("output"),
  plot = hook_foldable("plot")
)

## ---- results = 'asis', echo=FALSE--------------------------------------------
stratdat.lut <- data.frame(
  Variable = c("ESTN_UNIT", "STRATUMCD", "P1POINTCNT", "n.strata", "n.total", "ACRES", "strwt"),
  Description = c("Estimation unit", "Strata value", "Number of pixels by strata and estimation unit", 
                  "Number of plots in strata (and estimation unit)", "Number of plots for estimation unit",
                  "Total acres for estimation unit", "Summed proportions by strata and estimation unit"),
  stringsAsFactors = FALSE)

kable(stratdat.lut,
  format = "pandoc",   # default
  caption = "Description of variables in stratdat.",
  col.names = names(stratdat.lut),
  row.names = FALSE,
  align = c("l"),       # align = c("c", "c", "c", "r")
  # padding = 2         # inner spacing
)

pltdom.lut <- data.frame(
  Variable = c("ESTN_UNIT", "STRATUMCD", "plot_id", "category", "nbrpts.pltdom", "PtsPerPlot", "p.pltdom"),
  Description = c("Estimation unit", "Strata value", "Unique identifier for ICE plot", 
                  "Category (domain) for estimation", "Number of points by category (domain)",
                  "Number of points interpreted", "Proportion of plot by category"),
  stringsAsFactors = FALSE)

kable(pltdom.lut,
  format = "pandoc",   # default
  caption = "Description of variables in pltdom.",
  col.names = names(pltdom.lut),
  row.names = FALSE,
  align = c("l"),       # align = c("c", "c", "c", "r")
  # padding = 2         # inner spacing
)

## ---- results = 'asis', echo=FALSE--------------------------------------------

nonratio <- data.frame(
  Variable = c("phat", "phat.var", "phat.se", "phat.cv", "est", "est.var"), 
  Description = c("Estimated proportion of land", "Variance estimate of estimated proportion of land",
                  "Standard error of estimated proportion of land { sqrt(phat.var) }",
                  "Coefficient of variance of estimated proportion of land { phat.se/phat }",
                  "Estimated percent cover of land { phat*100 }",
                  "Variance of estimated percent cover of land { phat.var*100^2 }"),
  stringsAsFactors = FALSE)

ratio <- data.frame(
  Variable = c("phat.n", "phat.var.n", "phat.d", "phat.var.d", "covar", "rhat",
               "rhat.var", "rhat.se", "rhat.cv", "est", "est.var"),
  Description = c("Estimated proportion of land, for numerator",
                  "Variance of estimated proportion of land, for numerator",
                  "Estimated proportion of land, for denominator",
                  "Variance of estimated proportion of land, for denominator",
                  "Covariance of estimated proportion of numerator and denominator",
                  "Ratio of estimated proportions (numerator/denominator)",
                  "Variance of ratio of estimated proportions",
                  "Standard error of ratio of estimated proportions { rhat.se/rhat }",
                  "Coefficient of variation of ratio of estimated proportions { sqrt(rhat.var) }",
                  "Estimated percent cover of land { rhat*100 }",
                  "Variance of estimated percent cover of land { rhat.var*100^2 }"),
  stringsAsFactors = FALSE)

both <- data.frame(
  Variable = c("nbrpts", "ACRES", "est.se", "est.cv", "pse"),
  Description = 
    c("Number of points used in estimate",
      "Total acres for estimation unit (if tabtype='AREA')",
      "Standard error of estimated percent cover of land { sqrt(est.var) }",
      "Coefficient of variance of estimated percent cover of land { est.se/est }",
      "Percent sampling error of the estimated percent cover of land { est.cv*100 }"),
  stringsAsFactors = FALSE)
  
# gainloss <- data.frame(
#   Variable = c("gain.val", "loss.val", "gain.est", "gain.se", "loss.est", "loss.se", "diff.est", "diff.se"),
#   Description = 
#     c("Binary class for gain (Not-class to class). For each class, all other values are grouped to Not-class",
#       "Binary class for loss (Not-class to class). For each class, all other values are grouped to Not-class"),
#   "Estimated percent cover where the Class went from Not-class to Class", "Standard error of estimated gain",
#   "Estimated percent cover where the Class went from Class to Not-class",
#   "Standard error of estimated loss",
#   "Difference of estimated gain and estimate loss",
#   "Standard error of the difference of estimated gain and estimated loss")

all <- data.frame(
  Variable = c("CI99left", "CI99right", "CI95left", "CI95right", "CI68left", "CI68right"),
  Description = c("Left tail of 99% confidence interval for estimate { est - (2.58*est.se) }",
                  "Right tail of 99% confidence interval for estimate { est + (2.58*est.se) }", 
                  "Left tail of 95% confidence interval for estimate { est - (1.96*est.se) }",
                  "Right tail of 95% confidence interval for estimate { est + (1.96*est.se) }", 
                  "Left tail of 68% confidence interval for estimate { est - (0.97*est.se) }",
                  "Right tail of 68% confidence interval for estimate { est + (0.97*est.se) }"), 
  stringsAsFactors = FALSE)

kable(nonratio,
  format = "pandoc",   # default
  caption = "Description of variables in processing tables for nonratio estimates.",
  col.names = names(nonratio),
  row.names = FALSE,
  align = c("l"),       # align = c("c", "c", "c", "r")
  # padding = 2         # inner spacing
)

kable(ratio,
  format = "pandoc",   # default
  caption = "Description of variables in processing tables for ratio estimates.",
  col.names = names(ratio),
  row.names = FALSE,
  align = c("l"),       # align = c("c", "c", "c", "r")
  # padding = 2         # inner spacing
)

kable(both,
  format = "pandoc",   # default
  caption = "Description of variables in processing tables for both nonratio and ratio estimates.",
  col.names = names(both),
  row.names = FALSE,
  align = c("l"),       # align = c("c", "c", "c", "r")
  # padding = 2         # inner spacing
)

kable(all,
  format = "pandoc",   # default
  caption = "Description of variables in processing tables for all estimates.",
  col.names = names(all),
  row.names = FALSE,
  align = c("l"),       # align = c("c", "c", "c", "r")
  # padding = 2         # inner spacing
)

## ---- warning = F, message = F------------------------------------------------
library(FIESTA)

## -----------------------------------------------------------------------------
outfolder <- tempdir()

## -----------------------------------------------------------------------------
## Get external data file names
icepntfn <- system.file("extdata", "PB_data/icepnt_utco1135.csv", package = "FIESTA")
icepltfn <- system.file("extdata", "PB_data/icepltassgn_utco1135.csv", package = "FIESTA")
icepctcoverfn <- system.file("extdata", "PB_data/icepctcover_utco1135.csv", package = "FIESTA")
icechg_agfn <- system.file("extdata", "PB_data/chg_ag_LUT.csv", package = "FIESTA")
icecoverfn <- system.file("extdata", "PB_data/cover_LUT.csv", package = "FIESTA")
unitareafn <- system.file("extdata", "PB_data/unitarea_utco1135.csv", package = "FIESTA")
strlutfn <- system.file("extdata", "PB_data/strlut_utco1135.csv", package = "FIESTA")

icepnt <- read.csv(icepntfn)
iceplt <- read.csv(icepltfn)
icepctcover <- read.csv(icepctcoverfn)

icecover <- read.csv(icecoverfn)
icechg_ag <- read.csv(icechg_agfn)

## -----------------------------------------------------------------------------
str(icepnt, max.level = 1)

## -----------------------------------------------------------------------------
str(iceplt, max.level = 1)

## -----------------------------------------------------------------------------
str(icepctcover, max.level = 1)

## -----------------------------------------------------------------------------
icepltsp <- spMakeSpatialPoints(xyplt = iceplt,
                                xy.uniqueid = "plot_id",
                                xvar = "LON_PUBLIC",
                                yvar = "LAT_PUBLIC", 
                                prj = "longlat",
                                datum = "NAD83")
plot(icepltsp["ESTN_UNIT"])

## -----------------------------------------------------------------------------
icecover

## -----------------------------------------------------------------------------
icechg_ag

## -----------------------------------------------------------------------------
# Create look-up tables for Time 1 (cover_11) and Time 2 (cover_14) classes
icecover_1 <- icecover
names(icecover_1) <- sub("cover", "cover_1", names(icecover_1))
icecover_2 <- icecover
names(icecover_2) <- sub("cover", "cover_2", names(icecover_2))
icecover_1
icecover_2

## -----------------------------------------------------------------------------
## Area by estimation unit
unitarea <- read.csv(unitareafn)
unitarea

## -----------------------------------------------------------------------------
## Pixel counts by strata classes
strlut <- read.csv(strlutfn)
strlut

## -----------------------------------------------------------------------------
# Percent land cover at Time 1 (2011) for all land in Davis and Salt Lake Counties, UT
PBpopdat <- modPBpop(pnt = icepnt, 
                     pltassgn = iceplt,
                     pltassgnid = "plot_id",
                     pntid = "dot_cnt")
names(PBpopdat)


## -----------------------------------------------------------------------------
str(PBpopdat, max.level = 1)

## -----------------------------------------------------------------------------
# read in file 
unitarea <- read.csv(unitareafn)

# sum up the acres
sum(unitarea$ACRES)

## -----------------------------------------------------------------------------
PBpoparea <- modPBpop(pnt = icepnt, 
                      pltassgn = iceplt, 
                      pltassgnid = "plot_id", 
                      pntid = "dot_cnt", 
                      unitarea = sum(unitarea$ACRES)) # using the total number of acres

## -----------------------------------------------------------------------------
str(PBpoparea, max.level = 1)

## -----------------------------------------------------------------------------
PBpopunit <- modPBpop(pnt = icepnt, 
                      pltassgn = iceplt, 
                      pltassgnid = "plot_id", 
                      pntid = "dot_cnt",
                      unitarea = unitarea, 
                      unitvar = "ESTN_UNIT")
names(PBpopunit)

## -----------------------------------------------------------------------------
str(PBpopunit, max.level = 1)

## -----------------------------------------------------------------------------
head(icepctcover)
dim(icepctcover)

## -----------------------------------------------------------------------------
names11 <- names(icepctcover)[endsWith(names(icepctcover), "11")]
names14 <- names(icepctcover)[endsWith(names(icepctcover), "14")]
names11
names14

## -----------------------------------------------------------------------------
PBpctpop11 <- modPBpop(pltpct = icepctcover, 
                       pltpctvars = names11,
                       unitarea = sum(unitarea$ACRES))

## -----------------------------------------------------------------------------
str(PBpctpop11, max.level = 1)

## -----------------------------------------------------------------------------
PBpctpop14 <- modPBpop(pltpct = icepctcover, 
                       pltpctvars = names14,
                       unitarea = sum(unitarea$ACRES))

## -----------------------------------------------------------------------------
PBpctpop.veg <- modPBpop(pltpct = icepctcover, 
                       pltpctvars = "Veg.NonVeg",
                       unitarea = sum(unitarea$ACRES)
                       )

## -----------------------------------------------------------------------------
## Plot-level assignments
head(iceplt)

## Strata weights by estimation unit
head(read.csv(strlutfn))

## -----------------------------------------------------------------------------
PBpopareaPS <- modPBpop(pntdat = icepnt, 
                        pltassgn = iceplt, 
                        pltassgnid = "plot_id", 
                        pntid = "dot_cnt", 
                        strata = TRUE,
                        stratalut = strlutfn,
                        strvar = "STRATUMCD",
                        strata_opts = list(getwt=TRUE, 
                                           getwtvar="P1POINTCNT"),
                        unitarea = sum(unitarea$ACRES))

## -----------------------------------------------------------------------------
str(PBpopareaPS, max.level = 1)

## -----------------------------------------------------------------------------
PBpopareaPS$stratalut

## -----------------------------------------------------------------------------
PBpoparea_nonPS <- modPBpop(pntdat = icepnt, 
                        pltassgn = iceplt, 
                        pltassgnid = "plot_id", 
                        pntid = "dot_cnt", 
                        strata = FALSE,
                        unitarea = sum(unitarea$ACRES))

## -----------------------------------------------------------------------------
cover1 <- modPB(PBpopdat = PBpopdat, 
                rowvar = "cover_1", 
                table_opts = list(rowlut = icecover_1, 
                                  row.add0 = TRUE), 
                title_opts = list(title.rowvar = "Land Cover (2011)"))

## -----------------------------------------------------------------------------
str(cover1, max.level = 2)

## -----------------------------------------------------------------------------
str(cover1$est, max.level = 2)

## -----------------------------------------------------------------------------
cover1$raw$unit_rowest

## -----------------------------------------------------------------------------
head(cover1$raw$pltdom.row)

## -----------------------------------------------------------------------------
cover1 <- modPB(PBpopdat = PBpopdat,
		            rowvar = "cover_1", 
		            nonsamp.pntfilter = "cover_1 != 999", # added filter 
		            table_opts = list(rowlut = icecover_1), 
		            title_opts = list(title.rowvar = "Land Cover (2011)"),
		            returntitle = TRUE)

## -----------------------------------------------------------------------------
cover1$est

## -----------------------------------------------------------------------------
cover1.area <- modPB(PBpopdat = PBpoparea, 
                     tabtype = "AREA",
                     rowvar = "cover_1", 
                     nonsamp.pntfilter = "cover_1 != 999",
                     table_opts = list(rowlut = icecover_1), 
                     title_opts = list(title.rowvar = "Land Cover (2011)"))

## -----------------------------------------------------------------------------
str(cover1.area, max.level = 1)

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
cover1.area$est

## -----------------------------------------------------------------------------
cover1.pct <- modPB(PBpopdat = PBpoparea, 
                tabtype = "PCT", 
                rowvar = "cover_1", 
                nonsamp.pntfilter = "cover_1 != 999",
                table_opts = list(rowlut = icecover_1), 
                title_opts = list(title.rowvar = "Land Cover (2011)"),
                returntitle = TRUE, 
                savedata = TRUE, 
                savedata_opts = list(outfolder = outfolder))

## -----------------------------------------------------------------------------
str(cover1.pct, max.level = 1)

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
cover1.pct$est

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
cover1.pct$titlelst

## -----------------------------------------------------------------------------
cover2 <- modPB(PBpopdat = PBpoparea, 
                rowvar = "cover_2", 
                nonsamp.pntfilter = "cover_1 != 999",
                table_opts = list(rowlut = icecover_2), 
                title_opts = list(title.rowvar = "Land Cover (2014)"),
		            returntitle = TRUE)

## -----------------------------------------------------------------------------
str(cover2, max.level = 1)

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
cover2$est

## -----------------------------------------------------------------------------
netchg <- data.frame(Estimate1 = cover1$raw$unit_rowest$est, 
                     Estimate2 = cover2$raw$unit_rowest$est, 
                     NetChange.1to2 = cover1$raw$unit_rowest$est - cover2$raw$unit_rowest$est)
netchg

## -----------------------------------------------------------------------------
tabvars <- c("est", "est.se")
tab1 <- cover1$raw$unit_rowest[, c("cover_1", cover1$titlelst$title.rowvar, tabvars)]
data.table::setnames(tab1, tabvars, paste0(tabvars, ".1"))

tab2 <- cover2$raw$unit_rowest[, c("cover_2", cover2$titlelst$title.rowvar, tabvars)]
data.table::setnames(tab2, tabvars, paste0(tabvars, ".2"))

tabx <- merge(tab1, tab2, by.x="cover_1", by.y="cover_2")
tabx

## -----------------------------------------------------------------------------
sevar <- names(tabx)[grepl("est.se", names(tabx))]
yvar <- names(tabx)[grepl("est.", names(tabx)) & !names(tabx) %in% sevar]
xvar <- cover1$titlelst$title.rowvar

datBarplot(tabx, 
           yvar = yvar, 
           xvar = xvar,  
           errbars = TRUE,
           sevar = sevar, 
           ylabel = "Percent", 
           addlegend = TRUE, 
           args.legend = list(x = "topleft", 
                              bty = "n", 
                              cex = .8, 
                              legend = c("2011", "2014")), 
           main = substr(cover1$titlelst$title.row,
                         1,
                         nchar(cover1$titlelst$title.row)-7))

## -----------------------------------------------------------------------------
chg_ag <- modPB(PBpopdat = PBpoparea, 
                rowvar = "chg_ag_2", 
                table_opts = list(rowlut = icechg_ag), 
                title_opts = list(title.rowvar = "Agent of Change"),
		            returntitle=TRUE)

## -----------------------------------------------------------------------------
str(chg_ag, max.level = 1)

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
chg_ag$est

## -----------------------------------------------------------------------------
chg_ag.area <- modPB(PBpopdat = PBpoparea, 
                    tabtype = "AREA",
                    rowvar = "chg_ag_2", 
                    table_opts = list(rowlut = icechg_ag, metric=TRUE), 
                    title_opts = list(title.rowvar = "Agent of Change"),
		                returntitle=TRUE)

## -----------------------------------------------------------------------------
str(chg_ag.area, max.level = 1)

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
chg_ag.area$est

## -----------------------------------------------------------------------------
chg_ag.area$raw$areaunits

## -----------------------------------------------------------------------------
# Add a landarea filter to subset dataset to only plots with observed change.
landarea.filter <- "change_1_2 == 1"

chg_ag.plts <- modPB(PBpopdat = PBpoparea, 
                     rowvar = "chg_ag_2", 
                     table_opts = list(rowlut = icechg_ag), 
                     title_opts = list(title.rowvar = "Agent of Change"),
		                 landarea = "CHANGE", 
		                 landarea.filter = landarea.filter, 
		                 returntitle = TRUE)

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
chg_ag.plts$est

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
chg_ag.plts$titlelst

## -----------------------------------------------------------------------------
# Percent land changed by agent of change in Davis and Salt Lake Counties, UT
pntfilter <- "chg_ag_2 > 0"
chg_ag.pnts <- modPB(PBpopdat = PBpoparea, 
                     rowvar = "chg_ag_2", 
                     table_opts = list(rowlut = icechg_ag), 
                     title_opts = list(title.rowvar = "Agent of Change", 
                                       title.filter = "observed changed"),
		                 pntfilter = pntfilter, 
		                 returntitle = TRUE)

## -----------------------------------------------------------------------------
# All land
chg_ag$titlelst$title.estpse
chg_ag$est

# Land with observed change
chg_ag.plts$titlelst$title.estpse
chg_ag.plts$est

# Estimated change
chg_ag.pnts$titlelst$title.estpse
chg_ag.pnts$est

## -----------------------------------------------------------------------------
datBarplot(chg_ag.pnts$raw$unit_rowest, 
           xvar = "Agent of Change", 
           yvar = "est", 
           errbars = TRUE, 
           sevar = "est.se", 
           ylab = "Percent", 
           main = chg_ag.pnts$titlelst$title.row)

## -----------------------------------------------------------------------------
chg_ag_cover1 <- modPB(PBpopdat = PBpoparea, 
                       rowvar = "chg_ag_2", 
                       colvar = "cover_2", 
                       table_opts = list(rowlut = icechg_ag,
                                         collut = icecover_2), 
                       title_opts = list(title.rowvar = "Change agent",
                                         title.colvar = "Land cover (2011)"), 
                       returntitle = TRUE)

## -----------------------------------------------------------------------------
chg_ag_cover1$est

## -----------------------------------------------------------------------------
chg_ag_cover1$pse

## -----------------------------------------------------------------------------
cover1.unit.area <- modPB(PBpopdat = PBpopunit, 
                          tabtype = "AREA",
                          rowvar = "cover_1", 
                          nonsamp.pntfilter = "cover_1 != 999",
                          table_opts = list(rowlut=icecover_1), 
                          title_opts = list(title.rowvar="Land Cover (2011)"))

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
cover1.unit.area$est

## -----------------------------------------------------------------------------
## Percent sampling error of estimate
cover1.unit.area$pse

## -----------------------------------------------------------------------------
cover1.unitsum <- modPB(PBpopdat = PBpopunit, 
                        tabtype = "AREA",
                        sumunits = TRUE,
                        rowvar = "cover_1", 
                        nonsamp.pntfilter = "cover_1 != 999",
                        table_opts = list(rowlut=icecover_1), 
                        title_opts = list(title.rowvar="Land Cover (2011)"))

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
cover1.unitsum$est

## -----------------------------------------------------------------------------
str(cover1.unitsum$raw, max.level = 1)

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
cover1.unitsum$raw$unit_rowest

## -----------------------------------------------------------------------------
## Estimate and percent sampling error of estimate
cover1.unitsum$raw$rowest

## -----------------------------------------------------------------------------
cover12 <- modPB(PBpopdat = PBpoparea, 
                  rowvar = "cover_1", 
                  colvar = "cover_2", 
                  nonsamp.pntfilter = "cover_1 != 999",
                  table_opts = list(rowlut = icecover_1,
                                    collut = icecover_2), 
              		title_opts = list(title.rowvar = "Land Cover (2011)", 
              		                  title.colvar = "Land Cover (2014)"), 
    		          returntitle = TRUE)

## -----------------------------------------------------------------------------
cover12$est

## -----------------------------------------------------------------------------
cover12$pse

## -----------------------------------------------------------------------------
head(cover12$raw$pltdom.grp)

## -----------------------------------------------------------------------------
head(cover12$raw$unit_grpest)

## -----------------------------------------------------------------------------
cover12$raw$unit.grpest

## -----------------------------------------------------------------------------
cover12.area <- modPB(PBpopdat = PBpoparea, 
                 tabtype = "AREA",
                 rowvar = "cover_1", 
                 colvar = "cover_2", 
                 nonsamp.pntfilter="cover_1 != 999",
                 table_opts = list(rowlut = icecover_1,
                                   collut = icecover_2), 
              	 title_opts = list(title.rowvar = "Land Cover (2011)", 
              		                  title.colvar = "Land Cover (2014)"), 
    		         returntitle = TRUE)

## -----------------------------------------------------------------------------
head(cover12$pse)
head(cover12.area$pse)

## -----------------------------------------------------------------------------
PBpoparea2 <- PBpoparea
PBpoparea2$PBx <- merge(PBpoparea2$PBx, icecover_1, by = "cover_1")
PBpoparea2$PBx <- merge(PBpoparea2$PBx, icecover_2, by = "cover_2")
PBpoparea2$PBx$cover_12_nm <- paste(PBpoparea2$PBx$cover_1_nm, 
                                    PBpoparea2$PBx$cover_2_nm,
                                    sep = "-")
head(PBpoparea2$PBx)

## -----------------------------------------------------------------------------
cover12nm <- modPB(PBpopdat = PBpoparea2,
                   rowvar = "cover_12_nm", 
                   nonsamp.pntfilter = "cover_1 != 999", 
                   title_opts = list(title.rowvar = "Land Cover (2011-2014)"),  
		               returntitle = TRUE)

## -----------------------------------------------------------------------------
cover12nm$est
cover12$est

## -----------------------------------------------------------------------------
cover12.lt200 <- modPB(PBpopdat = PBpoparea,
                       rowvar = "cover_1", 
                       colvar = "cover_2", 
                       nonsamp.pntfilter = "cover_1 != 999", 
                       pntfilter = "cover_1 < 200",
                       table_opts = list(rowlut = icecover_1,
                                         collut = icecover_2), 
                    	 title_opts = list(title.rowvar = "Land Cover (2011)", 
                    	                   title.colvar = "Land Cover (2014)",
                    	                   title.filter = "Vegetated land"), 
                    	 returntitle = TRUE)


## -----------------------------------------------------------------------------
cover12.lt200$est

## -----------------------------------------------------------------------------
cover12.lt200$titlelst

## -----------------------------------------------------------------------------

cover12b <- modPB(PBpopdat = PBpoparea, 
                  rowvar = "cover_1", 
                  colvar = "cover_2", 
                  nonsamp.pntfilter="cover_1 != 999",
                  table_opts = list(rowlut = icecover_1,
                                    collut = icecover_2), 
              		title_opts = list(title.rowvar = "Land Cover (2011)", 
              		                  title.colvar = "Land Cover (2014"), 
    		          returntitle = TRUE,
    		          gainloss = TRUE)

## -----------------------------------------------------------------------------
str(cover12b$raw, max.level = 1)

## -----------------------------------------------------------------------------
cover12b$raw$est.gainloss

## -----------------------------------------------------------------------------
datPBplotchg(cover12b$raw$est.gainloss)

## -----------------------------------------------------------------------------
## We will first subset the raw data frame and set to an object
estcat <- "OtherVegetation"
othveg.gainloss <- cover12b$raw$est.gainloss[row.names(cover12b$raw$est.gainloss) == estcat,]

## -----------------------------------------------------------------------------
othveg.gainloss[, c("gain.CI95left", "gain.est", "gain.CI95right")]

## -----------------------------------------------------------------------------
othveg.gainloss[, c("loss.CI95left", "loss.est", "loss.CI95right")]

## -----------------------------------------------------------------------------
othveg.gainloss[, c("diff.CI95left", "diff.est", "diff.CI95right")]

## -----------------------------------------------------------------------------
changelut <- data.frame(change_1_2=c(0,1,2), 
                        change_1_2nm=c("No Change", "Change", "Expected Change"))
changelut

## -----------------------------------------------------------------------------
chgcover1 <- modPB(PBpopdat = PBpoparea, 
                   ratio = TRUE, 
                   rowvar = "change_1_2", 
                   colvar = "cover_1",
                   nonsamp.pntfilter = "cover_1 != 999",
                   table_opts = list(rowlut=changelut, collut=icecover_1),  
                   title_opts = list(title.rowvar="Change"))

## -----------------------------------------------------------------------------
chgcover1$est

## -----------------------------------------------------------------------------
chgcover1$pse

## -----------------------------------------------------------------------------
sum(as.numeric(chgcover1$est[1,-1]))
sum(as.numeric(chgcover1$est[2,-1]))

## -----------------------------------------------------------------------------

chg_ag_cover1.rat <- modPB(PBpopdat = PBpoparea, 
                           ratio = TRUE, 
                           rowvar = "chg_ag_2", 
                           colvar = "cover_1", 
                           nonsamp.pntfilter = "cover_1 != 999", 
                           table_opts = list(rowlut = icechg_ag, 
                                             collut = icecover_1), 
                        	 title_opts = list(title.rowvar = "Change agent", 
		                                         title.colvar = "Land cover (2011)"), 
                           returntitle = TRUE)

## -----------------------------------------------------------------------------
chg_ag_cover1.rat$est

## -----------------------------------------------------------------------------
chg_ag_cover1.rat$pse

## -----------------------------------------------------------------------------
chg_ag_cover1.rat$est$Total <- rowSums(apply(chg_ag_cover1.rat$est[,-1], 2, as.numeric),
                                       na.rm = TRUE)
chg_ag_cover1.rat$est

## -----------------------------------------------------------------------------
# Nonratio estimates
chg_ag_cover1$est

# Ratio to means estimates
chg_ag_cover1.rat$est

## -----------------------------------------------------------------------------
cover1_2.rat <- modPB(PBpopdat = PBpoparea,
		                  ratio = TRUE, 
		                  rowvar = "cover_1", 
		                  colvar = "cover_2", 
		                  nonsamp.pntfilter = "cover_1 != 999",
		                  table_opts = list(rowlut = icecover_1,
		                                    collut = icecover_2), 
                    	title_opts = list(title.rowvar = "Land cover (2011)", 
		                                    title.colvar = "Land cover (2014)"), 
		                  returntitle=TRUE)

## -----------------------------------------------------------------------------
cover1_2.rat$est

## -----------------------------------------------------------------------------
datBarStacked(x = cover1_2.rat$raw$unit_grpest, 
              main.attribute = "Land cover (2011)", 
              sub.attribute = "Land cover (2014)",
              response = "est", 
              xlabel = "Land Cover (2011)", 
              legend.title = "Land Cover (2014)")

## -----------------------------------------------------------------------------
x <- cover1_2.rat$raw$unit_grpest
x <- x[x$'Land cover (2011)' != x$'Land cover (2014)',]

datBarStacked(x = x, 
              main.attribute = "Land cover (2011)", 
              sub.attribute = "Land cover (2014)",
              response = "est", 
              xlabel = "Land Cover (2011)", 
              legend.title = "Land Cover (2014)",
              main.order = rev(c("Tree", "Shrub", "OtherVegetation",
                                 "Impervious", "Barren", "Water")))

## -----------------------------------------------------------------------------
pltpct11 <- modPB(PBpopdat = PBpctpop11, 
                  title_opts = list(title.rowvar="Land cover (2011)"),
                  returntitle = TRUE)
pltpct11$est

## -----------------------------------------------------------------------------
datBarplot(x = pltpct11$est, 
           xvar = "Land cover (2011)",
           yvar = "Estimate", 
           errbars = TRUE, 
           psevar = "Percent Sampling Error")

## -----------------------------------------------------------------------------
datBarplot(x = pltpct11$raw$unit_rowest, 
           xvar = "Land cover (2011)", 
           yvar = "est", 
           errbars = TRUE, 
           sevar = "est.se", 
           ylim = c(0,30), 
           ylabel = "Percent of land", 
           main = "Percent cover at Time 1 (2011)")

## -----------------------------------------------------------------------------
pltpct11.area <- modPB(PBpopdat = PBpctpop11, 
                       tabtype = "AREA",
                       returntitle = TRUE)
pltpct11.area$est

## -----------------------------------------------------------------------------
pltpct14 <- modPB(PBpopdat = PBpctpop14, 
                  returntitle = TRUE)
pltpct14$est

## -----------------------------------------------------------------------------
pltpct14.area <- modPB(PBpopdat = PBpctpop14, 
                       tabtype = "AREA",
                       returntitle = TRUE)
pltpct14.area$est

## -----------------------------------------------------------------------------
pltpct.veg <- modPB(PBpopdat = PBpctpop.veg, 
                    title_opts = list(title.rowvar = "Veg to NonVeg transition"),
                    returntitle = TRUE)
pltpct.veg$est

## -----------------------------------------------------------------------------
cover12ps <- modPB(PBpopdat = PBpopareaPS,
		               rowvar = "cover_1", 
		               colvar = "cover_2", 
		               nonsamp.pntfilter = "cover_1 != 999",
		               table_opts = list(rowlut = icecover_1,
		                                 collut = icecover_2), 
		               title_opts = list(title.rowvar = "Land Cover"))

## -----------------------------------------------------------------------------
cover12 <- modPB(PBpopdat = PBpoparea_nonPS,
		              rowvar = "cover_1", 
		              colvar = "cover_2", 
                  nonsamp.pntfilter = "cover_1 != 999",
		              table_opts = list(rowlut = icecover_1,
		                                collut = icecover_2), 
		              title_opts = list(title.rowvar = "Land Cover"))

## -----------------------------------------------------------------------------
cover12$est
cover12ps$est

cover12$pse
cover12ps$pse

