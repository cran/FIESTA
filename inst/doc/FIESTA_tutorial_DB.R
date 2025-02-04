## ----setup, include = F-------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(message = F, warning = F, eval = F)

## ----include=FALSE------------------------------------------------------------
# # Sets up output folding
# hooks = knitr::knit_hooks$get()
# hook_foldable = function(type) {
#   force(type)
#   function(x, options) {
#     res = hooks[[type]](x, options)
# 
#     if (isFALSE(options[[paste0("fold.", type)]])) return(res)
# 
#     paste0(
#       "<details><summary>", type, "</summary>\n\n",
#       res,
#       "\n\n</details>"
#     )
#   }
# }
# knitr::knit_hooks$set(
#   output = hook_foldable("output"),
#   plot = hook_foldable("plot")
# )

## ----echo=-1------------------------------------------------------------------
# data.table::setDTthreads(2)

## ----warning = F, message = F-------------------------------------------------
# library(FIESTA)

## -----------------------------------------------------------------------------
# outfolder <- tempdir()

## -----------------------------------------------------------------------------
# ## Get plot table for Wyoming
# WYplots <- DBgetCSV("PLOT", "Wyoming")
# dim(WYplots)
# 
# ## Get plot table for Wyoming and Utah
# WYUTplots <- DBgetCSV(DBtable = "PLOT",
#                       states = c("Wyoming", "Utah"))
# table(WYUTplots$STATECD)
# 
# ## Get survey table for Wyoming
# WYsurvey <- DBgetCSV("SURVEY", "Wyoming")
# WYsurvey

## -----------------------------------------------------------------------------
# # Get number of plots by inventory year for the state of Wyoming
# sql1 <- "SELECT INVYR, COUNT(*) AS NBRPLOTS
#          FROM PLOT
#          WHERE statecd = 56
#          GROUP BY INVYR"
# 
# nplots1 <- DBqryCSV(sql = sql1,
#                     states = "Wyoming",
#                     sqltables = "PLOT")
# 
# head(nplots1)
# 
# # Get number of plots by inventory year for Vermont and New Hampshire
# sql2 <- "SELECT STATECD, INVYR, COUNT(*) NBRPLOTS
#          FROM PLOT
#          WHERE statecd IN(50,33)
#          GROUP BY STATECD, INVYR"
# 
# nplots2 <- DBqryCSV(sql = sql2,
#                     states = c("Vermont", "New Hampshire"),
#                     sqltables = "PLOT")
# 
# head(nplots2)
# 
# # Get number of plots by inventory year for Iowa (stcd=19) that have silver maple (SPCD=317)
# sql3 <- "SELECT p.STATECD, p.INVYR, COUNT(*) NBRPLOTS
#          FROM PLOT p
#          JOIN TREE t ON p.CN = t.PLT_CN
#          WHERE p.statecd = 19 AND t.SPCD = 317
#          GROUP BY p.STATECD, p.INVYR"
# 
# nplots3 <- DBqryCSV(sql = sql3,
#                     states = "IOWA",
#                     sqltables = c("PLOT", "TREE"))
# 
# head(nplots3)

## -----------------------------------------------------------------------------
# WYeval <- DBgetEvalid(states = "Wyoming",
#                       evalCur = TRUE)
# 
# names(WYeval)
# WYeval$evalidlist
# WYeval$invyrs
# WYeval$invyrtab
# WYeval$invtype

## -----------------------------------------------------------------------------
# NYeval <- DBgetEvalid(states = c("New York"),
#                       evalType = c("VOL", "GRM"),
#                       evalCur = TRUE)
# 
# names(NYeval)
# NYeval$evalidlist
# NYeval$evalTypelist

## -----------------------------------------------------------------------------
# xydat1 <- DBgetXY(states = "Wyoming",
#                   datsource = "datamart",
#                   eval = "FIA",
#                   eval_opts = eval_options(Cur = TRUE))
# 
# names(xydat1)
# head(xydat1$xyCur_PUBLIC)

## -----------------------------------------------------------------------------
# xydat2 <- DBgetXY(states = "Wyoming",
#                   datsource = "datamart",
#                   eval = "FIA",
#                   eval_opts = eval_options(Cur = TRUE),
#                   pvars2keep = c("PLOT_STATUS_CD"),
#                   issp = TRUE)
# 
# spxy2 <- xydat2$spxy
# 
# ## Display points with by PLOT_STATUS_CD (1-light blue; 2-brown; 3-blue)
# spxy2$color <- ifelse(spxy2$PLOT_STATUS_CD == 2, "brown",
#                       ifelse(spxy2$PLOT_STATUS_CD == 3, "blue", "light blue"))
# 
# plot(sf::st_geometry(spxy2['PLOT_STATUS_CD']), pch = 16, cex = .5,
#                      col = spxy2$color)

## -----------------------------------------------------------------------------
# xydat3 <- DBgetXY(states = "Vermont",
#                   datsource = "datamart",
#                   eval = "custom",
#                   eval_opts = eval_options(invyrs = 2017:2019),
#                   issp = TRUE)
# 
# spxy3 <- xydat3$spxy
# 
# ## Display points
# plot(sf::st_geometry(spxy3), pch = 16, cex = .5, col="grey")
# 
# ## Now only include P2 plots only (intensity1 = TRUE)
# xydat3b <- DBgetXY(states = "Vermont",
#                    datsource = "datamart",
#                    eval = "custom",
#                    eval_opts = eval_options(invyrs = 2017:2019),
#                    intensity1 = TRUE,
#                    issp = TRUE)
# 
# spxy3b <- xydat3b$spxy
# 
# ## Display points
# plot(sf::st_geometry(spxy3b), pch = 16, cex = .5)

## -----------------------------------------------------------------------------
# dat1 <- DBgetPlots(states = "Rhode Island",
#                    datsource = "datamart",
#                    eval = "FIA",
#                    eval_opts = eval_options(Cur = TRUE,
#                                             Type = "ALL"),
#                    issp = TRUE)
# 
# names(dat1)
# plt1 <- dat1$tabs$plt
# spxy1 <- dat1$xyCur_PUBLIC
# table(plt1$INVYR)
# 
# # Display spatial output
# plot(sf::st_geometry(spxy1), pch = 16, cex = .5)

## -----------------------------------------------------------------------------
# # Add a filter to include only plots with Northern red oak forest type (FORTYPCD == 505)
# # Note: *allFilter* filters for plots and/or conditions for all states specified.
# 
# dat1b <- DBgetPlots(states = "Rhode Island",
#                     datsource = "datamart",
#                     eval = "FIA",
#                     eval_opts = eval_options(Cur = TRUE,
#                                              Type = "ALL"),
#                     issp = TRUE,
#                     allFilter = "FORTYPCD == 505")
# 
# names(dat1b)
# spxy1b <- dat1b$xyCur_PUBLIC
# dim(spxy1b)
# 
# # Display spatial output
# plot(sf::st_geometry(spxy1b), pch = 16, cex = .5, col="darkgreen")

## -----------------------------------------------------------------------------
# dat2 <- DBgetPlots(states = "Delaware",
#                    datsource = "datamart",
#                    eval = "FIA",
#                    eval_opts = eval_options(Cur = TRUE,
#                                             Type = "ALL"),
#                    issubp = TRUE,
#                    addplotgeom = TRUE)
# 
# names(dat2)
# tabs2 <- dat2$tabs
# plt2 <- tabs2$plt
# 
# ## subplot and subp_cond tables are added to tabs list
# names(tabs2)
# 
# ## PLOTGEOM data are appended to plt table (e.g., ALP_ADFORCD, FVS_VARIANT)
# head(plt2)

## -----------------------------------------------------------------------------
# 
# dat3 <- DBgetPlots(states = "Delaware",
#                    datsource = "datamart",
#                    eval = "FIA",
#                    eval_opts = eval_options(Cur = TRUE,
#                                             Type = "ALL"),
#                    savePOP = TRUE,
#                    othertables = c("POP_STRATUM", "POP_ESTN_UNIT"))
# 
# ## savePOP = TRUE, saves the POP_PLOT_STRATUM_ASSGN table used to select plots
# names(dat3)
# 
# ## pop_stratum and pop_estn_unit tables are added to tabs list
# tabs3 <- dat3$tabs
# names(tabs3)
# 

## -----------------------------------------------------------------------------
# DBgetPlots(states = "Rhode Island",
#            datsource = "datamart",
#            eval = "FIA",
#            eval_opts = eval_options(Cur = TRUE,
#                                     Type = "ALL"),
#            returndata = FALSE,
#            savedata = TRUE,
#            savedata_opts = savedata_options(outfolder = outfolder,
#                                             out_fmt = "csv",
#                                             overwrite_layer = TRUE))
# 
# ## Read in data from outfolder
# plt <- read.csv(file.path(outfolder, "plot.csv"), stringsAsFactors=FALSE)
# head(plt)

## -----------------------------------------------------------------------------
# dat5 <- DBgetPlots(states = "Rhode Island",
#                    datsource = "datamart",
#                    eval = "FIA",
#                    eval_opts = eval_options(Cur = TRUE,
#                                             Type = c("VOL", "CHNG", "P2VEG")))
# 
# names(dat5)
# tabs5 <- dat5$tabs
# names(tabs5)
# 
# ppsa5 <- dat5$pop_plot_stratum_assgn
# table(ppsa5$EVALID)
# 

## -----------------------------------------------------------------------------
# dat6 <- DBgetPlots(eval = "FIA",
#                    eval_opts = eval_options(Cur = TRUE,
#                                             evalid = c(101800, 101801, 101803)))
# 
# names(dat6)
# tabs6 <- dat6$tabs
# names(tabs6)
# 
# ppsa6 <- dat6$pop_plot_stratum_assgn
# table(ppsa6$EVALID)

## -----------------------------------------------------------------------------
# dat7 <- DBgetPlots(states = c("Connecticut"),
#                    eval = "FIA",
#                    eval_opts = eval_options(evalType = "ALL",
#                                             Endyr = 2017))
# 
# names(dat7)
# tabs7 <- dat7$tabs
# names(tabs7)
# 
# ppsa7 <- dat7$pop_plot_stratum_assgn
# table(ppsa7$EVALID)

## -----------------------------------------------------------------------------
# 
# dat8 <- DBgetPlots(states = "Vermont",
#                    eval = "custom",
#                    eval_opts = eval_options(invyrs = 2012:2014,
#                                             evalType = "ALL"))
# 
# names(dat8)
# tabs8 <- dat8$tabs
# names(tabs8)
# plt8 <- tabs8$plt
# 
# table(plt8$INVYR)

## -----------------------------------------------------------------------------
# dat9 <- DBgetPlots(states = "Wyoming",
#                     invtype = "PERIODIC",
#                     eval = "FIA",
#                     eval_opts = list(Cur = TRUE,
#                                      evalType = "VOL"))
# 
# names(dat9)
# tabs9 <- dat9$tabs
# names(tabs9)
# plt9 <- tabs9$plt
# 
# table(plt9$STATECD, plt9$INVYR)

## -----------------------------------------------------------------------------
# ## With only P2 plots (intensity1 = TRUE)
# dat10 <- DBgetPlots(states = "Vermont",
#                     eval = "FIA",
#                     eval_opts = list(Cur = TRUE,
#                                      Type = "ALL"),
#                     intensity1 = TRUE,
#                     issp = TRUE)
# 
# tabs10 <- dat10$tabs
# plt10 <- tabs10$plt
# 
# table(plt10$INVYR)
# spxy10 <- dat10$xyCur_PUBLIC
# 
# 
# # Display spatial output of public coordinates
# plot(sf::st_geometry(spxy10), pch = 16, cex = .5)

## -----------------------------------------------------------------------------
# strat1 <- DBgetStrata(states = "Wyoming",
#                       eval_opts = eval_options(Cur = TRUE))
# 
# names(strat1)
# 
# ## Look at plot assign data
# pltassgn1 <- strat1$pltassgn
# head(pltassgn1)
# 
# unique(pltassgn1$EVALID)
# strat1$evalid
# 
# ## Look at area data for estimation unit
# strat1$unitarea
# strat1$unitvar
# strat1$unitvar2
# strat1$areavar
# 
# ## Look at stratification data for estimation unit
# strat1$stratalut
# strat1$strvar
# strat1$getwtvar
# 

## -----------------------------------------------------------------------------
# strat2 <- DBgetStrata(eval_opts = eval_options(evalid = 561200))
# 
# unique(strat2$pltassgn$EVALID)
# strat2$evalid

## -----------------------------------------------------------------------------
# strat3 <- DBgetStrata(states = "Wyoming",
#                       eval_opts = eval_options(Endyr = 2014))
# 
# unique(strat3$pltassgn$EVALID)
# strat3$evalid

