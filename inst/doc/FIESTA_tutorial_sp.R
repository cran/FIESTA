## ----setup, include = F-------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(message = F, warning = F, eval = F)

## ---- include=FALSE-----------------------------------------------------------
#  # Sets up output folding
#  hooks = knitr::knit_hooks$get()
#  hook_foldable = function(type) {
#    force(type)
#    function(x, options) {
#      res = hooks[[type]](x, options)
#  
#      if (isFALSE(options[[paste0("fold.", type)]])) return(res)
#  
#      paste0(
#        "<details><summary>", type, "</summary>\n\n",
#        res,
#        "\n\n</details>"
#      )
#    }
#  }
#  knitr::knit_hooks$set(
#    output = hook_foldable("output"),
#    plot = hook_foldable("plot")
#  )

## ---- warning = F, message = F------------------------------------------------
#  library(FIESTA)

## -----------------------------------------------------------------------------
#  # File names for external spatial data
#  WYbhfn <- system.file("extdata",
#                        "sp_data/WYbighorn_adminbnd.shp",
#                        package = "FIESTA")
#  WYbhdistfn <- system.file("extdata",
#                            "sp_data/WYbighorn_districtbnd.shp",
#                            package = "FIESTA")
#  WYbhdist.att <- "DISTRICTNA"
#  
#  fornffn <- system.file("extdata",
#                         "sp_data/WYbighorn_forest_nonforest_250m.tif",
#                         package = "FIESTA")
#  demfn <- system.file("extdata",
#                       "sp_data/WYbighorn_dem_250m.img",
#                       package = "FIESTA")
#  
#  
#  # Other spatial layers used for examples, extracted using the raster package, getData function.
#  # County-level boundaries for USA and subset for Wyoming (Note: must have internet connection)
#  USAco <- raster::getData("GADM", country = "USA", level = 2)
#  WYco <- USAco[USAco$NAME_1 == "Wyoming",]
#  
#  head(USAco@data)

## ---- include = F-------------------------------------------------------------
#  outfolder <- tempdir()

## -----------------------------------------------------------------------------
#  ## Import external data shapefiles
#  WYbh <- spImportSpatial(WYbhfn)
#  WYbhdist <- spImportSpatial(WYbhdistfn)
#  
#  ## Display boundary
#  plot(sf::st_geometry(WYbhdist), border="blue")
#  plot(sf::st_geometry(WYbh), add=TRUE)
#  

## -----------------------------------------------------------------------------
#  ## Export Spatial Polygons layer to a shapefile
#  spExportSpatial(WYbh,
#                  savedata_opts = list(out_dsn = "WYbh.shp",
#                                       outfolder = outfolder,
#                                       overwrite_dsn = TRUE)
#                  )
#  

## -----------------------------------------------------------------------------
#  head(WYplt)
#  
#  WYspplt <- spMakeSpatialPoints(xyplt = WYplt,
#                                 xy.uniqueid = "CN",
#                                 xvar = "LON_PUBLIC",
#                                 yvar = "LAT_PUBLIC",
#                                 prj = "longlat",
#                                 datum = "NAD83")
#  head(WYspplt)
#  

## -----------------------------------------------------------------------------
#  WYspplt <- spMakeSpatialPoints(xyplt = WYplt,
#                                 xy.uniqueid = "CN",
#                                 xvar = "LON_PUBLIC",
#                                 yvar = "LAT_PUBLIC",
#                                 xy.crs = 4269
#                                 )
#  WYspplt
#  

## -----------------------------------------------------------------------------
#  ## Display output
#  plot(sf::st_geometry(WYbhdist))
#  plot(sf::st_geometry(WYspplt), add=TRUE)
#  
#  ## NOTE: To display multiple layers, all layers must be in the same coordinate system.
#  lapply(list(WYbh, WYbhdist, WYspplt), sf::st_crs)

## -----------------------------------------------------------------------------
#  WYspplt <- spMakeSpatialPoints(xyplt = WYplt,
#                                 xy.uniqueid = "CN",
#                                 xvar = "LON_PUBLIC",
#                                 yvar = "LAT_PUBLIC",
#                                 xy.crs = 4269,
#                                 exportsp = TRUE,
#                                 savedata_opts = list(
#                                      out_dsn = "spplt",
#                                      out_fmt = "shp",
#                                      outfolder = outfolder,
#                                      out_layer = "WYplots",
#                                      overwrite_layer = TRUE)
#                                 )
#  

## -----------------------------------------------------------------------------
#  sf::st_crs(WYspplt)
#  prj <- "+proj=utm +zone=12 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#  WYspplt.utm12 <- spReprojectVector(layer = WYspplt,
#                                     crs.new = prj)
#  
#  ## Check results
#  sf::st_crs(WYspplt.utm12)

## -----------------------------------------------------------------------------
#  ## Get points within Bighorn National Forest boundary (project on the fly)
#  
#  WYbhptslst <- spClipPoint(xyplt = WYplt,
#                            uniqueid = "CN",
#                            clippolyv = WYbh,
#                            spMakeSpatial_opts=list(xvar = "LON_PUBLIC",
#                                                    yvar = "LAT_PUBLIC",
#                                                    prj = "longlat",
#                                                    datum = "NAD83")
#                            )
#  
#  WYbhptslst <- spClipPoint(xyplt = WYplt,
#                            uniqueid = "CN",
#                            clippolyv = WYbh,
#                            savedata = TRUE,
#                            exportsp = TRUE,
#                            spMakeSpatial_opts=list(xvar = "LON_PUBLIC",
#                                                    yvar = "LAT_PUBLIC",
#                                                    prj = "longlat",
#                                                    datum = "NAD83"),
#                            savedata_opts = list(outfolder=outfolder,
#                                                 out_layer = "WYbh")
#                            )
#  
#  names(WYbhptslst)
#  WYbhspplt <- WYbhptslst$clip_xyplt
#  WYbhprj <- WYbhptslst$clip_polyv
#  

## -----------------------------------------------------------------------------
#  WYbhspplt
#  plot(sf::st_geometry(WYbhprj), border="red", lwd=2)
#  plot(sf::st_geometry(WYbhspplt), add=TRUE)

## -----------------------------------------------------------------------------
#  WYspplt <- spMakeSpatialPoints(xyplt = WYplt,
#                                 xy.uniqueid = "CN",
#                                 xvar = "LON_PUBLIC",
#                                 yvar = "LAT_PUBLIC",
#                                 xy.crs = 4269)
#  WYbhptslst <- spClipPoint(xyplt = WYspplt,
#                            uniqueid = "CN",
#                            clippolyv = WYbh)

## -----------------------------------------------------------------------------
#  WYbhptslst <- spClipPoint(xyplt = WYspplt,
#                            uniqueid = "CN",
#                            clippolyv = WYbh,
#                            othertabnms = c("WYcond", "WYtree"))
#  
#  names(WYbhptslst)
#  

## -----------------------------------------------------------------------------
#  ## Export clipped points
#  spExportSpatial(WYbhspplt,
#                  savedata_opts = list(out_layer = "WYbhpts",
#                                       outfolder = outfolder,
#                                       overwrite_layer = TRUE)
#                  )
#  

## -----------------------------------------------------------------------------
#  WYbhptslst <- spClipPoint(xyplt = WYplt,
#                            uniqueid = "CN",
#                            clippolyv = WYbh,
#                            othertabnms = c("WYcond", "WYtree"),
#                            exportsp = TRUE,
#                            spMakeSpatial_opts=list(xvar = "LON_PUBLIC",
#                                          yvar = "LAT_PUBLIC",
#                                          prj = "longlat",
#                                          datum = "NAD83"),
#                            savedata_opts = list(
#                                          outfolder = outfolder,
#                                          overwrite_layer = TRUE,
#                                          outfn.pre = "clip")
#                            )
#  names(WYbhptslst)
#  

## -----------------------------------------------------------------------------
#  WYbhco <- spClipPoly(polyv = USAco,
#                               clippolyv = WYbh)
#  

## -----------------------------------------------------------------------------
#  head(WYbhco)
#  plot(sf::st_geometry(WYbhco['NAME_2']), col = sf.colors(nrow(WYbhco)))
#  coords <- st_coordinates(st_centroid(st_geometry(WYbhco)))
#  text(coords[,"X"], coords[,"Y"], WYbhco[["NAME_2"]])

## -----------------------------------------------------------------------------
#  WYbhdist
#  WYbhMW <- WYbhdist[WYbhdist$DISTRICTNA == "Medicine Wheel Ranger District",]
#  
#  plot(st_geometry(WYbhdist))
#  plot(st_geometry(WYbhMW), border="red", add=TRUE)

## -----------------------------------------------------------------------------
#  WYbhMW.fornf <- spClipRast(fornffn,
#                             clippolyv = WYbhMW,
#                             outfolder = outfolder)

## -----------------------------------------------------------------------------
#  WYbhMWprj <- crsCompare(WYbhMW,
#                          rasterInfo(WYbhMW.fornf)$crs)$x
#  raster::plot(raster::raster(WYbhMW.fornf))
#  plot(st_geometry(WYbhMWprj),
#       border = "red",
#       add = TRUE)

## -----------------------------------------------------------------------------
#  extpolylst <- spExtractPoly(WYspplt,
#                              xy.uniqueid = "CN",
#                              polyvlst = WYbhdist)
#  WYspplt_bh <- extpolylst$spxyext
#  
#  dim(WYspplt)
#  dim(WYspplt_bh)
#  
#  head(WYspplt_bh)
#  plot(WYspplt_bh["DISTRICTNA"], pch = 8)
#  

## -----------------------------------------------------------------------------
#  extpolylst <- spExtractPoly(WYspplt,
#                              xy.uniqueid = "CN",
#                              polyvlst = WYbh,
#                              polyvarlst = c("FORESTNUMB", "FORESTNAME"),
#                              keepNA = FALSE)
#  WYspplt_bh2 <- extpolylst$spxyext
#  
#  dim(WYspplt)
#  dim(WYspplt_bh2)
#  
#  head(WYspplt_bh2)

## -----------------------------------------------------------------------------
#  extrastlst <- spExtractRast(WYspplt,
#                              rastlst = c(fornffn, demfn),
#                              xy.uniqueid = "CN",
#                              keepNA = FALSE)
#  WYspplt_dem <- extrastlst$sppltext
#  
#  dim(WYspplt)
#  dim(WYspplt_dem)
#  
#  head(WYspplt_dem)

## -----------------------------------------------------------------------------
#  WYspplt <- spMakeSpatialPoints(xyplt = WYplt,
#                                 xy.uniqueid = "CN",
#                                 xvar = "LON_PUBLIC",
#                                 yvar = "LAT_PUBLIC",
#                                 xy.crs = 4269)

## -----------------------------------------------------------------------------
#  library(raster)
#  dem <- raster(demfn)
#  slp <- raster::terrain(dem,
#                         opt = "slope",
#                         unit = "degrees",
#                         filename = paste0(outfolder, "/WYbh_slp.img"),
#                         overwrite = TRUE)
#  asp <- raster::terrain(dem,
#                         opt = "aspect",
#                         unit = "degrees",
#                         filename = paste0(outfolder, "/WYbh_asp.img"),
#                         overwrite = TRUE)

## -----------------------------------------------------------------------------
#  rastlst.cont <- c(demfn, slp, asp)
#  rastlst.cont.name <- c("dem", "slp", "asp")
#  rastlst.cat <- fornffn
#  rastlst.cat.name <- "fornf"
#  
#  modeldat <- spGetAuxiliary(xyplt = WYspplt,
#                             uniqueid = "CN",
#                             unit_layer = WYbhfn,
#                             unitvar = NULL,
#                             rastlst.cont = rastlst.cont,
#                             rastlst.cont.name = rastlst.cont.name,
#                             rastlst.cat = rastlst.cat,
#                             rastlst.cat.name = rastlst.cat.name,
#                             rastlst.cont.stat = "mean",
#                             asptransform = TRUE,
#                             rast.asp = asp,
#                             keepNA = FALSE,
#                             showext = FALSE,
#                             savedata = FALSE)
#  names(modeldat)
#  
#  pltassgn <- modeldat$pltassgn
#  unitzonal <- modeldat$unitzonal
#  unitvar <- modeldat$dunitvar
#  inputdf <- modeldat$inputdf
#  unitarea <- modeldat$unitarea
#  areavar <- modeldat$areavar
#  inputdf <- modeldat$inputdf
#  prednames <- modeldat$prednames
#  zonalnames <- modeldat$zonalnames
#  
#  unitvar
#  areavar
#  unitzonal
#  unitarea
#  head(pltassgn)
#  
#  prednames
#  zonalnames

## -----------------------------------------------------------------------------
#  WYbhxy <- spGetXY(bnd = WYbhfn,
#                    xy_datsource = "datamart",
#                    evalCur = TRUE,
#                    returnxy = TRUE)
#  names(WYbhxy)
#  
#  pltids <- WYbhxy$pltids
#  head(pltids)
#  
#  spxy <- WYbhxy$spxy
#  plot(st_geometry(spxy))
#  

## -----------------------------------------------------------------------------
#  
#  WYbhxyids <- spGetXY(bnd = WYbhfn,
#                       xy_datsource = "datamart",
#                       evalCur = TRUE,
#                       returnxy = FALSE)
#  names(WYbhxyids)
#  
#  pltids <- WYbhxyids$pltids
#  head(pltids)
#  
#  

## -----------------------------------------------------------------------------
#  WYbhdat <- spGetPlots(bnd = WYbhfn,
#                        states = "Wyoming",
#                        datsource = "datamart",
#                        evalCur = TRUE,
#                        istree = FALSE)
#  names(WYbhdat)

## -----------------------------------------------------------------------------
#  WYspplt <- spMakeSpatialPoints(xyplt = WYplt,
#                                 xy.uniqueid = "CN",
#                                 xvar = "LON_PUBLIC",
#                                 yvar = "LAT_PUBLIC",
#                                 xy.crs = 4269)
#  xyplt <- WYspplt

## -----------------------------------------------------------------------------
#  unitdat.bh <- spGetEstUnit(xyplt = WYplt,
#                             uniqueid = "CN",
#                             unit_layer = WYbhfn,
#                             spMakeSpatial_opts=list(xvar = "LON_PUBLIC",
#                                                     yvar = "LAT_PUBLIC",
#                                                     prj = "longlat",
#                                                     datum = "NAD83")
#                             )
#  
#  names(unitdat.bh)
#  unitarea.bh <- unitdat.bh$unitarea
#  unitvar.bh <- unitdat.bh$unitvar
#  areavar.bh <- unitdat.bh$areavar
#  
#  unitarea.bh
#  unitvar.bh
#  areavar.bh
#  

## -----------------------------------------------------------------------------
#  unitdat.bhdist <- spGetEstUnit(xyplt = WYplt,
#                                 uniqueid = "CN",
#                                 unit_layer = WYbhdistfn,
#                                 unitvar = "DISTRICTNA",
#                                 spMakeSpatial_opts=list(xvar = "LON_PUBLIC",
#                                                         yvar = "LAT_PUBLIC",
#                                                         prj = "longlat",
#                                                         datum = "NAD83")
#                                 )
#  
#  names(unitdat.bhdist)
#  unitarea.bhdist <- unitdat.bhdist$unitarea
#  unitvar.bhdist <- unitdat.bhdist$unitvar
#  areavar.bhdist <- unitdat.bhdist$areavar
#  
#  unitarea.bhdist
#  unitvar.bhdist
#  areavar.bhdist
#  

## -----------------------------------------------------------------------------
#  WYspplt <- spMakeSpatialPoints(xyplt = WYplt,
#                                 xy.uniqueid = "CN",
#                                 xvar = "LON_PUBLIC",
#                                 yvar = "LAT_PUBLIC",
#                                 xy.crs = 4269)

## -----------------------------------------------------------------------------
#  stratlst <- spGetStrata(WYspplt,
#                          uniqueid = "CN",
#                          unit_layer = WYbhfn,
#                          strattype = "RASTER",
#                          strat_layer = fornffn)
#  names(stratlst)
#  stratlst$stratalut

## -----------------------------------------------------------------------------
#  WYbh <- spImportSpatial(WYbhfn)
#  
#  polyUnion <- spUnionPoly(polyv1 = USAco[USAco$NAME_1 == "Wyoming",],
#                           polyv2 = WYbh,
#                           areacalc=TRUE)
#  
#  plot(st_geometry(polyUnion))
#  head(polyUnion)

## -----------------------------------------------------------------------------
#  WYbhdist <- spImportSpatial(WYbhdistfn)
#  
#  zonallst <- spZonalRast(polyv = WYbhdist,
#                          polyv.att = "DISTRICTNA",
#                          rastfn = demfn,
#                          zonalstat = c("mean", "sum"))
#  names(zonallst)
#  
#  zonalext <- zonallst$zonalext
#  outname <- zonallst$outname
#  rasterfile <- zonallst$rasterfile
#  
#  head(zonalext)
#  outname
#  rasterfile

## ---- include = FALSE---------------------------------------------------------
#  # deletes raster data
#  unlink("gadm36_USA_2_sp.rds")

