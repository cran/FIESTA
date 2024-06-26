#' Small area module - Compile population data for SA module.
#' 
#' Compile population data for input to the modSA* modules.
#' 
#' If variables are NULL, then it will prompt user to input variables.
#' 
#' Necessary variables:\cr \tabular{llll}{ \tab \bold{Data} \tab
#' \bold{Variable} \tab \bold{Description}\cr 
#' \tab tree \tab tuniqueid \tab Unique identifier for each plot, to link to 
#' pltassgn (e.g. PLT_CN).\cr 
#' \tab \tab CONDID \tab Unique identifier of each condition on plot, to link to
#' cond.  Set CONDID=1, if only 1 condition per plot.\cr 
#' \tab \tab TPA_UNADJ \tab Number of trees per acre each sample tree 
#' represents (e.g. DESIGNCD=1: TPA_UNADJ=6.018046 for trees on subplot; 
#' 74.965282 for trees on microplot).\cr 
#' \tab cond \tab cuniqueid \tab Unique identifier for each plot, to link to 
#' pltassgn (e.g. PLT_CN).\cr 
#' \tab \tab CONDID \tab Unique identifier of each condition on plot. Set 
#' CONDID=1, if only 1 condition per plot.\cr 
#' \tab \tab CONDPROP_UNADJ \tab Unadjusted proportion of condition on
#' each plot.  Set CONDPROP_UNADJ=1, if only 1 condition per plot.\cr 
#' \tab \tab COND_STATUS_CD \tab Status of each forested condition on plot (i.e.
#' accessible forest, nonforest, water, etc.)\cr 
#' \tab \tab NF_COND_STATUS_CD \tab If ACI=TRUE. Status of each nonforest 
#' condition on plot (i.e. accessible nonforest, nonsampled nonforest)\cr 
#' \tab \tab SITECLCD \tab If landarea=TIMBERLAND. Measure of site productivity.\cr 
#' \tab \tab RESERVCD \tab If landarea=TIMBERLAND. Reserved status.\cr 
#' \tab \tab SUBPROP_UNADJ \tab Unadjusted proportion of subplot conditions 
#' on each plot. Set SUBPROP_UNADJ=1, if only 1 condition per subplot.\cr 
#' \tab \tab MICRPROP_UNADJ \tab If microplot tree attributes. Unadjusted 
#' proportion of microplot conditions on each plot. Set MICRPROP_UNADJ=1, 
#' if only 1 condition per microplot.\cr 
#' \tab \tab MACRPROP_UNADJ \tab If macroplot tree attributes. Unadjusted 
#' proportion of macroplot conditions on each plot. Set MACRPROP_UNADJ=1, 
#' if only 1 condition per macroplot.\cr 
#' \tab pltassgn \tab puniqueid \tab Unique identifier for each plot, to 
#' link to cond (e.g. CN).\cr 
#' \tab \tab STATECD \tab Identifies state each plot is located in.\cr 
#' \tab \tab INVYR \tab Identifies inventory year of each plot.\cr 
#' \tab \tab PLOT_STATUS_CD \tab Status of each plot (i.e. sampled, nonsampled).  
#' If not included, all plots are assumed as sampled.\cr }
#' 
#' For available reference tables: sort(unique(FIESTAutils::ref_codes$VARIABLE)) \cr
#' 
#' @param popType String. Type of evaluation(s) to include in population data.
#' Note: currently only c('CURR', 'VOL', 'LULC') are available. See details
#' below for descriptions of each.
#' @param popTabs List of population tables the user would like returned.
#'  See help(popTables) for a list of options.
#' @param popTabIDs List of unique IDs corresponding to the population tables
#' that the user has requested. See help(popTableIDs) for a list of
#' options.
#' @param popFilter List of population filters. See help(popFilters) for a 
#' list of options. 
#' @param SAdoms sf object. SA domains with attributes for joining.
#' @param smallbnd sf object. small bound.
#' @param smallbnd.domain String. Name of attribute defining domain attribute.
#' @param largebnd.unique String. Name of the large boundary unique identifer
#' to define plots within a model extent. If NULL, all plots are used for model
#' extent.
#' @param pltassgn DF/DT, comma-separated values (CSV) file(*.csv), or layer in
#' dsn, Can also be a shapefile(*.shp) with one record per plot, a spatial
#' layer in dsn, or a sf R object. Plot-level assignment of estimation unit
#' and/or strata. Optional.
#' @param dsn String. Name of database where tree, cond, and pltassgn tables
#' reside.  The dsn varies by driver. See gdal OGR vector formats
#' (https://www.gdal.org/ogr_formats.html).
#' @param pltassgnid String. Unique identifier of plot in pltassgn.
#' @param pjoinid String. Join variable in plot to match pltassgnid. Does not
#' need to be uniqueid. If using most current XY coordinates for plot
#' assignments, use identifier for plot (e.g., PLOT_ID).
#' @param areawt String. Name of variable for summarizing area weights (e.g.,
#' CONDPROP_UNADJ).
#' @param adjplot Logical. If TRUE, adjusts for nonresponse at plot-level.
#' @param dunitvar String. Name of the domain unit variable in cond, plt, or
#' pltassgn with domain unit assignment for each plot.
#' @param dunitarea Numeric or DF. Total area by domain unit.
#' @param areavar String. Name of area variable in unitarea. Default="ACRES".
#' @param dunitzonal DF/DT. Data frame with zonal auxiliary information by
#' domain unit. For continuous data, means by domain unit; for categorical
#' data, proportion of class by domain unit.
#' @param prednames String vector. Name(s) of predictor variables to use in
#' model. 
#' @param predfac String vector. Name(s) of factor predictor variables to use
#' in model. Names will change in output depending on number of categories.
#' @param addxy Logical. If TRUE, adds X/Y attributes to pltassgn.
#' @param returndata Logical. If TRUE, returns data objects.
#' @param savedata Logical. If TRUE, saves table(s) to outfolder. 
#' @param saveobj Logical. If TRUE, saves returned list object to outfolder.
#' @param objnm String. Name of *.rds object.
#' @param unit_opts List. See help(unit_options()) for a list of options.
#' @param savedata_opts List. See help(savedata_options()) for a list
#' of options.
#' @param SAdata R List object. Output data list components from
#' FIESTA::SAdata().
#' @param pltdat R List object. Output data list components from
#' FIESTA::spGetPlots().
#' @param auxdat R List object. Output data list components from
#' FIESTA::spGetAuxiliary().
#' @param gui Logical. If gui, user is prompted for parameters.
#' @param ... For extendibility. 
#' @return A list with population data for Small-Area estimates.
#' 
#' \item{SAdomsdf}{ Data frame. Attribute table from SAdoms spatial layer.
#' Includes DOMAIN and AOI attributes. DOMAIN represents modeling domains.  AOI
#' identifies the small area of interest. } 
#' \item{condx}{ Data frame. Condition-level data with condition proportions, 
#' domain and predictor assignments, and adjusted condition proportions, 
#' if adjplot = TRUE. } 
#' \item{pltcondx}{ Data frame. Plot/Condition data used for estimation. }
#' \item{cuniqueid}{ String. Unique identifier of plot in condx and pltcondx. }
#' \item{condid}{ String. Unique identifier of condition in condx and pltcondx. } 
#' \item{treex}{ Data frame. If esttype='TREE', tree-level data, including
#' adjustment factors, if adjplot = TRUE. } 
#' \item{tuniqueid}{ String. If esttype='TREE', unique identifier of plot in 
#' treex. } 
#' \item{ACI.filter}{ String. If ACI=FALSE, ACI.filter="COND_STATUS_CD == 1" . } 
#' \item{dunitarea}{ Data frame. Area by model domain unit. } 
#' \item{areavar}{ String. Name of area variable in dunitarea. } 
#' \item{dunitvar}{ String. Name of variable defining model domain units in 
#' dunitarea. } 
#' \item{dunitlut}{ Data frame. Table of model domain units with zonal statistics 
#' of predictor values, number of plots by domain unit. } 
#' \item{prednames}{ String vector. Name of variables in dunitlut and condx 
#' defining potential predictors for small area estimation. } 
#' \item{plotsampcnt}{ Data frame. Number of plots by PLOT_STATUS_CD. } 
#' \item{condsampcnt}{ Data frame. Number of conditions by COND_STATUS_CD. } 
#' \item{states}{ String. State names in dataset. }
#' \item{invyrs}{ String. Range of inventory years in dataset. }
#' \item{adjtree}{ Logical. If TRUE, treex includes adjustment factors. }
#' @note
#' 
#' ADJUSTMENT FACTOR:\cr The adjustment factor is necessary to account for
#' nonsampled conditions.  For model-based estimation, we calculate adjustment
#' factors by plot.
#' 
#' It is calculated by dividing 1 / summed condition proportions by plot. An
#' adjustment factor is determined for each tree based on the size of the plot
#' it was measured on. This is identified using TPA_UNADJ as follows:
#' 
#' \tabular{llr}{ \tab \bold{PLOT SIZE} \tab \bold{TPA_UNADJ} \cr 
#' \tab SUBPLOT \tab 6.018046 \cr 
#' \tab MICROPLOT \tab 74.965282 \cr 
#' \tab MACROPLOT \tab 0.999188 \cr }
#' 
#' If ACI=FALSE, only nonsampled forest conditions are accounted for in the
#' adjustment factor. \cr 
#' If ACI=TRUE, the nonsampled nonforest conditions are removed as well and 
#' accounted for in adjustment factor.  This is if you are interested in 
#' estimates for all lands or nonforest lands in the All-Condition-Inventory.
#' @author Tracey S. Frescino, Paul L. Patterson
#' @keywords data
#' @examples 
#' \donttest{
#' # NOTE: FIA data objects used in these examples are stored in `FIESTA`, but
#' # can be generated for populations of interest by the user with functions in
#' # `FIESTA` such as `spGetPlots()`, `spGetAuxiliary()`, etc. For more
#' # information, see `FIESTA`'s extensive vignettes.
#'  
#' # Population data for counties in Wyoming
#' modSApop(popTabs = list(tree = FIESTA::WYtree,
#'                         cond = FIESTA::WYcond),
#'          pltassgn = FIESTA::WYpltassgn,
#'          pltassgnid = "CN",
#'          dunitarea = FIESTA::WYunitarea,
#'          dunitvar = "ESTN_UNIT",
#'          dunitzonal = FIESTA::WYunitzonal,
#'          prednames = c("dem", "tcc", "tpi", "tnt"),
#'          predfac = "tnt")
#'          
#' # Adding seedling data as well
#' modSApop(popTabs = list(tree = FIESTA::WYtree,
#'                         cond = FIESTA::WYcond,
#'                         seed = FIESTA::WYseed),
#'          pltassgn = FIESTA::WYpltassgn,
#'          pltassgnid = "CN",
#'          dunitarea = FIESTA::WYunitarea,
#'          dunitvar = "ESTN_UNIT",
#'          dunitzonal = FIESTA::WYunitzonal,
#'          prednames = c("dem", "tcc", "tpi", "tnt"),
#'          predfac = "tnt")
#' }
#' @export modSApop
modSApop <- function(popType = "VOL",
                     popTabs = popTables(),
                     popTabIDs = popTableIDs(), 
                     popFilter = popFilters(),
                     pltassgn = NULL,
                     pltassgnid = "PLT_CN", 
                     dsn = NULL, 
                     pjoinid = "CN", 
                     areawt = "CONDPROP_UNADJ", 
                     adjplot = TRUE,
                     dunitvar = NULL, 
                     dunitarea = NULL, 
                     areavar = "ACRES",
                     dunitzonal = NULL, 
                     prednames = NULL, 
                     predfac = NULL, 
                     addxy = FALSE,
                     returndata = TRUE,
                     savedata = FALSE, 
                     saveobj = FALSE, 
                     objnm = "SApopdat", 
                     unit_opts = NULL, 
                     savedata_opts = NULL, 
                     SAdoms = NULL, 
                     smallbnd = NULL, 
                     smallbnd.domain = NULL, 
                     largebnd.unique = NULL,
                     SAdata = NULL, 
                     pltdat = NULL, 
                     auxdat = NULL, 
                     gui = FALSE, 
                     ...) {

  ##################################################################################
  ## DESCRIPTION:
  ## Generates population data for small area estimation
  ## - check population data
  ## - check dunitarea data
  ## - calculate plot-level adjustment factors by dividing 1 by summed proportions in plot
  ## VALUE:
  ## - return all data needed for input to modSAest() function
  ##################################################################################

  ## CHECK GUI - IF NO ARGUMENTS SPECIFIED, ASSUME GUI=TRUE
  if (nargs() == 0) gui <- TRUE

  ## If gui.. set variables to NULL
  if (gui) {
    areavar=strata=strvar=getwt=cuniqueid=ACI=tuniqueid=savedata=unitvar <- NULL
  }
  
  ## Set parameters
  adjtree <- FALSE
  returnSApopdat <- FALSE
  nonsamp.pfilter=nonsamp.cfilter <- NULL 
  returnlst <- list(module = "SA")
  adj <- ifelse(adjplot, "plot", "none")
  areawt2 <- NULL
  pvars2keep <- NULL

  # dunitvar2=NULL
  # pvars2keep=NULL
  # cvars2keep=NULL
  # adj="plot"
  # ACI=FALSE
  # gui <- FALSE 
 
  ## Set global variables
  ONEUNIT=n.total=n.strata=strwt=TOTAL=stratcombinelut <- NULL
 
  ##################################################################
  ## CHECK PARAMETER NAMES
  ##################################################################

  ## Check input parameters
  input.params <- names(as.list(match.call()))[-1]
  if (!all(input.params %in% names(formals(modSApop)))) {
    miss <- input.params[!input.params %in% formals(modSApop)]
    stop("invalid parameter: ", toString(miss))
  }

  
  ## Check parameter lists
  pcheck.params(input.params, unit_opts=unit_opts, savedata_opts=savedata_opts)
  
  ## Set unit defaults
  unit_defaults_list <- formals(unit_options)[-length(formals(unit_options))]
  
  for (i in 1:length(unit_defaults_list)) {
    assign(names(unit_defaults_list)[[i]], unit_defaults_list[[i]])
  }
  
  ## Set user-supplied unit values
  if (length(unit_opts) > 0) {
    for (i in 1:length(unit_opts)) {
      if (names(unit_opts)[[i]] %in% names(unit_defaults_list)) {
        assign(names(unit_opts)[[i]], unit_opts[[i]])
      } else {
        stop(paste("Invalid parameter: ", names(unit_opts)[[i]]))
      }
    }
  }
  
  ## Set popFilters defaults
  popFilters_defaults_list <- formals(popFilters)[-length(formals(popFilters))]
  
  for (i in 1:length(popFilters_defaults_list)) {
    assign(names(popFilters_defaults_list)[[i]], popFilters_defaults_list[[i]])
  }
  
  ## Set user-supplied popFilters values
  popFilter2 <- popFilters_defaults_list
  if (length(popFilter) > 0) {
    for (i in 1:length(popFilter)) {
      if (names(popFilter)[[i]] %in% names(popFilters_defaults_list)) {
		    popFilter2[[names(popFilter)[[i]]]] <- popFilter[[i]]
      } else {
        stop(paste("Invalid parameter: ", names(popFilter)[[i]]))
      }
    }
  }
  
  ## Set savedata defaults
  savedata_defaults_list <- formals(savedata_options)[-length(formals(savedata_options))]
  
  for (i in 1:length(savedata_defaults_list)) {
    assign(names(savedata_defaults_list)[[i]], savedata_defaults_list[[i]])
  }
  
  ## Set user-supplied savedata values
  if (length(savedata_opts) > 0) {
    if (!savedata && !saveobj) {
      message("savedata=FALSE with savedata parameters... no data are saved")
    }
    for (i in 1:length(savedata_opts)) {
      if (names(savedata_opts)[[i]] %in% names(savedata_defaults_list)) {
        assign(names(savedata_opts)[[i]], savedata_opts[[i]])
      } else {
        stop(paste("Invalid parameter: ", names(savedata_opts)[[i]]))
      }
    }
  }

  
  ##################################################################
  ## CHECK PARAMETER INPUTS
  ##################################################################
  
  ## Check addxy 
  addxy <- pcheck.logical(addxy, varnm="addxy", 
    title="Add XY?", first="NO", gui=gui, stopifnull=TRUE)

    ## Check savedata 
  savedata <- pcheck.logical(savedata, varnm="savedata", 
		title="Save data tables?", first="YES", gui=gui, stopifnull=TRUE)

  ## Check saveobj 
  saveobj <- pcheck.logical(saveobj, varnm="saveobj", 
		title="Save SApopdat object?", first="YES", gui=gui, stopifnull=TRUE)
  

  ## Check output
  ########################################################
  if (savedata) {
    outlst <- pcheck.output(outfolder=outfolder, out_dsn=out_dsn, 
                  out_fmt=out_fmt, outfn.pre=outfn.pre, outfn.date=outfn.date, 
                  overwrite_dsn=overwrite_dsn, overwrite_layer=overwrite_layer,
                  add_layer=add_layer, append_layer=append_layer, gui=gui)
    outfolder <- outlst$outfolder
    out_dsn <- outlst$out_dsn
    out_fmt <- outlst$out_fmt
    overwrite_layer <- outlst$overwrite_layer
    append_layer <- outlst$append_layer
    outfn.date <- outlst$outfn.date
    outfn.pre <- outlst$outfn.pre
  } 

  if (saveobj) {
    outobj_fmtlst <- c('rds', 'rda')
    outobj_fmt <- pcheck.varchar(var2check=outobj_fmt, varnm="outobj_fmt", gui=gui,
		checklst=outobj_fmtlst, caption="outobj_fmt", multiple=FALSE, stopifnull=TRUE)

    if (is.null(objnm)) {
      objnm <- "SApopdat"
    }
    #if (append_layer) overwrite_layer <- FALSE
    if (append_layer) message("currently cannot append to object lists")
    objfn <- getoutfn(outfn=objnm, ext=outobj_fmt, outfolder=outfolder, 
		overwrite=overwrite_layer, outfn.pre=outfn.pre, outfn.date=outfn.date)
  }


  ## Check popType
  ########################################################
  #evalTyplst <- c("ALL", "CURR", "VOL", "LULC", "P2VEG", "INV", "GRM", "DWM")
  DWM_types <- c("CWD", "FWD_SM", "FWD_LG", "DUFF")
  evalTyplst <- c("ALL", "CURR", "VOL", "LULC", "P2VEG", "INV", "DWM", "CHNG", "GRM")
  popType <- pcheck.varchar(var2check=popType, varnm="popType", gui=gui,
		checklst=evalTyplst, caption="popType", multiple=FALSE, stopifinvalid=FALSE)
  if (is.null(popType)) {
    message("popType is invalid... must be from following list:\n", toString(evalTyplst))
  }
  popevalid <- as.character(popFilter2$evalid)
  if (!is.null(popevalid)) {
    substr(popevalid, nchar(popevalid)-1, nchar(popevalid)) <- 
		formatC(FIESTAutils::ref_popType[FIESTAutils::ref_popType$popType %in% popType, "EVAL_TYP_CD"], 
		width=2, flag="0")
    #evalid <- as.character(evalid)
    #substr(evalid, nchar(evalid)-1, nchar(evalid)) <- "01"
  } 
  if (popType %in% c("GROW", "MORT", "REMV")) {
    popType <- "GRM"
  }

  ###################################################################################
  ## Load data
  ###################################################################################
  if (!is.null(SAdata)) {
    list.items <- c("tabs", "pltassgn", "pltassgnid", "pjoinid", "unitarea",
		                "unitvar", "areavar", "unitzonal")
    SAdata <- pcheck.object(SAdata, "SAdata", list.items=list.items)
    SAdoms <- SAdata$bnd
    #smallbnd <- SAdata$smallbnd
    popTabs <- SAdata$tabs
    popTabIDs <- SAdata$tabIDs
    pltassgn <- SAdata$pltassgn
    pltassgnid <- SAdata$pltassgnid
    dunitarea <- SAdata$unitarea
    areavar <- SAdata$areavar
    dunitvar <- SAdata$unitvar
    areavar <- SAdata$areavar
    dunitzonal <- SAdata$unitzonal
    puniqueid <- SAdata$puniqueid
    pjoinid <- SAdata$pjoinid
    predfac <- SAdata$predfac
    spxy <- SAdata$spxy
    xy.uniqueid <- SAdata$xy.uniqueid
    pvars2keep <- SAdata$vars2keep

    if (is.null(prednames)) {
      prednames <- SAdata$prednames
    } else {
      if (!all(prednames %in% SAdata$prednames))
        stop("invalid prednames: ", 
		            toString(prednames[!prednames %in% SAdata$prednames]))
      predfac <- predfac[predfac %in% prednames]
    }
  } else {
    if (!is.null(pltdat)) {
      popTabs <- pltdat$tabs
      popTabIDs <- pltdat$tabIDs
      pjoinid <- pltdat$pjoinid
      spxy <- pltdat$spxy
      xy.uniqueid <- pltdat$xy.uniqueid
    }
    if (!is.null(auxdat)) {
      list.items <- c("pltassgn", "unitzonal", "unitvar", "prednames", "unitarea")
      auxdat <- pcheck.object(auxdat, "auxdat", list.items=list.items)
      pltassgn <- auxdat$pltassgn
      pltassgnid <- auxdat$pltassgnid
      dunitarea <- auxdat$unitarea
      dunitvar <- auxdat$unitvar
      areavar <- auxdat$areavar
      dunitzonal <- auxdat$unitzonal
      zonalnames <- auxdat$zonalnames
      predfac <- auxdat$predfac
      areavar <- auxdat$areavar

      if (is.null(prednames)) {
        prednames <- auxdat$prednames
      } else {
        if (!all(prednames %in% auxdat$prednames))
          stop("invalid prednames: ", toString(prednames[!prednames %in% auxdat$prednames]))
        predfac <- predfac[predfac %in% prednames]
      }
    } 
  }
  ## check addxy
  if (addxy && is.null(spxy)) {
    message("no xy available...  use returnxy=TRUE when getting plot data")
  }
  
  ## Set user-supplied popTable values 
  popTables_defaults_list <- formals(popTables)[-length(formals(popTables))]
  if (length(popTabs) > 0) {
    for (i in 1:length(popTabs)) {
      if (names(popTabs)[[i]] %in% names(popTables_defaults_list)) {
        assign(names(popTabs)[[i]], popTabs[[i]])
      } else {
        stop(paste("Invalid parameter: ", names(popTabs)[[i]]))
      }
    }
  } else {
    stop("need to include popTabs")
  }
  
  list.items <- {}
  if (popType == "LULC") {
    list.items <- c(list.items, "lulcx")
  }
  if (popType == "P2VEG") {
    list.items <- c(list.items, "vsubpspp", "vsubpstr", "subplot", "subp_cond")
  }
  popTabs <- pcheck.object(popTabs, "popTabs", list.items=list.items)
  
  
  ## Set user-supplied popTabIDs values
  ### Check for invalid parameters first
  popTableIDs_defaults_list <- formals(popTableIDs)[-length(formals(popTableIDs))]
  #for (i in 1:length(popTabIDs)) {
  #  if (!(names(popTabIDs)[[i]] %in% names(popTableIDs_defaults_list))) {
  #    stop(paste("Invalid parameter: ", names(popTabIDs)[[i]]))
  #  }
  #}
  for (i in 1:length(popTableIDs_defaults_list)) {
    if (names(popTableIDs_defaults_list)[[i]] == "cond") {
      assign("cuniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "plt") {
      assign("puniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "tree") {
      assign("tuniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "seed") {
      assign("suniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "vsubpspp") {
      assign("vsppuniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "vsubpstr") {
      assign("vstruniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "invsubp") {
      assign("invuniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "subplot") {
      assign("subpuniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "subp_cond") {
      assign("subcuniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "cond_dwm_calc") {
      assign("dwmuniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "grm") {
      assign("grmuniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "plot_pplot") {
      assign("pplotuniqueid", popTableIDs_defaults_list[[i]])
    }
    if (names(popTableIDs_defaults_list)[[i]] == "cond_pcond") {
      assign("pconduniqueid", popTableIDs_defaults_list[[i]])
    }
  }

  ## Check SAdoms
  if (!is.null(SAdoms) && !"sf" %in% class(SAdoms)) {
    stop("invalid SAdoms")
  }
  
  
  pvars2keep <- unique(c(largebnd.unique, "AOI", pvars2keep))
  ###################################################################################
  ## CHECK PLOT PARAMETERS AND DATA
  ## Generate table of sampled/nonsampled plots and conditions
  ## Remove nonsampled plots (if nonsamp.pfilter != "NONE")
  ## Applies plot filters
  ###################################################################################
  pltcheck <- check.popdataPLT(dsn=dsn, tabs=popTabs, tabIDs=popTabIDs, 
      pltassgn=pltassgn, pltassgnid=pltassgnid, pjoinid=pjoinid, 
      module="SA", popType=popType, popevalid=popevalid, adj=adj, 
	    popFilter=popFilter2, nonsamp.pfilter=nonsamp.pfilter, 
	    unitarea=dunitarea, areavar=areavar, unitvar=dunitvar, 
	    unitvar2=unitvar2, areaunits=areaunits, unit.action=unit.action, 
      prednames=prednames, predfac=predfac, pvars2keep=pvars2keep)
  if (is.null(pltcheck)) return(NULL)
  pltassgnx <- pltcheck$pltassgnx
  pltassgnid <- pltcheck$pltassgnid
  pfromqry <- pltcheck$pfromqry
  palias <- pltcheck$palias
  pjoinid <- pltcheck$pjoinid
  whereqry <- pltcheck$whereqry
  ACI <- pltcheck$ACI
  pltx <- pltcheck$pltx
  puniqueid <- pltcheck$puniqueid
  dunitvar <- pltcheck$unitvar
  dunitvar2 <- pltcheck$unitvar2
  dunitarea <- pltcheck$unitarea
  areavar <- pltcheck$areavar
  areaunits <- pltcheck$areaunits
  dunit.action <- pltcheck$unit.action
  prednames <- pltcheck$prednames
  predfac <- pltcheck$predfac
  plotsampcnt <- pltcheck$plotsampcnt
  states <- pltcheck$states
  invyrs <- pltcheck$invyrs
  dbconn <- pltcheck$dbconn
  
  # subset pvars2keep 
  if (!"AOI" %in% names(pltassgnx)) {
    pltassgnx$AOI <- 1
  }

  if (!is.null(pvars2keep)) {
    pvars2keep <- pvars2keep[pvars2keep %in% names(pltx) & !pvars2keep %in% names(pltassgnx)]
    if (length(pvars2keep) > 0) {
      pltassgnx <- merge(pltassgnx, pltx[, c(puniqueid, pvars2keep), with=FALSE], 
                       by.x=pltassgnid, by.y=puniqueid)
      pltx <- pltx[, names(pltx)[!names(pltx) %in% pvars2keep], with=FALSE]
      setcolorder(pltassgnx, c(pltassgnid, pvars2keep, 
                          names(pltassgnx)[!names(pltassgnx) %in% c(pltassgnid, pvars2keep)]))
    }
  }

  if (ACI) {
    nfplotsampcnt <- pltcheck$nfplotsampcnt
  }

  if (popType %in% c("ALL", "CURR", "AREA", "VOL")) {
    ###################################################################################
    ## Check parameters and data for popType AREA/VOL
    ###################################################################################
    popcheck <- check.popdataVOL(gui=gui, 
        tabs=popTabs, tabIDs=popTabIDs, pltassgnx=pltassgnx, 
        pfromqry=pfromqry, palias=palias, pjoinid=pjoinid, 
			  whereqry=whereqry, adj=adj, ACI=ACI, 
			  pltx=pltx, puniqueid=puniqueid, dsn=dsn, dbconn=dbconn,
        condid="CONDID", nonsamp.cfilter=nonsamp.cfilter, 
			  areawt=areawt, areawt2=areawt2, cvars2keep="AOI")
    if (is.null(popcheck)) return(NULL)
    condx <- popcheck$condx
    pltcondx <- popcheck$pltcondx
    treef <- popcheck$treef
    seedf <- popcheck$seedf
    cuniqueid <- popcheck$cuniqueid
    condid <- popcheck$condid
    tuniqueid <- popcheck$tuniqueid
    ACI.filter <- popcheck$ACI.filter
    condsampcnt <- popcheck$condsampcnt
    areawt <- popcheck$areawt
    tpropvars <- popcheck$tpropvars
  }
  
  if (is.null(treef) && is.null(seedf)) {
    stop("must include tree data")
  }

  ###################################################################################
  ## Check auxiliary data
  ###################################################################################
  auxdat <- check.auxiliary(module = "SA",
                    pltx = pltassgnx, puniqueid = pltassgnid, 
                    unitvar = dunitvar, 
                    unitarea = dunitarea, 
                    areavar = areavar, 
                    minplotnum.unit = minplotnum.unit, 
                    unit.action=dunit.action,
                    auxlut = dunitzonal, 
                    prednames = prednames, 
                    predfac = predfac, 
                    makedummy = TRUE,  
                    standardize = TRUE,                 
                    auxtext = "dunitlut", removetext = "dunitarea")  
  pltassgnx <- setDT(auxdat$pltx)
  dunitarea <- auxdat$unitarea
  dunitvar <- auxdat$unitvar
  dunitlut <- auxdat$auxlut
  prednames <- auxdat$prednames
  predfac <- auxdat$predfac
  if (is.null(key(pltassgnx))) setkeyv(pltassgnx, pltassgnid)
  
  # subset pvars2keep 
  if (!"AOI" %in% names(dunitlut)) {
    dunitlut$AOI <- 1
  }
  
  ## Change names based on data.frame names
  dunitlutcols <- which(names(dunitlut) %in% prednames)
  pltassgnxcols <- which(names(pltassgnx) %in% prednames)
  dfnames <- colnames(data.frame(dunitlut[, dunitlutcols, with=FALSE]))
  names(dunitlut)[dunitlutcols] <- dfnames
  names(pltassgnx)[pltassgnxcols] <- dfnames
  prednames <- dfnames
  
  ## Change name of dunitvar to DOMAIN if not already
  if (dunitvar != "DOMAIN") {
    setnames(dunitlut, dunitvar, "DOMAIN")
    setnames(dunitarea, dunitvar, "DOMAIN")
    setnames(pltassgnx, dunitvar, "DOMAIN")
    dunitvar <- "DOMAIN"
  }

  ###################################################################################
  ## CALCULATE ADJUSTMENT FACTORS FOR NONSAMPLED CONDITIONS
  ## If adj="samp", calculate adjustment factors by strata and estimation unit
  ## If adj="plot", calculate adjustment factors by plot
  ## adjfac = 1 / summed condition proportions (by plot-size), for area and trees
  ###################################################################################
  ## Returns:
  ##  Adjusted condition proportion (CONDPROP_ADJ) appended to condx
  ###################################################################################
  ## Merge plot strata info to condx
  if (is.null(key(condx))) setkeyv(condx, c(cuniqueid, condid))
  condx <- condx[pltassgnx[, c(pltassgnid, dunitvar, prednames), with=FALSE]]
  setkeyv(condx, c(cuniqueid, condid))

  if (adj == "none") {
    setkeyv(condx, c(cuniqueid, condid))
  } else {
    
    if (popType %in% c("ALL", "VOL", "CURR")) {
      adjfacdata <- getadjfactorVOL(adj="plot", 
                        condx = condx, 
                        treex = treef, 
                        seedx = seedf, 
                        cuniqueid = cuniqueid, 
                        condid = condid,
                        areavar = areavar, 
                        areawt = areawt,
                        tpropvars = tpropvars
                        )
      condx <- adjfacdata$condx
      treef <- adjfacdata$treex
      seedf <- adjfacdata$seedx
      varadjlst <- adjfacdata$varadjlst
      areawtnm <- adjfacdata$areawtnm
    }
  }
 

  ###################################################################################
  ## Return population data objects
  ###################################################################################
  estvar.area <- ifelse(adj == "none", "CONDPROP_UNADJ", "CONDPROP_ADJ")

  if (!is.null(SAdoms)) {
    returnlst$SAdomsdf <- sf::st_drop_geometry(SAdoms)
  }
  if (!is.null(SAdoms) && is.null(smallbnd)) {
    if (!"AOI" %in% names(SAdoms)) {
      stop("missing AOI attribute in SAdoms")
    }
    smallbnd <- SAdoms[SAdoms$AOI == 1, ]
  } else {
    smallbnd <- pcheck.spatial(layer=smallbnd, caption="smallbnd")
    if (!"AOI" %in% names(smallbnd)) {
      smallbnd$AOI <- 1
    }
  }  
  if (is.null(smallbnd.domain)) {
    if ("DOMAIN" %in% names(smallbnd)) {
      smallbnd.domain <- "DOMAIN"
    } else if (length(names(sf::st_drop_geometry(smallbnd))) == 1) {
      smallbnd.domain <- names(sf::st_drop_geometry(smallbnd))
    } else {
      stop("must include smallbnd.domain for smallbnd")
    }
  } 
  returnlst$smallbnd <- smallbnd
  returnlst$smallbnd.domain <- smallbnd.domain


  ## Add xy attributes to pltassgnx
  ###############################################################################
  if (addxy) {
    if (is.null(smallbnd)) {
      message("need smallbnd for addxy")
      stop()
    }
    spxycompare <- crsCompare(spxy, smallbnd, nolonglat=TRUE)
    spxy <- spxycompare$x
    smallbnd <- spxycompare$y
    
    if (is.null(xy.uniqueid)) {
      if (pltassgnid %in% names(spxy)) {
        xy.uniqueid <- pltassgnid
      }
    }
    if (!is.null(xy.uniqueid) && xy.uniqueid %in% names(spxy)) {
      xy.coords <- data.frame(sf::st_drop_geometry(spxy[, xy.uniqueid]), sf::st_coordinates(spxy))
      pltassgnx <- merge(pltassgnx, xy.coords, by.x=pltassgnid, by.y=xy.uniqueid)
    } else {
      message(pltassgnid, " not in spxy names... cannot merge")
      stop()
    }
  }

  ###################################################################################
  ## Add new variables to pltcondx for estimation
  ###################################################################################
  ## Get order of pltcondx columns
  pltcondxcols <- names(pltcondx)
  newcols <- {}
  
  if (!"LANDSTATUSCD" %in% names(pltcondx) && "LANDSTATUSCD" %in% names(pltcondx)) {
    ## Add LANDSTATUSCD based on the following lookup table
    LANDSTATUSlut <- data.frame(LANDSTATUS = c(101:108, 111:117),
                    LANDSTATUSCD = c(rep(1, 6), rep(2, 2), rep(3, 6), 4),
                    LANDSTATUSNM = c(rep("Timberland", 6), 
                                     rep("Other forestland", 2), 
                                     rep("Reserved productive forestland", 6),
                                         "Reserved other forestland"))
    pltcondx$LANDSTATUS <- with(pltcondx, COND_STATUS_CD * 100 + RESERVCD * 10 + SITECLCD)
    pltcondx <- merge(pltcondx, LANDSTATUSlut, by="LANDSTATUS", all.x=TRUE)
    pltcondx$LANDSTATUS <- NULL
    newcols <- c("LANDSTATUSCD", "LANDSTATUSNM")
  }

  if (!"FORTYPGRPCD" %in% names(pltcondx) && "FORTYPCD" %in% names(pltcondx)) {
    ## Add FORTYPGRPCD to pltcondx if not already in dataset
    #pltcondx <- addFORTYPGRPCD(pltcondx)
    ref_fortyp <- ref_codes[ref_codes$VARIABLE == "FORTYPCD", c("VALUE", "GROUPCD")]
    names(ref_fortyp) <- c("FORTYPCD", "FORTYPGRPCD")
    pltcondx <- merge(pltcondx, ref_fortyp, by="FORTYPCD", all.x=TRUE)
    newcols <- c(newcols, "FORTYPGRPCD")
  }

  if (!"DSTRBGRP" %in% names(pltcondx) && "DSTRBCD1" %in% names(pltcondx)) {
    ## Add FORTYPGRPCD to pltcondx if not already in dataset
    #pltcondx <- addFORTYPGRPCD(pltcondx)
    ref_dstrbcd <- ref_codes[ref_codes$VARIABLE == "DSTRBCD", c("VALUE", "GROUPCD")]
    names(ref_dstrbcd) <- c("DSTRBCD1", "DSTRBGRP")
    pltcondx <- merge(pltcondx, ref_dstrbcd, by="DSTRBCD1", all.x=TRUE)
    newcols <- c(newcols, "DSTRBGRP")
  }
   
  ## Move new columns to end of table
  setcolorder(pltcondx, c(pltcondxcols, newcols))

  ## Subset condx
  condx <- condx[, c("DOMAIN", cuniqueid, condid, "CONDPROP_ADJ"), with=FALSE]
  
  ## Build list of data to return
  ###################################################################################
  returnlst <- append(returnlst, list(condx=condx, pltcondx=pltcondx, 
             pltassgnx = pltassgnx, pltassgnid = pltassgnid,
             cuniqueid=cuniqueid, condid=condid, ACI.filter=ACI.filter, 
             dunitarea=dunitarea, areavar=areavar, areaunits=areaunits, 
             dunitvar=dunitvar, dunitlut=data.table(dunitlut), 
             plotsampcnt=plotsampcnt, condsampcnt=condsampcnt, 
             states=states, invyrs=invyrs, estvar.area=estvar.area, adj=adj))

  if (!is.null(treef)) {
    returnlst$treex <- as.data.frame(treef)
    returnlst$tuniqueid <- tuniqueid
    returnlst$adjtree <- adjtree
  }
  if (!is.null(seedf)) {
    returnlst$seedx <- as.data.frame(seedf)
  }
  returnlst$prednames <- prednames
  returnlst$predfac <- predfac
  returnlst$largebnd.unique <- largebnd.unique


  ## Save list object
  ##################################################################
  if (saveobj) {
    if (getext(objfn) == "rds") {
      message("saving list object to: ", objfn)
      saveRDS(returnlst, objfn)
    } else if (getext(objfn) == "rda") {
      message("saving list object to: ", objfn)
      save(returnlst, objfn)
    } else {
      message("invalid object name... must end in: ", toString(c("rds", "rda")))
    } 
  } 

  ## Save data frames
  ##################################################################
  if (savedata) {
    datExportData(condx, 
          savedata_opts=list(outfolder=outfolder, 
                              out_fmt=out_fmt, 
                              out_dsn=out_dsn, 
                              out_layer="condx",
                              outfn.pre=outfn.pre, 
                              outfn.date=outfn.date, 
                              overwrite_layer=overwrite_layer,
                              append_layer=append_layer,
                              add_layer=TRUE))
    datExportData(pltcondx, 
          savedata_opts=list(outfolder=outfolder, 
                              out_fmt=out_fmt, 
                              out_dsn=out_dsn, 
                              out_layer="pltcondx",
                              outfn.pre=outfn.pre, 
                              outfn.date=outfn.date, 
                              overwrite_layer=overwrite_layer,
                              append_layer=append_layer,
                              add_layer=TRUE))

    if (!is.null(treef)) {
      datExportData(treef, 
            savedata_opts=list(outfolder=outfolder, 
                                out_fmt=out_fmt, 
                                out_dsn=out_dsn, 
                                out_layer="treex",
                                outfn.pre=outfn.pre, 
                                outfn.date=outfn.date, 
                                overwrite_layer=overwrite_layer,
                                append_layer=append_layer,
                                add_layer=TRUE))
    }
    if (!is.null(seedf)) {
      datExportData(seedf, 
            savedata_opts=list(outfolder=outfolder, 
                                out_fmt=out_fmt, 
                                out_dsn=out_dsn, 
                                out_layer="seedx",
                                outfn.pre=outfn.pre, 
                                outfn.date=outfn.date, 
                                overwrite_layer=overwrite_layer,
                                append_layer=append_layer,
                                add_layer=TRUE))
    }
    
    datExportData(pltassgnx, 
          savedata_opts=list(outfolder=outfolder, 
                              out_fmt=out_fmt, 
                              out_dsn=out_dsn, 
                              out_layer="pltassgn",
                              outfn.pre=outfn.pre, 
                              outfn.date=outfn.date, 
                              overwrite_layer=overwrite_layer,
                              append_layer=append_layer,
                              add_layer=TRUE))
    datExportData(dunitarea, 
          savedata_opts=list(outfolder=outfolder, 
                              out_fmt=out_fmt, 
                              out_dsn=out_dsn, 
                              out_layer="dunitarea",
                              outfn.pre=outfn.pre, 
                              outfn.date=outfn.date, 
                              overwrite_layer=overwrite_layer,
                              append_layer=append_layer,
                              add_layer=TRUE))
    datExportData(dunitlut, 
          savedata_opts=list(outfolder=outfolder, 
                              out_fmt=out_fmt, 
                              out_dsn=out_dsn, 
                              out_layer="dunitlut",
                              outfn.pre=outfn.pre, 
                              outfn.date=outfn.date, 
                              overwrite_layer=overwrite_layer,
                              append_layer=append_layer,
                              add_layer=TRUE))
  }

  if (returndata) {
    return(returnlst)
  } 
  rm(returnlst)
  # gc()
}
