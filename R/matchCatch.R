#' Match observed catch, yield and production with selectable sigmoid type
#'
#' This function adjusts the gear‐selectivity and mortality parameters for one
#' species so that a steady‐state model reproduces the observed catch size
#' distribution, the observed yield and the observed production, if available.
#'
#' @param params A `MizerParams` object
#' @param species The species for which to match the catch. Optional. By default
#'   all target species are selected. A vector of species names, or a numeric
#'   vector with the species indices, or a logical vector indicating for each
#'   species whether it is to be selected (TRUE) or not.
#' @param catch A data frame containing the observed binned catch data. It must
#'   contain the following columns:
#'   * `length`: The start of each bin.
#'   * `dl`: The width of each bin.
#'   * `count`: The observed count for each bin.
#' @param sel_type Character; either `"double"` (default) for dome-shaped
#'   selectivity, or `"single"` for a monotonic logistic.
#' @param lambda The slope of the community spectrum. Default is 2.05.
#' @param yield_lambda A parameter that controls the strength of the penalty for
#'   deviation from the observed yield.
#' @param production_lambda A parameter that controls the strength of the penalty
#'   for deviation from the observed production.
#'
#' @return A MizerParams object with the adjusted external mortality, gear
#'   selectivity, catchability and steady‐state spectrum for the selected
#'   species.
#' @family match functions
#' @examples
#' params <- matchCatch(celtic_params, species = "Hake", catch = celtic_catch)
#' plot_catch(params, species = "Hake", catch = celtic_catch)
#' @export
matchCatch <- function(params,
                       species = NULL,
                       catch,
                       sel_type = c("double", "single"),
                       lambda = 2.05,
                       yield_lambda = 1,
                       production_lambda = 1) {
    ## Handle catch/count column
    if (!"count" %in% names(catch) && "catch" %in% names(catch)) {
        catch$count <- catch$catch
    }

    ## Validate inputs
    species <- valid_species_arg(params, species = species, error_on_empty = TRUE)
    params <- validParams(params)

    ## If multiple species, recurse
    if (length(species) > 1) {
        for (s in species) {
            params <- matchCatch(params,
                                 species            = s,
                                 catch              = catch,
                                 sel_type           = sel_type,
                                 lambda             = lambda,
                                 yield_lambda       = yield_lambda,
                                 production_lambda  = production_lambda)
        }
        return(params)
    }

    ## Prepare data including sel_type flag
    sel_type           <- match.arg(sel_type)
    use_double_sigmoid <- as.integer(sel_type == "double")

    data <- prepare_data(params,
                         species            = species,
                         catch               = catch,
                         yield_lambda        = yield_lambda,
                         production_lambda   = production_lambda)
    if (is.null(data)) {
        warning(species, " cannot be matched because neither catches nor production are given.")
        return(params)
    }
    data$use_double <- use_double_sigmoid

    ## Extract species and gear params
    sp <- species_params(params)
    gp <- gear_params(params)
    sps <- sp[sp$species == species, ]
    gps <- gp[gp$species == species, ]

    ## Determine mu_mat
    mat_idx    <- sum(params@w < sps$w_mat)
    w_mat      <- params@w[mat_idx]
    g_mat      <- getEReproAndGrowth(params)[sp$species == species, mat_idx]
    mu_mat_max <- g_mat / w_mat * (lambda - sps$n)
    mu_mat     <- if (is.na(sps$mu_mat)) ext_mort(params)[sp$species == species, mat_idx] else sps$mu_mat

    ## Initial parameters for optimization
    initial_params <- list(
        l50          = gps$l50,
        ratio        = gps$l25 / gps$l50,
        d50          = gps$l50_right - gps$l50,
        mu_mat       = mu_mat,
        catchability = pmax(gps$catchability, 1e-8),
        r_right      = gps$l25_right / gps$l50_right
    )

    ## Bounds
    default_bounds <- list(
        l50          = c(5,    Inf),
        ratio        = c(0.1,  0.8),
        d50          = c(5,   80),
        mu_mat       = c(0.2,  mu_mat_max),
        catchability = c(1e-8, Inf),
        r_right      = c(1.3,  4)
    )
    lower_bounds <- sapply(default_bounds, `[`, 1)
    upper_bounds <- sapply(default_bounds, `[`, 2)

    ## If single sigmoid, fix descending params
    map <- list()
    if (!use_double_sigmoid) {
        initial_params["d50"]     <- default_bounds$d50[1]
        initial_params["r_right"] <- default_bounds$r_right[1]
        map$d50     <- factor(NA)
        map$r_right <- factor(NA)
        keep <- !names(lower_bounds) %in% c("d50", "r_right")
        lower_bounds <- lower_bounds[keep]
        upper_bounds <- upper_bounds[keep]
    }

    ## Lock selectivity if no count data
    if (!data$use_counts) {
        map$l50    <- factor(NA)
        map$ratio  <- factor(NA)
        map$d50    <- factor(NA)
        map$r_right<- factor(NA)
        keep <- !names(lower_bounds) %in% c("l50", "ratio", "d50", "r_right")
        lower_bounds <- lower_bounds[keep]
        upper_bounds <- upper_bounds[keep]
    }

    ## Lock catchability if yield not matched
    if (data$yield_lambda == 0) {
        map$catchability <- factor(NA)
        keep <- names(lower_bounds) != "catchability"
        lower_bounds <- lower_bounds[keep]
        upper_bounds <- upper_bounds[keep]
    }

    ## Build and run optimization
    obj <- TMB::MakeADFun(data       = data,
                          parameters = initial_params,
                          map        = map,
                          DLL        = "mizerEcopath",
                          silent     = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  lower = lower_bounds,
                  upper = upper_bounds)

    ## Update params with optimal values
    w_select      <- w(params) %in% data$w
    params_opt    <- update_params(params,
                                   species,
                                   opt$par,
                                   data,
                                   w_select)

    return(params_opt)
}
