#!/usr/bin/Rscript
# ASWKolaChokeGameSolver.R
# 2023-03-27
# curtis
# This is a one-line description of the file.

# argparser: A package for handling command line arguments
if (!suppressPackageStartupMessages(require("argparser"))) {
  install.packages("argparser")
  require("argparser")
}

# FUNCTIONS --------------------------------------------------------------------

#' Linear Optimizer Coefficients for Kola Scenario ASW Game
#'
#' Obtains Kola ASW scenario game coefficients that can be used for linear
#' programming problem
#'
#' The coefficients returned are probabilities that the arrangement would result
#' in the detection of an intruding target submarine (i.e. an SSBN in this
#' scenario). The scenario includes three waypoints: near the Kola peninsula,
#' transit through the Iceland-UK gap, and transit through the Greenland-Iceland
#' gap. Ships and ASW airplanes can be allocated to the Iceland gaps, as well as
#' submarines, but only submarines can wait outside the Kola peninsula, as ships
#' or planes patrolling the area would be intruding into Soviet waters and
#' airspace, provoking an incident. Nevertheless the Soviets may have delousing
#' assets that chase away submarines waiting for the target to exit the harbor.
#'
#' The mathematics for constructing the coefficients are given in the source
#' code.
#'
#' @param ships Number of ASW ships
#' @param subs Number of ASW submarines
#' @param planes Number of ASW planes
#' @param ships_green_ice Number of ASW ships allocated to Greenland-Iceland gap
#'                        (rest in Greenland-UK gap)
#' @param planes_green_ice Number of ASW planes allocated to Greenland-Iceland
#'                         gap (rest to Greenland-UK gap)
#' @param subs_green_ice Number of ASW submarines allocated to Greenland-Iceland
#'                       gap
#' @param subs_kola ASW submarines patroling near Kola peninsula
#' @param ships_sweep Effective sweep width of ships with constant speed
#' @param subs_sweep Effective sweep width of submarines with constant speed
#' @param planes_sweep Effective sweep width of planes with constant speed
#' @param kola_width Width of the Kola patrol area
#' @param green_ice_width Width of the Greenland-Iceland gap
#' @param ice_uk_width Width of the Iceland-UK gap
#' @param search_const Multiplicative constant with all effective widths
#' @param delouse_total_sweep Effective sweep width of all delousing assets
#' @return Vector with coefficients for the two routes the target could take,
#'         which are probabilities of detection
#' @examples
#' sub_strategy_kola_coef()
sub_strategy_kola_coef <- function(ships = 0,
                                   subs = 0,
                                   planes = 0,
                                   ships_green_ice = 0,
                                   planes_green_ice = 0,
                                   subs_green_ice = 0,
                                   subs_kola = 0,
                                   ships_sweep = 1,
                                   subs_sweep = 3,
                                   planes_sweep = 15,
                                   kola_width = 1/4,
                                   green_ice_width = 1,
                                   ice_uk_width = 2,
                                   search_const = 1,
                                   delouse_total_sweep = 0) {
  ships_ice_uk <- ships - ships_green_ice
  planes_ice_uk <- planes - planes_green_ice
  subs_ice_uk <- subs - subs_green_ice - subs_kola
  stopifnot(ships_ice_uk >= 0 && planes_ice_uk >= 0 && subs_ice_uk >= 0)
  stopifnot(ships_sweep >= 0 && subs_sweep >= 0 && planes_sweep >= 0 &&
            search_const > 0 && delouse_total_sweep >= 0)
  stopifnot(ships >= 0 && subs >= 0 && planes >= 0)
  stopifnot(kola_width >0 && ice_uk_width > 0 && green_ice_width >= 0)

  ice_uk_detect_prob <- 1 - exp(- search_const * ice_uk_width *
    (ships_ice_uk * ships_sweep + subs_ice_uk * subs_sweep +
      planes_ice_uk * planes_sweep))
  green_ice_detect_prob <- 1 - exp(- search_const * green_ice_width *
    (ships_green_ice * ships_sweep + subs_green_ice * subs_sweep +
      planes_green_ice * planes_sweep))
  delouse_prob <- 1 - exp(- search_const * kola_width * delouse_total_sweep)
  kola_detect_prob <- 1 - (delouse_prob + (1 - delouse_prob) *
    exp(-search_const * kola_width * subs_sweep))^subs_kola

  1 - c("green_ice" = 1 - green_ice_detect_prob,
        "ice_uk" = 1 - ice_uk_detect_prob) * (1 - kola_detect_prob)
}

#' @param max_time_on_station Maximum time on station for target
#' @param time_lost_green_ice Time lost for target when taking Greenland-Iceland
#'        route

# EXECUTABLE SCRIPT MAIN FUNCTIONALITY -----------------------------------------

main <- function() {
}

# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION -------------------------

if (sys.nframe() == 0) {
  p <- arg_parser("This is a template for executable R scripts.")
  p <- add_argument(p, "foo", type = "character", nargs = 1,
                    help = "A positional command-line argument")
  p <- add_argument(p, "--bar", type = "integer", default = 0,
                    nargs = 1, help = "A command-line option")
  p <- add_argument(p, "--baz", flag = TRUE, help = "A command-line flag")

  cl_args <- parse_args(p)
  cl_args <- cl_args[!(names(cl_args) %in% c("help", "opts"))]
  if (any(sapply(cl_args, is.na))) {
    # User did not specify all inputs; print help message
    print(p)
    cat("\n\nNot all needed inputs were given.\n")
    quit()
  }

  do.call(main, cl_args[2:length(cl_args)])
}

