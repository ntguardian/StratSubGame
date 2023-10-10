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

suppressPackageStartupMessages(library(lpSolve))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))

# COMMAND LINE INTERFACE -------------------------------------------------------

p <- arg_parser("Calculate Stackelberg strategic solution to Kola peninsula theater anti-submarine warfare (ASW) scenario")
p <- add_argument(p, "--ships", type = "integer", nargs = 1, default = 0,
                  help = "Numbver of ASW ships")
p <- add_argument(p, "--subs", type = "integer", nargs = 1, default = 0,
                  help = "Numbver of ASW submarines")
p <- add_argument(p, "--planes", type = "integer", nargs = 1, default = 0,
                  help = "Numbver of ASW airplanes")
p <- add_argument(p, "--shipssweep", type = "numeric", default = 1,
                  nargs = 1, help = "Sweep width of ASW ships")
p <- add_argument(p, "--subssweep", type = "numeric", default = 1,
                  nargs = 1, help = "Sweep width of ASW submarines")
p <- add_argument(p, "--planessweep", type = "numeric", default = 15,
                  nargs = 1, help = "Sweep width of ASW airplanes")
p <- add_argument(p, "--kolawidth", type = "numeric", default = 1/4,
                  nargs = 1, help = "Width of the Kola peninsula gap")
p <- add_argument(p, "--greenicewidth", type = "numeric", default = 1,
                  nargs = 1, help = "Width of the Greenland-Iceland gap")
p <- add_argument(p, "--iceukwidth", type = "numeric", default = 2,
                  nargs = 1, help = "Width of the Iceland-UK gap")
p <- add_argument(p, "--searchconst", type = "numeric", default = 1,
                  nargs = 1, help = "Multiplicative constant with all effective widths")
p <- add_argument(p, "--delouse", type = "numeric", default = 1,
                  nargs = 1, help = "Effective sweep width of all delousing assets")
p <- add_argument(p, "--maxtime", type = "numeric", default = 5,
                  nargs = 1, help = "Maximum time on station for intruder submarine")
p <- add_argument(p, "--greenicetimeloss", type = "numeric", default = 1,
                  nargs = 1, help = "Time lost by transiting through the Greenland-Iceland gap relative to traveling through the Iceland-UK gap")

# FUNCTIONS --------------------------------------------------------------------

#' ASW Linear Optimizer Coefficients for Kola Scenario ASW Game
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
#' asw_strategy_kola_coef()
asw_strategy_kola_coef <- function(ships = 0,
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

#' Sub Payoff Linear Optimizer Coefficients for Kola Scenario ASW Game
#'
#' Payoff column for the submarine in the ASW Kola peninsula scenario
#'
#' This returns the payoff function of the submarine. It requires running
#' \code{\link{asw_intercept_probs}} first to generate the interception
#' probabilities. Once done this will appropriately account for time on station
#' due to the route taken.
#'
#' @param asw_intercept_probs A vector returned by \code{asw_strategy_kola_coef}
#' @param max_time_on_station Maximum time on station for target
#' @param time_lost_green_ice Time lost for target when taking Greenland-Iceland
#'                            route
#' @return Vector of payoff for submarine
#' @examples
#' sub_strategy_kola_coef(asw_strategy_kola_coef())
sub_strategy_kola_coef <- function(asw_intercept_probs,
                                   max_time_on_station = 5,
                                   time_lost_green_ice = 1) {
  sub_payoff <- max_time_on_station * (1 - asw_intercept_probs)
  sub_payoff - c("green_ice" = time_lost_green_ice, "ice_uk" = 0)
}

#' Generate ASW Strategies and Payoffs for Kola ASW Scenario
#' 
#' Generate a \code{\link[base]{matrix}} containing the strategies and the
#' payoffs for the ASW side of the Kola scenario
#'
#' The last columns contain the actual game matrix while the first columns are
#' descriptors of the strategies.
#'
#' @param ships,subs,planes Number of ships, subs, and planes for ASW force in
#'                          scenario
#' @param ... Parameters to pass to \code{\link{asw_strategy_kola_coef}}; note
#'            that  parameters like \code{ships_green_ice},
#'            \code{subs_green_ice}, and \code{subs_kola} are ignored
#' @return Matrix of ASW payoff (probability of interception depending on
#'         submarine route)
#' @examples
#' generate_asw_kola_strats(5, 4, 3)
#' @export
generate_asw_kola_strats <- function(ships, subs, planes, ...) {
  ship_counts <- 0:ships
  plane_counts <- 0:planes
  sub_count_prelim <- combn(subs + 2, 2)
  sub_counts <- cbind(sub_count_prelim[1,] - 1,
                      sub_count_prelim[2,] - sub_count_prelim[1,] - 1)
  colnames(sub_counts) <- c("subs_green_ice", "subs_kola")

  game_mat <- Reduce(rbind, lapply(ship_counts,
                                   FUN = function(ships_green_ice) {
    Reduce(rbind, lapply(plane_counts, FUN = function(planes_green_ice) {
      Reduce(rbind, apply(sub_counts, MARGIN = 1, FUN = function(sub_alloc) {
        c("ships_green_ice" = ships_green_ice,
          "ships_ice_uk" = ships - ships_green_ice,
          "planes_green_ice" = planes_green_ice,
          "planes_ice_uk" = planes - planes_green_ice,
          sub_alloc,
          "subs_ice_uk" = subs - sum(sub_alloc),
          asw_strategy_kola_coef(ships = ships,
                                 subs = subs,
                                 planes = planes,
                                 subs_green_ice = sub_alloc["subs_green_ice"],
                                 subs_kola = sub_alloc["subs_kola"],
                                 planes_green_ice = planes_green_ice,
                                 ships_green_ice = ships_green_ice,
                                 ...))
  }, simplify = FALSE))}))}))
  if (!is.matrix(game_mat)) {
    game_mat <- t(as.matrix(game_mat))
  }
  colnames(game_mat) <- gsub("green_ice.subs_green_ice", "pay_green_ice",
    gsub("ice_uk.subs_green_ice", "pay_ice_uk", colnames(game_mat)))
  rownames(game_mat) <- NULL
  game_mat
}

#' Generate Sub Strategies and Payoffs for Kola ASW Scenario
#'
#' Generate a \code{\link[base]{matrix}} containing the strategies and the
#' payoffs for the submarine side of the Kola scenario
#'
#' The last columns contain the actual game matrix while the first columns are
#' descriptors of the strategies.
#'
#' @param asw_strat_matrix Matrix returned by \code{\link{generate_asw_strats}}
#'                         containing the ASW detection probabilities
#' @param ... Additional parameters to pass to
#'            \code{\link{sub_strategy_kola_coef}}
#' @return Matrix of strategies with payoffs for the intruding submarine
#' @examples
#' generate_sub_kola_strats(generate_asw_kola_strats(5, 4, 3))
#' @export
generate_sub_kola_strats <- function(asw_strat_matrix, ...) {
  asw_probs <- asw_strat_matrix[, c("pay_green_ice", "pay_ice_uk"), drop = FALSE]
  colnames(asw_probs) <- c("green_ice", "ice_uk")
  sub_payoff <- apply(asw_probs, MARGIN = 1,
    FUN = function(asw_intercept_probs) {
      sub_strategy_kola_coef(asw_intercept_probs, ...)
  })
  rownames(sub_payoff) <- c("pay_green_ice", "pay_ice_uk")
  cbind(asw_strat_matrix[,
                         !(colnames(asw_strat_matrix) %in% c("pay_green_ice",
                                                             "pay_ice_uk")),
                         drop = FALSE],
    t(sub_payoff))
}

#' Find Optimal Submarine Strategy Given ASW Strategy
#'
#' Find the submarine's optimal response to a given ASW strategy
#'
#' This obtains the ASW force's optimal strategy given the strategy of the
#' submarine, a key step in obtaining the game's optimal solution.
#'
#' @param asw_coef Row of the ASW payoff matrix corresponding to the strategy
#'                 for which to solve the program
#' @param sub_strat_row Row of the submarine payoff matrix that corresponds to
#'                      the strategy given in \code{asw_coef}
#' @param sub_strat_matrix Matrix of the submarine payoff matrix
#' @return Value and strategic allocation for given strategies
#' @examples
#' asw_strat <- generate_asw_kola_strats(5, 4, 3)
#' sub_strat <- generate_sub_kola_strats(asw_strat)
#' strat_linear_prog_solve(asw_strat[, "pay_green_ice"],
#'                         sub_strat[, "pay_green_ice"],
#'                         sub_strat[, c("pay_green_ice", "pay_ice_uk")])
strat_linear_prog_solve <- function(asw_coef, sub_strat_col, sub_strat_matrix) {
  sub_coefs_matrix <- sub_strat_col - sub_strat_matrix
  constraint_mat <- rbind(rep(1, length(asw_coef)),
                          t(sub_coefs_matrix))
  constraint_type <- c("==", rep(">=", ncol(sub_coefs_matrix)))
  constraint_number <- c(1, rep(0, ncol(sub_coefs_matrix)))
  lp_obj <- lp(direction = "max",
     objective.in = asw_coef,
     const.mat = constraint_mat,
     const.dir = constraint_type,
     const.rhs = constraint_number
  )[c("objval", "solution", "status")]
  lp_obj
}

#' Solve Submarine Theater Chokepoint ASW Game
#'
#' Completely solve the theater ASW game
#'
#' This solves the game by running \code{\link{strat_linear_prog_solve}} for
#' every pure strategy available to the intruding submarine. The optimal
#' strategies for the ASW force are returned and the optimal strategy can be
#' selected.
#'
#' @param asw_strat_matrix \code{link[base]{matrix}} with game payoff matrix for
#'                         asw force
#' @param sub_strat_matrix \code{link[base]{matrix}} with game payoff matrix for
#'                         submarine force
#' @param game_cols \code{\link[base]{vector}} describing which columns of both
#'                  \code{asw_strat_matrix} and \code{sub_strat_matrix} contain
#'                  the actual game payoff matrix; the others are presumed to be
#'                  useful but not presently relevant information, and while not
#'                  a part of the solver, will be a part of the returned data
#'                  frame
#' @param full Return the full data frame; otherwise, only the rows with the
#'             best value for the ASW force
#' @return A list with the solution matrix (\code{solution_df}), an identifier
#'         for the optimal strategies (\code{solution_codes}), and the optimal
#'         strategy payoff (\code{solution_value}); see
#'         \code{\link{strat_linear_prog_solve}} for interpretation of number
#'         associated with \code{solution_codes}
#' @examples
#' asw_strat <- generate_asw_kola_strats(5, 4, 3)
#' sub_strat <- generate_sub_kola_strats(asw_strat)
#' asw_chokepoint_game_solution_df(asw_strat, sub_strat,
#'                                 c("pay_green_ice", "pay_ice_uk"))
#' @export
asw_chokepoint_game_solution_df <- function(asw_strat_matrix, sub_strat_matrix,
                                            game_cols = 1:ncol(asw_strat_matrix),
                                            full = FALSE) {
  if (is.character(game_cols)) {
    game_cols <- which(colnames(asw_strat_matrix) %in% game_cols)
  }
  not_game_cols <- which(!(1:ncol(asw_strat_matrix) %in% game_cols))
  asw_game_mat <- asw_strat_matrix[, game_cols, drop = FALSE]
  sub_game_mat <- sub_strat_matrix[, game_cols, drop = FALSE]

  solution_df <- as_tibble(asw_strat_matrix[, not_game_cols, drop = FALSE]) %>%
    bind_cols(as_tibble(asw_game_mat) %>%
      rename_with(.fn = ~ paste0(.x, "_asw"))) %>%
    bind_cols(as_tibble(sub_game_mat) %>%
      rename_with(.fn = ~ paste0(.x, "_sub")))

  solution_cols <- lapply(1:ncol(asw_game_mat), FUN = function(strat_col) {
      strat_linear_prog_solve(asw_game_mat[, strat_col],
        sub_game_mat[, strat_col],
        sub_game_mat)})
  game_solution_list <- map(solution_cols, ~ .x$solution)
  names(game_solution_list) <- paste0(colnames(asw_game_mat), "_optim")
  game_solution_df <- as_tibble(game_solution_list)
  solution_codes <- map_int(solution_cols, ~ .x$status)
  names(solution_codes) <- colnames(asw_game_mat)
  solution_value <- map_dbl(solution_cols, ~ .x$objval)
  names(solution_value) <- colnames(asw_game_mat)

  if (full) {
    solution_df <- bind_cols(solution_df, game_solution_df)
  } else {
    solution_df <- bind_cols(solution_df,
      game_solution_df[which.max(solution_value)])
    solution_codes <- solution_codes[which.max(solution_value)]
    solution_value <- max(solution_value)
  }

  list("solution_df" = solution_df,
       "solution_codes" = solution_codes,
       "solution_value" = solution_value)
}

# EXECUTABLE SCRIPT MAIN FUNCTIONALITY -----------------------------------------

main <- function(ships = 0,
                 subs = 0,
                 planes = 0,
                 shipssweep = 1,
                 subssweep = 3,
                 planessweep = 15,
                 kolawidth = 1/4,
                 greenicewidth = 1,
                 iceukwidth = 2,
                 searchconst = 1,
                 delouse = 0,
                 maxtime = 5,
                 greenicetimeloss = 1) {
  asw_strats <- generate_asw_kola_strats(
    ships = ships,
    subs = subs,
    planes = planes,
    ships_sweep = shipssweep,
    subs_sweep = subssweep,
    planes_sweep = planessweep,
    kola_width = kolawidth,
    ice_uk_width = iceukwidth,
    green_ice_width = greenicewidth,
    search_const = searchconst,
    delouse_total_sweep = delouse
  )
  sub_strats <- generate_sub_kola_strats(
    asw_strats,
    max_time_on_station = maxtime,
    time_lost_green_ice = greenicetimeloss
  )
  solution <- asw_chokepoint_game_solution_df(asw_strats, sub_strats,
                                              c("pay_green_ice", "pay_ice_uk"))
  solution_df <- solution$solution_df
  solution_col <- paste0(names(solution$solution_codes), "_optim")

  cat("Probability of interception:", solution$solution_value, "\n")
  cat("Sub transit path:", substr(names(solution$solution_codes), 5, 999), "\n")
}

# INTERFACE DEFINITION AND COMMAND LINE IMPLEMENTATION -------------------------

if (sys.nframe() == 0) {
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

