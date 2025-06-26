#' Prepare Cox Model Output for Plotting
#'
#' Creates a tidy data frame from two Cox model outputs, including cumulative hazard estimates, standard errors, and confidence intervals.
#'
#' @param cox.out1 Output from a Cox model (e.g., `survfit`) for event 1.
#' @param cox.out2 Output from a Cox model for event 2.
#'
#' @return A `data.frame` containing combined cumulative hazard data for both transitions.
#' @export makeCox2eventsDf
makeCox2eventsDf <- function(cox.out1, cox.out2){

  df1 <- data.frame(
    time = cox.out1$time,
    estimate = cox.out1$cumhaz, # Cumulative hazard estimate
    std.err = cox.out1$std.err,
    n.risk = cox.out1$n.risk,
    n.event = cox.out1$n.event,
    Transition = "0 1"
  )

  df2 <- data.frame(
    time = cox.out2$time,
    estimate = cox.out2$cumhaz, # Cumulative hazard estimate
    std.err = cox.out2$std.err,
    n.risk = cox.out2$n.risk,
    n.event = cox.out2$n.event,
    Transition = "0 2"
  )

  cox_df <- bind_rows(df1,df2)

  # point-wise CI from variance
  cox_df$conf.low <- cox_df$estimate - 1.96 * cox_df$std.err # For 95% CI, critical value in normal distribution = 1.96
  cox_df$conf.high <- cox_df$estimate + 1.96 * cox_df$std.err

  return(cox_df)

}

#' Plot Cause-Specific Cumulative Hazards (Semi-Parametric)
#'
#' Plots the cumulative hazard estimates from two Cox model outputs using ggplot2, including 95% confidence intervals.
#'
#' @param cox_e1 Output from a Cox model (e.g., `survfit`) for event 1.
#' @param cox_e2 Output from a Cox model for event 2.
#' @param cols Vector of colors for each transition line. Default is `c("black", "darkgrey")`.
#'
#' @return A `ggplot` object visualizing the cumulative hazard functions.
#' @import ggplot2
#' @export plot_AbsRiskCox2Events
plot_AbsRiskCox2Events <- function(cox_e1,cox_e2, cols = default_cols){

  default_cols <- c("black", "darkgrey")

  cox_df <- makeCox2eventsDf(cox_e1,cox_e2)

  p <- ggplot(data = cox_df, aes(x = .data$time, y = .data$estimate, colour = .data$Transition)) +
    geom_line()  +
    geom_ribbon(
      aes(
        ymin = .data$conf.low,
        ymax = .data$conf.high,
        colour = .data$Transition
      ),
      fill = NA,
      linetype = "dotted"
    ) +
    scale_color_manual(values = cols) +
    ylab("Cause-specific cumulative hazard\n(semi-parametric)") +
    theme_classic()

  return(p)
}

#' Prepare Non-Parametric Cumulative Hazard Data
#'
#' Converts cumulative hazard estimates (e.g., from Aalen-Johansen estimator) into a tidy format with confidence intervals.
#'
#' @param cumhaz.out A list of two data frames for the cumulative hazard of each transition, typically output from `prodlim::prodlim`.
#'
#' @return A `data.frame` with hazard estimates and confidence bounds.
#' @import tidyr
#' @export makeCumhazDf
makeCumhazDf <- function(cumhaz.out){

  df1 <- cumhaz.out[[1]]
  colnames(df1)[1] <- "hazard"
  df1$transition <- names(cumhaz.out)[1]

  df2 <- cumhaz.out[[2]]
  colnames(df2)[1] <- "hazard"
  df2$transition <- names(cumhaz.out)[2]

  cumhaz_df <- bind_rows(df1,df2)

  # point-wise CI from variance given in eq. (4.1.6) of Andersen et al. (1993)
  # Andersen, P.K., Borgan, O., Gill, R.D. and Keiding, N. (1993). Statistical models based on counting processes. Springer Series in Statistics. New York, NY: Springer.
  cumhaz_df$lower <- cumhaz_df$hazard - 1.96 * sqrt(cumhaz_df$var.aalen)
  cumhaz_df$upper <- cumhaz_df$hazard + 1.96 * sqrt(cumhaz_df$var.aalen)

  return(cumhaz_df)

}

#' Plot Cause-Specific Cumulative Hazards (Non-Parametric)
#'
#' Generates a non-parametric plot of cause-specific cumulative hazards with pointwise 95% confidence intervals.
#'
#' @param comhaz.out List output of cumulative hazard estimates from a non-parametric method (e.g., `prodlim` or `survfit`).
#' @param cols Vector of colors for each transition. Default is `c("black", "darkgrey")`.
#'
#' @return A `ggplot` object.
#' @import ggplot2
#' @export plot_AbsRiskNA2Events
plot_AbsRiskNA2Events <- function(comhaz.out, cols = default_cols){

  default_cols <- c("black", "darkgrey")

  cumhaz_df <- makeCumhazDf(comhaz.out)

  p <- ggplot(cumhaz_df, aes(x = .data$time, y = .data$hazard, color = .data$transition)) +
    geom_line() +
    geom_ribbon(
      aes(ymin = .data$lower, ymax = .data$upper, color = .data$transition),
      fill = NA,
      linetype = "dotted"
    ) +
    scale_color_manual(values = cols) +
    labs(x = "Time", y = "Cause-specific cumulative hazard\n(non-parametric)", color = "Transition") +
    theme_classic()

  return(p)
}

#' Prepare Non-Parametric CIF Estimates
#'
#' Converts the output from a non-parametric CIF estimation into a tidy format for plotting.
#'
#' @param transprob.out Output from `etm::etm` or similar, representing transition probabilities.
#'
#' @return A `data.frame` with CIF estimates and confidence intervals.
#' @import tidyr
#' @export makeCifNpDf
makeCifNpDf <- function(transprob.out) {

  transprob.out <- summary(transprob.out)[[1]] # Gives CI in df format

  df1 <- transprob.out[[1]]
  df1$transition <- names(transprob.out)[1]

  df2 <- transprob.out[[2]]
  df2$transition <- names(transprob.out)[2]

  cif_df <- bind_rows(df1,df2)

  cif_df$transition <- gsub("CIF","0", cif_df$transition)

  return(cif_df)
}

#' Plot Non-Parametric Transition Probabilities (Aalen-Johansen)
#'
#' Generates a plot of non-parametric Aalen-Johansen estimates of transition probabilities over time.
#'
#' @param transprob Output from a multistate non-parametric estimator such as `etm::etm`.
#' @param cols Vector of colors for each transition. Default is `c("black", "darkgrey")`.
#'
#' @return A `ggplot` object.
#' @import ggplot2
#' @export plot_AbsRiskNpaj2Events
plot_AbsRiskNpaj2Events <- function(transprob, cols = default_cols){

  default_cols <- c("black", "darkgrey")

  cif_df <- makeCifNpDf(transprob)

  p <- ggplot(cif_df,
               aes(
                 x = .data$time,
                 y = .data$P,
                 color = .data$transition
               )) +
    geom_line() +
    geom_ribbon(
      aes(ymin = .data$lower, ymax = .data$upper, color = .data$transition),
      fill = NA,
      linetype = "dotted"
    ) +
    scale_color_manual(values = cols) +
    labs(x = "Time", y = "Transition probability\n(non-parametric)", color = "Transition") +
    theme_classic()

  return(p)

}

#' Prepare Semi-Parametric CIF Estimates
#'
#' Converts semi-parametric transition probability output into a tidy format for plotting.
#'
#' @param transprob.out A list with time-varying probabilities for each state, typically from `mstate` or similar.
#'
#' @return A `data.frame` suitable for plotting.
#' @import dplyr
#' @import tidyr
#' @export makeCifSpDf
makeCifSpDf <- function(transprob.out) {

  cif_df <- transprob.out[[1]] %>%
    dplyr::select(-c("pstate1")) %>%
    tidyr::pivot_longer(cols = c("pstate2","pstate3"), names_to = "transition", values_to = "P")

  cif_df$transition <- gsub("pstate2","0 1", cif_df$transition)
  cif_df$transition <- gsub("pstate3","0 2", cif_df$transition)

  return(cif_df)
}

#' Plot Semi-Parametric Transition Probabilities (Aalen-Johansen)
#'
#' Plots semi-parametric transition probability estimates over time using a Cox-Aalen-Johansen estimator.
#'
#' @param transprob Output from a semi-parametric transition probability estimator.
#' @param cols Vector of colors for each transition. Default is `c("black", "darkgrey")`.
#'
#' @return A `ggplot` object.
#' @import ggplot2
#' @export plot_AbsRiskSpaj2Events
plot_AbsRiskSpaj2Events <- function(transprob, cols = default_cols){

  default_cols <- c("black", "darkgrey")

  cif_df <- makeCifSpDf(transprob)

  p <- ggplot(cif_df,
              aes(
                x = .data$time,
                y = .data$P,
                color = .data$transition
              )) +
    geom_line() +
    # geom_ribbon(
    #   aes(ymin = lower, ymax = upper, color = transition),
    #   fill = NA,
    #   linetype = "dotted"
    # ) +
    scale_color_manual(values = cols) +
    labs(x = "Time", y = "Transition probability\n(semi-parametric)", color = "Transition") +
    theme_classic()

  return(p)

}









