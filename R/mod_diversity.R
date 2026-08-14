#' diversity UI Function
#'
#' @description Alpha and beta diversity with design-aware guards. Every
#'   procedure here has a sample size below which its output is not merely
#'   under-powered but undefined, and those limits are enforced rather than
#'   left to the reader: classical MDS has at most `n - 1` axes, a
#'   two-dimensional NMDS configuration is saturated at small `n` so its stress
#'   is trivially zero, and a permutation test cannot reach `p < 0.05` unless
#'   more than 20 distinguishable label permutations exist.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_diversity_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        3,
        vpf_card(
          "Design",
          uiOutput(ns("group_ui")),
          uiOutput(ns("design_report"))
        ),
        vpf_card(
          "Alpha diversity",
          selectInput(ns("alpha_index"), "Index",
                      choices = c("Shannon" = "shannon",
                                  "Simpson (Gini-Simpson)" = "gini_simpson",
                                  "Inverse Simpson" = "inverse_simpson",
                                  "Observed richness" = "observed_richness",
                                  "Chao1" = "chao1"),
                      selected = "shannon"),
          uiOutput(ns("alpha_note"))
        ),
        vpf_card(
          "Beta diversity",
          selectInput(ns("beta_method"), "Dissimilarity",
                      choices = c("Bray-Curtis" = "bray",
                                  "Jaccard (presence/absence)" = "jaccard"),
                      selected = "bray"),
          selectInput(ns("ordination"), "Ordination",
                      choices = c("PCoA (classical MDS)" = "pcoa",
                                  "NMDS" = "nmds"),
                      selected = "pcoa"),
          numericInput(ns("permutations"), "PERMANOVA permutations",
                       value = 999, min = 99, max = 9999, step = 100)
        )
      ),
      column(
        9,
        uiOutput(ns("subset_warning")),
        tabsetPanel(
          id = ns("tabs"),
          tabPanel(
            "Alpha diversity",
            plotly::plotlyOutput(ns("plt_alpha"), height = "420px"),
            uiOutput(ns("alpha_test")),
            reactable::reactableOutput(ns("tbl_alpha"))
          ),
          tabPanel(
            "Ordination",
            plotly::plotlyOutput(ns("plt_ordination"), height = "460px"),
            uiOutput(ns("ordination_note"))
          ),
          tabPanel(
            "Distances",
            plotly::plotlyOutput(ns("plt_distance"), height = "420px"),
            vpf_caption(
              "Pairwise dissimilarity between samples. With two samples this is a",
              "single descriptive number: 0 means identical composition, 1 means",
              "no shared abundance."
            )
          ),
          tabPanel(
            "PERMANOVA",
            uiOutput(ns("permanova_out"))
          )
        )
      )
    )
  )
}

#' diversity Server Functions
#'
#' @param id Module id.
#' @param r_filter A list of reactives from [mod_vpfilter_server()].
#'
#' @noRd
mod_diversity_server <- function(id, r_filter) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$group_ui <- renderUI({
      groups <- vpf_group_candidates(r_filter$tse())
      if (length(groups) == 0) {
        return(vpf_notice(
          type = "warning",
          title = "No grouping variable",
          "Group comparisons and PERMANOVA need a colData column with at least",
          "two levels and at least two samples per level. Sample identifiers do",
          "not qualify: they create one group per sample. Upload a metadata",
          "table on the Data tab to supply one."
        ))
      }
      selectInput(ns("group"), "Grouping variable",
                  choices = c("(none)" = "", groups), selected = "")
    })

    active_group <- reactive({
      g <- input$group
      if (is.null(g) || !nzchar(g)) NULL else g
    })

    design <- reactive({
      vpf_design_check(r_filter$tse(), active_group(), r_filter$assay())
    })

    output$design_report <- renderUI({
      d <- design()
      if (d$n_samples == 0) {
        return(vpf_caption("Load a dataset first."))
      }
      sizes <- if (d$n_groups > 0) {
        paste0(" (", paste(paste0(names(d$group_sizes), ": ", d$group_sizes),
                           collapse = ", "), ")")
      } else {
        ""
      }
      tagList(
        vpf_caption(
          sprintf("%d samples, %d contigs, %d group(s)%s.",
                  d$n_samples, d$n_features, d$n_groups, sizes)
        ),
        tags$ul(
          style = "font-size: 0.85em; padding-left: 18px;",
          tags$li(paste0("Group tests: ", if (d$can_group_test) "available" else "disabled")),
          tags$li(paste0("PCoA (2 axes): ", if (d$can_pcoa) "available" else "disabled")),
          tags$li(paste0("NMDS (2 axes): ", if (d$can_nmds) "available" else "disabled")),
          tags$li(paste0("PERMANOVA: ", if (d$can_permanova) "available" else "disabled"))
        )
      )
    })

    output$subset_warning <- renderUI({
      tse <- r_filter$tse()
      raw <- r_filter$raw()
      if (is.null(tse) || is.null(raw)) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      removed <- nrow(raw) - nrow(tse)
      votes <- vpf_viral_votes(tse)
      vote_line <- if (!is.null(votes) && length(votes) > 0) {
        tab <- table(votes)
        singletons <- sum(votes <= 1, na.rm = TRUE)
        paste0(
          " Of the retained contigs, ",
          paste(sprintf("%d have %s vote(s)", as.integer(tab), names(tab)),
                collapse = ", "),
          ". The ", singletons,
          " admitted on a single detector contribute to these indices exactly as",
          " much as the fully supported ones; the Filter tab can exclude them."
        )
      } else {
        ""
      }
      vpf_notice(
        type = "info",
        title = "These numbers describe a selected subset",
        paste0(
          sprintf(paste("%d of %d contigs are currently retained; %d were removed by",
                        "the filters on the Filter tab, and that object is itself the",
                        "viral subset of the assembly. Selecting contigs by viral",
                        "evidence removes the low-scoring tail of the community, so",
                        "the diversity values below are comparable between samples",
                        "processed identically but are not absolute community",
                        "diversities. The Filter tab's audit lists which rule removed",
                        "what."),
                  nrow(tse), nrow(raw), removed),
          vote_line
        )
      )
    })

    output$alpha_note <- renderUI({
      d <- design()
      tse <- r_filter$tse()
      assay_name <- r_filter$assay()
      if (is.null(tse) || is.null(assay_name)) {
        return(NULL)
      }
      msgs <- list()
      if (identical(input$alpha_index, "chao1") && !d$can_chao1) {
        msgs <- c(msgs, list(vpf_notice(
          type = "danger",
          paste0("Chao1 is estimated from singleton and doubleton features and is ",
                 "only defined for raw integer counts. The assay '", assay_name,
                 "' does not hold raw counts.")
        )))
      }
      if (identical(input$alpha_index, "chao1") && d$can_chao1) {
        msgs <- c(msgs, list(vpf_notice(
          type = "warning",
          paste("Chao1 is exploratory for assembled viral contigs: assembly",
                "fragmentation, detection thresholds and multi-mapping violate the",
                "unseen-species model it assumes.")
        )))
      }
      if (identical(input$alpha_index, "observed_richness")) {
        msgs <- c(msgs, list(vpf_caption(
          "Observed richness counts detected contigs and therefore tracks",
          "sequencing depth and the coverage threshold as much as biology."
        )))
      }
      if (length(d$empty_libraries) > 0) {
        msgs <- c(msgs, list(vpf_notice(
          type = "warning",
          paste0("Samples with no signal after filtering are excluded: ",
                 paste(d$empty_libraries, collapse = ", "), ".")
        )))
      }
      if (length(msgs) == 0) NULL else do.call(tagList, msgs)
    })

    alpha_values <- reactive({
      tse <- r_filter$tse()
      assay_name <- r_filter$assay()
      index <- input$alpha_index %||% "shannon"
      d <- design()
      if (is.null(tse) || is.null(assay_name) || nrow(tse) == 0 || ncol(tse) == 0) {
        return(NULL)
      }
      if (identical(index, "chao1") && !d$can_chao1) {
        return(NULL)
      }
      keep <- !vpf_sample_names(tse) %in% d$empty_libraries
      if (!any(keep)) {
        return(NULL)
      }
      sub <- tse[, keep, drop = FALSE]
      res <- withProgress(message = "Computing alpha diversity", value = 0.5, {
        tryCatch(
          mia::addAlpha(sub, assay.type = assay_name, index = index, name = "alpha_value"),
          error = function(e) e
        )
      })
      if (inherits(res, "error")) {
        return(list(error = conditionMessage(res)))
      }
      cd <- SummarizedExperiment::colData(res)
      df <- data.frame(
        Sample = vpf_sample_names(res),
        Value = as.numeric(cd[["alpha_value"]]),
        stringsAsFactors = FALSE
      )
      grp <- active_group()
      if (!is.null(grp) && grp %in% colnames(cd)) {
        df$Group <- vpf_blank_na(cd[[grp]])
        df$Group[is.na(df$Group)] <- "not annotated"
      }
      df$Sample <- factor(df$Sample, levels = df$Sample)
      list(data = df, index = index, assay = assay_name)
    })

    output$plt_alpha <- plotly::renderPlotly({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(vpf_message_plot("Load a dataset on the Data tab first."))
      }
      if (nrow(tse) == 0) {
        return(vpf_message_plot("No contigs pass the current filters."))
      }
      av <- alpha_values()
      if (is.null(av)) {
        if (identical(input$alpha_index, "chao1")) {
          return(vpf_message_plot(paste(
            "Chao1 requires raw integer counts because it is estimated from",
            "singleton and doubleton features. Choose the counts assay on the",
            "Filter tab, or select a different index."
          )))
        }
        return(vpf_message_plot("No sample retains any signal after filtering."))
      }
      if (!is.null(av$error)) {
        return(vpf_message_plot(paste("Alpha diversity failed:", av$error)))
      }
      df <- av$data
      has_group <- "Group" %in% colnames(df)
      # Samples are always labelled: with two anonymous points a reader cannot
      # tell which sample is which.
      p <- if (has_group) {
        ggplot2::ggplot(df, ggplot2::aes(x = .data$Group, y = .data$Value,
                                         colour = .data$Group,
                                         text = paste0(.data$Sample, "<br>",
                                                       signif(.data$Value, 4)))) +
          ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.2) +
          ggplot2::geom_point(size = 3, position = ggplot2::position_jitter(width = 0.12, height = 0))
      } else {
        ggplot2::ggplot(df, ggplot2::aes(x = .data$Sample, y = .data$Value,
                                         text = paste0(.data$Sample, "<br>",
                                                       signif(.data$Value, 4)))) +
          ggplot2::geom_col(fill = "#4E79A7")
      }
      p <- p +
        ggplot2::labs(x = NULL, y = av$index, colour = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      vpf_ggplotly(p, tooltip = "text")
    })

    output$tbl_alpha <- reactable::renderReactable({
      av <- alpha_values()
      if (is.null(av) || !is.null(av$error)) {
        return(reactable::reactable(data.frame(
          Message = "No alpha-diversity values available for the current selection."
        )))
      }
      reactable::reactable(vpf_round_df(av$data), striped = TRUE, highlight = TRUE,
                           bordered = TRUE, pagination = FALSE)
    })

    output$alpha_test <- renderUI({
      d <- design()
      av <- alpha_values()
      if (is.null(av) || !is.null(av$error) || !"Group" %in% colnames(av$data)) {
        if (d$n_samples > 0 && d$n_groups < 2) {
          return(vpf_notice(
            type = "info",
            "No group test: a comparison needs a grouping variable with at least",
            "two levels."
          ))
        }
        return(NULL)
      }
      if (!d$can_group_test) {
        return(vpf_notice(
          type = "warning",
          title = "Group test disabled",
          sprintf(paste("Smallest group has %d sample(s). A group test requires at",
                        "least two biological samples per group; with one sample per",
                        "group there is no within-group variability and a p-value",
                        "would not be evidence."), d$min_group_size)
        ))
      }
      df <- av$data
      df$Group <- factor(df$Group)
      res <- tryCatch({
        if (nlevels(df$Group) == 2) {
          ht <- stats::wilcox.test(Value ~ Group, data = df, exact = FALSE)
          list(name = "Wilcoxon rank-sum test", p = ht$p.value)
        } else {
          ht <- stats::kruskal.test(Value ~ Group, data = df)
          list(name = "Kruskal-Wallis test", p = ht$p.value)
        }
      }, error = function(e) list(name = NA, p = NA, error = conditionMessage(e)))

      if (!is.null(res$error)) {
        return(vpf_notice(type = "danger", paste("Test failed:", res$error)))
      }

      warnings <- list()
      if (nlevels(df$Group) == 2 && !isTRUE(d$wilcoxon_can_reach_05)) {
        sz <- as.integer(d$group_sizes)
        warnings <- c(warnings, sprintf(
          paste("With group sizes %s the smallest attainable two-sided exact",
                "p-value is 2 / choose(%d, %d) = %.3f, so this test cannot reach",
                "p < 0.05 whatever the data show. Report the values themselves."),
          paste(sz, collapse = " vs "), sum(sz), sz[[1]],
          2 / choose(sum(sz), sz[[1]])
        ))
      }
      if (nlevels(df$Group) > 2 && !d$kw_asymptotic_ok) {
        warnings <- c(warnings, paste(
          "Group sizes are below five, so the asymptotic chi-square p-value of",
          "the Kruskal-Wallis test is unreliable. Treat it as descriptive."
        ))
      }
      tagList(
        vpf_notice(
          type = if (length(warnings) > 0) "warning" else "info",
          title = res$name,
          sprintf("p = %.4g", res$p),
          if (length(warnings) > 0) {
            tagList(tags$br(), lapply(warnings, tags$div))
          }
        )
      )
    })

    # -- Beta diversity ------------------------------------------------------
    dissimilarity <- reactive({
      tse <- r_filter$tse()
      assay_name <- r_filter$assay()
      d <- design()
      if (is.null(tse) || is.null(assay_name) || nrow(tse) == 0 || ncol(tse) < 2) {
        return(NULL)
      }
      keep <- !vpf_sample_names(tse) %in% d$empty_libraries
      if (sum(keep) < 2) {
        return(NULL)
      }
      sub <- tse[, keep, drop = FALSE]
      mat <- t(as.matrix(SummarizedExperiment::assay(sub, assay_name)))
      mat[is.na(mat)] <- 0
      rownames(mat) <- vpf_sample_names(sub)
      method <- input$beta_method %||% "bray"
      dist <- tryCatch({
        if (identical(method, "jaccard")) {
          vegan::vegdist(mat, method = "jaccard", binary = TRUE)
        } else {
          vegan::vegdist(mat, method = "bray")
        }
      }, error = function(e) e)
      if (inherits(dist, "error")) {
        return(list(error = conditionMessage(dist)))
      }
      if (anyNA(dist)) {
        return(list(error = paste(
          "The dissimilarity matrix contains missing values, which happens when a",
          "sample has no non-zero abundance. Loosen the filters or the coverage",
          "threshold."
        )))
      }
      list(dist = dist, method = method, samples = rownames(mat),
           groups = if (!is.null(active_group())) {
             vpf_blank_na(SummarizedExperiment::colData(sub)[[active_group()]])
           } else {
             NULL
           })
    })

    output$plt_ordination <- plotly::renderPlotly({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(vpf_message_plot("Load a dataset on the Data tab first."))
      }
      d <- design()
      kind <- input$ordination %||% "pcoa"
      if (identical(kind, "pcoa") && !d$can_pcoa) {
        return(vpf_message_plot(paste(
          "Two-dimensional PCoA requires at least three samples: classical MDS",
          "has at most n - 1 axes. With", d$n_samples,
          "sample(s) only the pairwise distance itself is available, shown on the",
          "Distances tab."
        )))
      }
      if (identical(kind, "nmds") && !d$can_nmds) {
        return(vpf_message_plot(paste0(
          "NMDS is not shown for ", d$n_samples, " samples. Two dimensions are ",
          "saturated at this sample size, so the stress would be zero or near ",
          "zero regardless of the data and would not indicate a good fit. Six ",
          "samples are needed before a two-dimensional NMDS constrains anything."
        )))
      }
      dd <- dissimilarity()
      if (is.null(dd)) {
        return(vpf_message_plot("At least two samples with non-zero signal are needed."))
      }
      if (!is.null(dd$error)) {
        return(vpf_message_plot(dd$error))
      }

      res <- withProgress(message = "Computing ordination", value = 0.5, {
        tryCatch({
          if (identical(kind, "nmds")) {
            fit <- vegan::metaMDS(dd$dist, k = 2, trymax = 50, trace = FALSE,
                                  autotransform = FALSE)
            list(points = vegan::scores(fit, display = "sites"),
                 xlab = "NMDS 1", ylab = "NMDS 2",
                 subtitle = sprintf("stress = %.4f, converged = %s",
                                    fit$stress, isTRUE(fit$converged >= 1)),
                 stress = fit$stress, converged = fit$converged)
          } else {
            # cmdscale warns when fewer than k eigenvalues are positive. That
            # condition is detected and reported below in plainer terms, so the
            # warning would only duplicate it in the console.
            fit <- suppressWarnings(stats::cmdscale(dd$dist, k = 2, eig = TRUE))
            eig <- fit$eig
            pos <- eig[eig > 0]
            if (ncol(as.matrix(fit$points)) < 2) {
              return(list(fail = paste(
                "This dissimilarity matrix has fewer than two positive PCoA axes,",
                "so a two-dimensional classical MDS plot is not available."
              )))
            }
            rel <- 100 * eig / sum(pos)
            list(points = fit$points,
                 xlab = sprintf("Axis 1 (%.1f%%)", rel[[1]]),
                 ylab = sprintf("Axis 2 (%.1f%%)", rel[[2]]),
                 subtitle = if (any(eig < -1e-8)) {
                   "Negative eigenvalues present: this dissimilarity is not Euclidean."
                 } else {
                   ""
                 })
          }
        }, error = function(e) list(fail = conditionMessage(e)))
      })
      if (!is.null(res$fail)) {
        return(vpf_message_plot(res$fail))
      }

      pts <- as.data.frame(res$points)
      colnames(pts)[1:2] <- c("Axis1", "Axis2")
      pts$Sample <- dd$samples
      if (!is.null(dd$groups)) {
        pts$Group <- dd$groups
        pts$Group[is.na(pts$Group)] <- "not annotated"
      }
      p <- ggplot2::ggplot(pts, ggplot2::aes(
        x = .data$Axis1, y = .data$Axis2,
        text = paste0(.data$Sample, "<br>", round(.data$Axis1, 4), ", ",
                      round(.data$Axis2, 4))
      ))
      if ("Group" %in% colnames(pts)) {
        p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data$Group), size = 4)
      } else {
        p <- p + ggplot2::geom_point(size = 4, colour = "#4E79A7")
      }
      # plotly cannot convert ggrepel layers, so labels are placed directly.
      p <- p +
        ggplot2::geom_text(ggplot2::aes(label = .data$Sample), size = 3,
                           vjust = -1.1, show.legend = FALSE) +
        ggplot2::labs(x = res$xlab, y = res$ylab, title = res$subtitle) +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    output$ordination_note <- renderUI({
      d <- design()
      kind <- input$ordination %||% "pcoa"
      notes <- list()
      if (identical(kind, "pcoa") && d$pcoa_saturated) {
        notes <- c(notes, list(vpf_notice(
          type = "warning",
          paste("With three samples a two-dimensional PCoA is saturated: the plot",
                "simply re-expresses the three pairwise distances and reveals no",
                "lower-dimensional structure.")
        )))
      }
      if (d$can_pcoa && d$n_samples <= 5) {
        notes <- c(notes, list(vpf_caption(
          "Confidence ellipses and clustering claims are deliberately omitted at",
          "this sample size."
        )))
      }
      if (length(notes) == 0) NULL else do.call(tagList, notes)
    })

    output$plt_distance <- plotly::renderPlotly({
      dd <- dissimilarity()
      if (is.null(dd)) {
        return(vpf_message_plot(paste(
          "At least two samples with non-zero signal are needed for a",
          "dissimilarity."
        )))
      }
      if (!is.null(dd$error)) {
        return(vpf_message_plot(dd$error))
      }
      m <- as.matrix(dd$dist)
      df <- data.frame(
        A = rep(rownames(m), times = ncol(m)),
        B = rep(colnames(m), each = nrow(m)),
        Value = as.numeric(m),
        stringsAsFactors = FALSE
      )
      df$A <- factor(df$A, levels = rownames(m))
      df$B <- factor(df$B, levels = rev(rownames(m)))
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$A, y = .data$B, fill = .data$Value,
        text = paste0(.data$A, " vs ", .data$B, "<br>",
                      dd$method, " = ", round(.data$Value, 4))
      )) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c",
                                     limits = c(0, 1)) +
        ggplot2::labs(x = NULL, y = NULL, fill = dd$method) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      vpf_ggplotly(p, tooltip = "text")
    })

    output$permanova_out <- renderUI({
      d <- design()
      if (d$n_samples == 0) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      if (d$n_groups < 2) {
        return(vpf_notice(
          type = "warning",
          title = "PERMANOVA is not defined here",
          "It requires a grouping variable with at least two levels.",
          "Sample names are identifiers, not groups."
        ))
      }
      if (d$min_group_size < 2) {
        return(vpf_notice(
          type = "warning",
          title = "PERMANOVA is not defined here",
          "At least two biological samples per group are required; a singleton",
          "group provides no within-group replication."
        ))
      }
      if (!d$can_permanova) {
        return(vpf_notice(
          type = "warning",
          title = "PERMANOVA is not defined here",
          sprintf("At least four samples are required; this selection has %d.",
                  d$n_samples)
        ))
      }
      dd <- dissimilarity()
      if (is.null(dd) || !is.null(dd$error)) {
        return(vpf_notice(type = "danger",
                          dd$error %||% "No usable dissimilarity matrix."))
      }
      perms <- suppressWarnings(as.integer(input$permutations))
      if (is.na(perms) || perms < 99) perms <- 999L

      grp <- dd$groups
      if (is.null(grp) || length(unique(stats::na.omit(grp))) < 2) {
        return(vpf_notice(type = "warning",
                          "The grouping variable has fewer than two levels among the samples retained."))
      }
      meta <- data.frame(Group = factor(grp), stringsAsFactors = FALSE)
      res <- withProgress(message = "Running PERMANOVA", value = 0.5, {
        tryCatch(
          vegan::adonis2(dd$dist ~ Group, data = meta, permutations = perms),
          error = function(e) e
        )
      })
      if (inherits(res, "error")) {
        return(vpf_notice(type = "danger",
                          paste("PERMANOVA failed:", conditionMessage(res))))
      }
      tab <- as.data.frame(res)
      tab <- cbind(Term = rownames(tab), tab)
      p_col <- grep("^Pr", colnames(tab), value = TRUE)
      p_val <- if (length(p_col) > 0) tab[[p_col[[1]]]][[1]] else NA_real_

      resolution <- if (!d$permanova_resolution_ok) {
        vpf_notice(
          type = "warning",
          title = "This p-value cannot fall below 0.05",
          sprintf(paste("Only %.0f distinguishable group-label permutations exist for",
                        "group sizes %s, so the smallest attainable p-value is %.3f.",
                        "Interpret the pseudo-F and the distances descriptively."),
                  d$n_permutations,
                  paste(as.integer(d$group_sizes), collapse = " vs "),
                  1 / d$n_permutations)
        )
      } else {
        NULL
      }

      tagList(
        vpf_notice(
          type = "info",
          title = sprintf("PERMANOVA on %s dissimilarity, %d permutations",
                          dd$method, perms),
          sprintf("Grouping variable: %s. p = %s.", d$group,
                  if (is.na(p_val)) "not returned" else format(p_val, digits = 4))
        ),
        resolution,
        vpf_notice(
          type = "info",
          paste("PERMANOVA confounds differences in group centroids with",
                "differences in group dispersion. At these group sizes a dispersion",
                "test is itself poorly supported, so a significant result should not",
                "be read as evidence of a location shift specifically.")
        ),
        reactable::reactable(vpf_round_df(tab), striped = TRUE, bordered = TRUE,
                             pagination = FALSE)
      )
    })
  })
}

## To be copied in the UI
# mod_diversity_ui("diversity_1")

## To be copied in the server
# mod_diversity_server("diversity_1", r_filter)
