#' genes UI Function
#'
#' @description Gene-level annotations from `metadata(tse)$gene_annotations`:
#'   DRAM-v auxiliary metabolic genes, pharokka PHROG functional categories and
#'   CARD/VFDB hits. The table is restricted to the contigs that pass the
#'   current filters, so the counts always match the rest of the app.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_genes_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("status")),
    tabsetPanel(
      id = ns("tabs"),
      tabPanel(
        "Functional categories",
        fluidRow(
          column(6, plotly::plotlyOutput(ns("plt_phrog"), height = "400px")),
          column(6, plotly::plotlyOutput(ns("plt_amg"), height = "400px"))
        ),
        fluidRow(
          column(6, plotly::plotlyOutput(ns("plt_card"), height = "400px")),
          column(6, plotly::plotlyOutput(ns("plt_vfdb"), height = "400px"))
        ),
        uiOutput(ns("evidence_caveat"))
      ),
      tabPanel(
        "CheckAMG",
        fluidRow(
          column(6, plotly::plotlyOutput(ns("plt_checkamg_class"), height = "400px")),
          column(6, plotly::plotlyOutput(ns("plt_checkamg_function"), height = "400px"))
        ),
        vpf_caption(
          "CheckAMG separates auxiliary genes into metabolic (AMG), physiological",
          "(APG) and regulatory (AReG) classes and reports how confident it is",
          "that each protein is of viral origin. Genes outside a strict viral",
          "region deserve the most scepticism."
        ),
        uiOutput(ns("checkamg_status"))
      ),
      tabPanel(
        "Gene table",
        vpf_caption("Every annotated gene on the currently selected contigs."),
        mod_bigtable_ui(ns("tbl"), height = "520px")
      ),
      tabPanel(
        "Per-contig gene counts",
        plotly::plotlyOutput(ns("plt_per_contig"), height = "420px"),
        vpf_caption(
          "Number of annotated genes per contig, which scales with contig length",
          "and completeness rather than with biological gene density."
        )
      )
    )
  )
}

#' genes Server Functions
#'
#' @param id Module id.
#' @param r_filter A list of reactives from [mod_vpfilter_server()].
#'
#' @noRd
mod_genes_server <- function(id, r_filter) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    gene_table <- reactive({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(NULL)
      }
      ga <- vpf_gene_annotations(r_filter$raw())
      if (is.null(ga)) {
        return(NULL)
      }
      key <- vpf_gene_contig_column(ga)
      if (is.null(key)) {
        return(ga)
      }
      ga[as.character(ga[[key]]) %in% rownames(tse), , drop = FALSE]
    })

    output$status <- renderUI({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      if (is.null(vpf_gene_annotations(r_filter$raw()))) {
        return(vpf_missing_annotation(
          "Gene-level annotations",
          "--use_dram (DRAM-v) and pharokka",
          paste("This panel reads metadata(tse)$gene_annotations, a per-gene table",
                "holding DRAM-v categories and KO/Pfam hits, pharokka PHROG",
                "functional categories, and CARD and VFDB matches. The loaded",
                "object has no such table, which means the run did not produce",
                "gene-level annotations or they were not merged into the result",
                "object by create_vpftse().")
        ))
      }
      ga <- gene_table()
      if (is.null(ga) || nrow(ga) == 0) {
        return(vpf_notice(
          type = "warning",
          "No annotated gene belongs to a contig that passes the current filters."
        ))
      }
      vpf_notice(
        type = "success",
        {
          key <- vpf_gene_contig_column(ga)
          n_contigs <- if (is.null(key)) NA_integer_ else length(unique(ga[[key]]))
          sprintf("%s annotated genes on %s of the selected contigs.",
                  format(nrow(ga), big.mark = ","),
                  if (is.na(n_contigs)) "an unknown number" else
                    format(n_contigs, big.mark = ","))
        }
      )
    })

    #' Bar chart of the most frequent values of one annotation column
    #' @noRd
    category_plot <- function(column, title, colour, absent_message,
                              filter_column = NULL, filter_value = NULL,
                              top_n = 20) {
      ga <- gene_table()
      if (is.null(ga)) {
        return(vpf_message_plot(paste(
          "No gene-level annotations in this dataset. They are produced by the",
          "ViroProfiler options --use_dram and pharokka."
        )))
      }
      if (!column %in% colnames(ga)) {
        return(vpf_message_plot(absent_message))
      }
      df <- ga
      if (!is.null(filter_column) && filter_column %in% colnames(df)) {
        keep <- !is.na(vpf_blank_na(df[[filter_column]])) &
          vpf_blank_na(df[[filter_column]]) == filter_value
        df <- df[keep, , drop = FALSE]
      }
      vals <- vpf_blank_na(df[[column]])
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) {
        return(vpf_message_plot(absent_message))
      }
      counts <- as.data.frame(table(Label = vals), stringsAsFactors = FALSE)
      counts <- counts[order(-counts$Freq), , drop = FALSE]
      counts <- utils::head(counts, top_n)
      counts$Label <- factor(counts$Label, levels = rev(counts$Label))
      p <- ggplot2::ggplot(counts, ggplot2::aes(
        x = .data$Label, y = .data$Freq,
        text = paste0(.data$Label, ": ", .data$Freq, " genes")
      )) +
        ggplot2::geom_col(fill = colour) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Genes", title = title) +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    }

    # Which gene annotator produced this table is not fixed: a run may carry
    # pharokka columns, DRAM-v columns, CheckAMG columns, or a merge of
    # several. Each panel therefore takes the first column that is present
    # rather than assuming one naming scheme.
    first_column <- function(...) {
      ga <- gene_table()
      candidates <- c(...)
      if (is.null(ga)) {
        return(candidates[[1]])
      }
      hit <- intersect(candidates, colnames(ga))
      if (length(hit) == 0) candidates[[1]] else hit[[1]]
    }

    output$plt_phrog <- plotly::renderPlotly({
      column <- first_column("pharokka_category", "checkamg_phrog", "phrog_category")
      title <- if (identical(column, "checkamg_phrog")) {
        "PHROG groups (CheckAMG)"
      } else {
        "PHROG functional categories (pharokka)"
      }
      category_plot(
        column, title, "#4E79A7",
        paste("No PHROG functional category in the gene table. It is produced by",
              "the pharokka annotation step, and CheckAMG also reports PHROG",
              "numbers when --use_checkamg is enabled.")
      )
    })

    output$plt_amg <- plotly::renderPlotly({
      ga <- gene_table()
      if (!is.null(ga) && "dramv_category" %in% colnames(ga)) {
        return(category_plot(
          "dramv_ko", "DRAM-v auxiliary metabolic genes (KO)", "#E15759",
          paste("No DRAM-v AMG hit among the selected contigs. AMG calls come",
                "from the ViroProfiler option --use_dram."),
          filter_column = "dramv_category", filter_value = "AMG"
        ))
      }
      column <- first_column("dramv_ko", "checkamg_kegg_ko")
      title <- if (identical(column, "checkamg_kegg_ko")) {
        "KEGG orthologs (CheckAMG)"
      } else {
        "DRAM-v KEGG orthologs"
      }
      category_plot(
        column, title, "#E15759",
        paste("No KEGG ortholog assignment in the gene table. It comes from the",
              "ViroProfiler options --use_dram (DRAM-v) or --use_checkamg.")
      )
    })

    output$plt_card <- plotly::renderPlotly({
      category_plot(
        "pharokka_card", "CARD antibiotic-resistance hits", "#EDC948",
        paste("No CARD hit among the selected contigs. CARD matches are reported",
              "by pharokka, and by the ViroProfiler option --use_abricate.")
      )
    })

    output$plt_vfdb <- plotly::renderPlotly({
      category_plot(
        "pharokka_vfdb", "VFDB virulence-factor hits", "#B07AA1",
        paste("No VFDB hit among the selected contigs. VFDB matches are reported",
              "by pharokka, and by the ViroProfiler option --use_abricate.")
      )
    })

    output$plt_checkamg_class <- plotly::renderPlotly({
      category_plot(
        "checkamg_class", "CheckAMG auxiliary-gene classes", "#4E79A7",
        paste("No CheckAMG annotation in the gene table. It is produced by the",
              "ViroProfiler option --use_checkamg.")
      )
    })

    output$plt_checkamg_function <- plotly::renderPlotly({
      column <- first_column("checkamg_function", "checkamg_kegg_ko",
                             "checkamg_pfam", "checkamg_cazy")
      category_plot(
        column, paste0("CheckAMG annotations (", column, ")"), "#F28E2B",
        paste("No CheckAMG function annotation in the gene table. It is produced",
              "by the ViroProfiler option --use_checkamg.")
      )
    })

    output$checkamg_status <- renderUI({
      ga <- gene_table()
      if (is.null(ga) || !"checkamg_class" %in% colnames(ga)) {
        return(NULL)
      }
      if (!"checkamg_in_viral_region" %in% colnames(ga)) {
        return(NULL)
      }
      flag <- vpf_blank_na(ga$checkamg_in_viral_region)
      outside <- sum(!is.na(flag) & tolower(flag) %in% c("false", "no", "0"))
      if (outside == 0) {
        return(NULL)
      }
      vpf_notice(
        type = "warning",
        sprintf(paste("%d of %d annotated proteins lie outside a strict viral",
                      "region. An auxiliary-gene call there may belong to residual",
                      "host sequence rather than to the virus."),
                outside, nrow(ga))
      )
    })

    output$evidence_caveat <- renderUI({
      if (is.null(gene_table())) {
        return(NULL)
      }
      vpf_notice(
        type = "warning",
        title = "How to read these hits",
        tags$ul(
          tags$li(paste("A CARD or VFDB match is a sequence similarity hit, not a",
                        "demonstrated resistance or virulence phenotype.")),
          tags$li(paste("A hit on a contig with residual host sequence may belong",
                        "to the host rather than to the virus. Check the contig's",
                        "CheckV contamination on the Contigs tab before drawing a",
                        "conclusion.")),
          tags$li(paste("DRAM-v labels a gene as auxiliary from its genomic context;",
                        "the call is far weaker on a fragmentary contig than on a",
                        "complete genome."))
        )
      )
    })

    output$plt_per_contig <- plotly::renderPlotly({
      ga <- gene_table()
      if (is.null(ga)) {
        return(vpf_message_plot(paste(
          "No gene-level annotations in this dataset. They are produced by the",
          "ViroProfiler options --use_dram and pharokka."
        )))
      }
      key <- vpf_gene_contig_column(ga)
      if (is.null(key)) {
        return(vpf_message_plot("The gene table has no contig identifier column."))
      }
      counts <- as.data.frame(table(Contig = as.character(ga[[key]])),
                              stringsAsFactors = FALSE)
      if (nrow(counts) == 0) {
        return(vpf_message_plot("No genes on the currently selected contigs."))
      }
      counts <- counts[order(-counts$Freq), , drop = FALSE]
      counts <- utils::head(counts, 40)
      counts$Contig <- factor(counts$Contig, levels = rev(counts$Contig))
      p <- ggplot2::ggplot(counts, ggplot2::aes(
        x = .data$Contig, y = .data$Freq,
        text = paste0(.data$Contig, ": ", .data$Freq, " annotated genes")
      )) +
        ggplot2::geom_col(fill = "#59A14F") +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Annotated genes") +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    mod_bigtable_server("tbl", r_data = gene_table, filename = "gene_annotations")
  })
}

## To be copied in the UI
# mod_genes_ui("genes_1")

## To be copied in the server
# mod_genes_server("genes_1", r_filter)
