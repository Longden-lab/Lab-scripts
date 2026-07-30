# app.R — Mural cell gene expression, side-by-side dataset comparison
#
# Left column  = our data (hippocampus)
# Right column = published data (whole brain)
# Rows are matched by plot type so the two can be read one against the other.
#
# Secrets (set in Connect Cloud, never committed):
#   APP_PASSWORD     — password users type on the login screen
#   DATA_REPO_TOKEN  — GitHub PAT with read access to the private data repo
#
# Neither .rds is part of the deployed bundle. Both are pulled from the private
# repo at startup into tempdir(), read into memory, then deleted from disk.

library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(httr)

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit this block, leave the rest alone
# ══════════════════════════════════════════════════════════════════════════════

data_repo <- "Qbottle/data"

# Shared label vocabulary. Both datasets are mapped onto these levels, so a
# label and a colour mean the same cell type across the whole page.
panel_levels <- c("aSMC", "C_PC", "Ts_PC", "vSMC")
panel_cols   <- c(aSMC  = "firebrick",
                  C_PC  = "darkorange2",
                  Ts_PC = "goldenrod3",
                  vSMC  = "steelblue")

# ggplot accepts R colour names; CSS does not. "darkorange2" and "goldenrod3"
# are R-only, so HTML swatches using them render blank. Convert once here and
# use the hex form for every HTML element.
panel_cols_hex <- apply(grDevices::col2rgb(panel_cols), 2,
                        function(x) sprintf("#%02X%02X%02X", x[1], x[2], x[3]))
names(panel_cols_hex) <- names(panel_cols)

# Both objects carry `mural_final` with the same five levels, so a single
# shared map serves both and the two configs differ only in provenance.
shared_label_map <- c(aSMC = "aSMC", aaSMC = "aSMC",
                      C_PC = "C_PC", Ts_PC = "Ts_PC", vSMC = "vSMC")

# Column tints. Cool blue-grey on the left echoes the vSMC steelblue; warm sand
# on the right echoes the goldenrod/orange end. Both are desaturated enough that
# plot ink stays dominant. Change these four pairs to restyle the whole page.
tint <- list(
  left  = list(bg = "#eef2f7", edge = "#b9cbdf", head = "#dae5f1"),
  right = list(bg = "#faf5ec", edge = "#dfcdaf", head = "#f2e7d3")
)

# ── Left panel: our data ──
cfg_left <- list(
  key        = "left",
  name       = "Our data",
  tissue     = "Hippocampus",
  note       = NULL,          # small text beside the name; NULL to omit
  platform   = NULL,          # EDIT: e.g. "10x, 3' UMI" — NULL to omit
  file       = "data/mural_obj_app.rds",
  source_col = "mural_final",
  label_map  = shared_label_map
)

# ── Right panel: published data ──
cfg_right <- list(
  key        = "right",
  name       = "Published (Betsholtz)",
  tissue     = "Whole brain",
  note       = "re-annotated by Gyu",
  platform   = "Plate-based, full-length reads",
  file       = "data/Betsholtz_mural_app.rds",
  source_col = "mural_final",
  label_map  = shared_label_map
)

max_genes <- 8L   # 8 genes x 2 panels

# Facet columns, tuned for half-width columns rather than full width
ncol_bar     <- 3L
ncol_violin  <- 2L
ncol_feature <- 2L

# ══════════════════════════════════════════════════════════════════════════════
# Secrets — fail loudly at startup rather than silently at login
# ══════════════════════════════════════════════════════════════════════════════

app_password <- Sys.getenv("APP_PASSWORD")
if (!nzchar(app_password)) {
  stop("APP_PASSWORD is not set. Add it as a secret in Connect Cloud.")
}

data_token <- Sys.getenv("DATA_REPO_TOKEN")
if (!nzchar(data_token)) {
  stop("DATA_REPO_TOKEN is not set. Add it as a secret in Connect Cloud.")
}

# ══════════════════════════════════════════════════════════════════════════════
# Fetch and load
# ══════════════════════════════════════════════════════════════════════════════

fetch_from_repo <- function(repo_path, dest) {
  api_url <- sprintf("https://api.github.com/repos/%s/contents/%s", data_repo, repo_path)
  partial <- paste0(dest, ".part")
  on.exit(unlink(partial), add = TRUE)

  resp <- httr::GET(
    api_url,
    httr::add_headers(
      Authorization          = paste("Bearer", data_token),
      Accept                 = "application/vnd.github.raw",
      `X-GitHub-Api-Version` = "2022-11-28"
    ),
    httr::write_disk(partial, overwrite = TRUE),
    httr::timeout(900)
  )
  httr::stop_for_status(resp)

  # A failed request can still leave a small JSON error body on disk, and an
  # interrupted transfer leaves a truncated file. Either would load as if valid.
  if (!file.exists(partial) || file.size(partial) < 1e5) {
    stop("Download of ", repo_path, " failed or was truncated (",
         if (file.exists(partial)) file.size(partial) else 0, " bytes).")
  }
  if (!file.rename(partial, dest)) {
    stop("Could not move ", repo_path, " into place.")
  }
  invisible(dest)
}

# Plot backgrounds must be see-through for the column tint to read as one band
theme_clear <- theme(
  plot.background       = element_rect(fill = "transparent", colour = NA),
  panel.background      = element_rect(fill = "transparent", colour = NA),
  legend.background     = element_rect(fill = "transparent", colour = NA),
  legend.box.background = element_rect(fill = "transparent", colour = NA),
  legend.key            = element_rect(fill = "transparent", colour = NA),
  strip.background      = element_rect(fill = "transparent", colour = NA)
)

# Turns a config entry into everything the plot builders and UI need.
load_dataset <- function(cfg) {
  message("Loading ", cfg$name, " ...")
  dest <- file.path(tempdir(), paste0(cfg$key, ".rds"))
  fetch_from_repo(cfg$file, dest)
  obj <- readRDS(dest)
  unlink(dest)   # private data does not linger on disk

  DefaultAssay(obj) <- "RNA"
  if (!"umap" %in% Reductions(obj)) {
    stop(cfg$name, ": no 'umap' reduction. Add one in prepare_objects.R.")
  }
  if (!cfg$source_col %in% colnames(obj@meta.data)) {
    stop(cfg$name, ": metadata column '", cfg$source_col, "' not found. Found: ",
         paste(colnames(obj@meta.data), collapse = ", "))
  }

  # Map source labels onto the shared vocabulary, drop anything unmapped
  mapped   <- unname(cfg$label_map[as.character(obj@meta.data[[cfg$source_col]])])
  unmapped <- setdiff(as.character(obj@meta.data[[cfg$source_col]]), names(cfg$label_map))
  if (length(unmapped)) {
    message(cfg$name, ": dropping unmapped labels: ", paste(unmapped, collapse = ", "))
  }
  obj$panel_class <- factor(mapped, levels = panel_levels)
  obj <- obj[, !is.na(obj$panel_class)]

  present <- panel_levels[panel_levels %in% levels(droplevels(obj$panel_class))]
  genes   <- sort(rownames(obj))

  # Case-insensitive index. Both datasets use mouse symbols with identical
  # capitalisation, so no ortholog aliasing is needed; this only forgives typing
  # "kcnj8" for "Kcnj8". Add an alias table here if a cross-species set is added.
  lookup <- setNames(genes, toupper(genes))

  # Genes that exist as rownames but are zero in every cell. These would render
  # a blank panel with no explanation, so they are reported separately from
  # genes that are genuinely absent.
  expressed <- genes[Matrix::rowSums(LayerData(obj, assay = "RNA", layer = "data")) > 0]

  # Counts per cell type, over the full shared vocabulary so the left and right
  # tables always have the same number of rows and the columns stay aligned.
  counts <- table(factor(as.character(obj$panel_class), levels = panel_levels))

  # No on-plot labels and no legend: the count table sits where the legend was
  # and doubles as the colour key.
  ref_plot <- DimPlot(obj, group.by = "panel_class", reduction = "umap",
                      label = FALSE, pt.size = 0.5,
                      cols = panel_cols[present]) +
    labs(title = cfg$name, subtitle = cfg$tissue) +
    theme(plot.title      = element_text(face = "bold", size = 12),
          plot.subtitle   = element_text(size = 10, color = "grey30"),
          legend.position = "none") +
    theme_clear

  list(
    key       = cfg$key,
    name      = cfg$name,
    tissue    = cfg$tissue,
    note      = cfg$note,
    platform  = cfg$platform,
    obj       = obj,
    present   = present,
    genes     = genes,
    lookup    = lookup,
    expressed = expressed,
    counts    = counts,
    n_cells   = ncol(obj),
    ref_plot  = ref_plot
  )
}

datasets <- list(left = load_dataset(cfg_left), right = load_dataset(cfg_right))

message("Loaded both datasets. Total in-memory size: ",
        format(object.size(datasets), units = "MB"))

# ══════════════════════════════════════════════════════════════════════════════
# Gene resolution
# ══════════════════════════════════════════════════════════════════════════════

# One query symbol can be absent from one dataset, so resolve and report per
# dataset rather than assuming a shared gene space.
resolve_genes <- function(txt) {
  if (is.null(txt)) txt <- ""
  raw <- trimws(unlist(strsplit(txt, "[,;[:space:]]+")))
  raw <- unique(raw[nzchar(raw)])

  in_any <- vapply(raw, function(q) {
    any(vapply(datasets, function(d) !is.na(d$lookup[toupper(q)]), logical(1)))
  }, logical(1))

  queries <- head(raw[in_any], max_genes)

  per_dataset <- lapply(datasets, function(d) {
    hit   <- d$lookup[toupper(queries)]
    found <- unname(hit[!is.na(hit)])
    list(found      = found,
         missing    = queries[is.na(hit)],           # not in this dataset
         undetected = setdiff(found, d$expressed))   # present but all-zero
  })

  list(
    queries     = queries,
    unknown     = raw[!in_any],
    dropped     = max(0L, sum(in_any) - max_genes),
    per_dataset = per_dataset,
    # Row heights come from the larger side so the two columns stay aligned
    n_rows_for  = max(vapply(per_dataset, function(p) length(p$found), integer(1)))
  )
}

grid_height <- function(n, ncol_max, row_h) {
  if (n < 1) return(row_h)
  row_h * ceiling(n / min(ncol_max, n))
}

# ══════════════════════════════════════════════════════════════════════════════
# Plot builders — every one takes a dataset as its first argument
# ══════════════════════════════════════════════════════════════════════════════

make_bar <- function(ds, genes) {
  df <- FetchData(ds$obj, vars = c(genes, "panel_class"), layer = "data") %>%
    pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") %>%
    group_by(gene, panel_class) %>%
    summarise(
      mean = mean(expr),
      sem  = if (n() > 1) sd(expr) / sqrt(n()) else 0,   # sd() is NA for n = 1
      .groups = "drop"
    ) %>%
    mutate(gene = factor(gene, levels = genes))

  ggplot(df, aes(panel_class, mean, fill = panel_class)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = mean, ymax = mean + sem), width = 0.2, color = "grey30") +
    geom_text(aes(label = round(mean, 2), y = mean + sem), vjust = -0.4,
              size = 2.6, color = "grey20") +
    facet_wrap(~ gene, ncol = min(ncol_bar, length(genes)), scales = "free_y") +
    scale_fill_manual(values = panel_cols, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = ds$name, x = NULL, y = "Mean expression (log-norm)") +
    theme_classic(base_size = 11) +
    theme(plot.title  = element_text(face = "bold", size = 12),
          strip.text  = element_text(face = "bold.italic"),
          legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    theme_clear
}

make_violin <- function(ds, genes) {
  VlnPlot(ds$obj, features = genes, group.by = "panel_class", pt.size = 0,
          cols = panel_cols[ds$present],
          ncol = min(ncol_violin, length(genes))) &
    theme(plot.title   = element_text(size = 11, face = "bold.italic"),
          axis.title.x = element_blank(),
          axis.text.x  = element_text(angle = 30, hjust = 1)) &
    theme_clear
}

make_feature <- function(ds, genes) {
  FeaturePlot(ds$obj, features = genes, reduction = "umap", order = TRUE,
              pt.size = 0.5, cols = c("lightgrey", "firebrick"),
              ncol = min(ncol_feature, length(genes))) &
    theme(plot.title = element_text(size = 11, face = "bold.italic"),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank()) &
    theme_clear
}

# ══════════════════════════════════════════════════════════════════════════════
# UI pieces
# ══════════════════════════════════════════════════════════════════════════════

# Cell counts per type, rendered beside the UMAP that shows the distribution.
# Static, so it is built once and embedded directly rather than via an output.
count_table <- function(ds) {
  n     <- as.integer(ds$counts)
  names(n) <- names(ds$counts)
  total <- sum(n)

  body <- lapply(panel_levels, function(k) {
    tags$tr(
      tags$td(
        style = "padding:2px 6px 2px 0;",
        tags$span(style = sprintf(
          "display:inline-block; width:9px; height:9px; border-radius:2px;
           background:%s; margin-right:5px;", panel_cols_hex[[k]])),
        k
      ),
      tags$td(style = "padding:2px 7px; text-align:right;
                       font-variant-numeric:tabular-nums;",
              format(n[[k]], big.mark = ",")),
      tags$td(style = "padding:2px 0; text-align:right; color:#777;
                       font-variant-numeric:tabular-nums;",
              if (total > 0) sprintf("%.1f%%", 100 * n[[k]] / total) else "—")
    )
  })

  tags$table(
    style = "font-size:11.5px; margin:0 auto; border-collapse:collapse;
             white-space:nowrap;",
    tags$thead(tags$tr(
      tags$th(style = "text-align:left; padding-bottom:3px; border-bottom:1px solid #ccc;",
              "Cell type"),
      tags$th(style = "text-align:right; padding:0 7px 3px 7px;
                       border-bottom:1px solid #ccc;",
              "Cells"),
      tags$th(style = "text-align:right; padding-bottom:3px; border-bottom:1px solid #ccc;",
              "%")
    )),
    tags$tbody(
      body,
      tags$tr(
        tags$td(style = "padding-top:4px; border-top:1px solid #ccc; font-weight:600;",
                "Total"),
        tags$td(style = "padding:4px 7px 0 7px; text-align:right; font-weight:600;
                         border-top:1px solid #ccc; font-variant-numeric:tabular-nums;",
                format(total, big.mark = ",")),
        tags$td(style = "border-top:1px solid #ccc;")
      )
    )
  )
}

# One comparison row: a full-width heading strip, then the two tinted columns.
# `extra` is a named list of additional content to place under each plot.
compare_row <- function(heading, note, output_prefix, ref_height = NULL,
                        extra = NULL, extra_at = c("below", "right")) {
  extra_at <- match.arg(extra_at)

  side_col <- function(side) {
    id <- paste0(output_prefix, "_", side)

    # height="auto" matters: the default plotOutput container is 400px, and any
    # renderPlot image taller than that overflows the tinted band.
    plt <- if (is.null(ref_height)) {
      plotOutput(id, height = "auto")
    } else {
      plotOutput(id, height = ref_height)
    }

    inner <- if (!is.null(extra[[side]]) && extra_at == "right") {
      fluidRow(
        column(7, plt),
        column(5, div(
          style = sprintf("display:flex; align-items:center; justify-content:center;
                           height:%s;", if (is.null(ref_height)) "auto" else ref_height),
          extra[[side]]
        ))
      )
    } else {
      tagList(plt, if (!is.null(extra[[side]])) extra[[side]])
    }

    column(6, div(class = paste0("cmp-band cmp-", side), inner))
  }
  tagList(
    div(class = "cmp-head",
        span(class = "cmp-head-title", heading),
        if (!is.null(note)) span(class = "cmp-head-note", note)),
    fluidRow(side_col("left"), side_col("right"))
  )
}

# Column header: name, small provenance note beside it, tissue underneath
column_header <- function(ds, side) {
  div(
    class = paste0("cmp-colhead cmp-colhead-", side),
    div(
      span(class = "cmp-colhead-name", ds$name),
      if (!is.null(ds$note)) span(class = "cmp-colhead-note", ds$note)
    ),
    div(class = "cmp-colhead-tissue", ds$tissue)
  )
}

app_css <- sprintf("
  .cmp-head {
    background:#3f4756; border-radius:3px 3px 0 0;
    padding:7px 12px; margin:14px 15px 0 15px;
  }
  .cmp-head-title { font-weight:600; font-size:15px; color:#ffffff;
                    letter-spacing:0.2px; }
  .cmp-head-note  { color:#c3cad6; font-size:12px; margin-left:10px; }

  /* Consecutive bands touch, so each column reads as one continuous strip */
  .cmp-band { padding:10px 12px 12px 12px; }
  .cmp-left  { background:%s; border-left:3px solid %s; }
  .cmp-right { background:%s; border-left:3px solid %s; }

  .cmp-colhead { text-align:center; padding:6px 4px; border-radius:3px 3px 0 0; }
  .cmp-colhead-left  { background:%s; border-bottom:2px solid %s; }
  .cmp-colhead-right { background:%s; border-bottom:2px solid %s; }
  .cmp-colhead-name   { font-weight:700; font-size:15px; }
  .cmp-colhead-note   { font-size:11.5px; color:#5c5c5c; font-style:italic;
                        margin-left:8px; }
  .cmp-colhead-tissue { font-size:12px; color:#4a4a4a; margin-top:1px; }

  .cmp-legend-swatch {
    display:inline-block; width:10px; height:10px; border-radius:2px;
    margin-right:6px; vertical-align:middle;
  }
  .cmp-notebox {
    background:#f1f3f5; border-left:3px solid #8a94a6;
    padding:8px 12px; margin:0 15px 12px 15px; font-size:13px;
  }
",
  tint$left$bg,   tint$left$edge,
  tint$right$bg,  tint$right$edge,
  tint$left$head, tint$left$edge,
  tint$right$head, tint$right$edge
)

main_ui <- function() {
  tagList(
    titlePanel("Mural cell gene expression — dataset comparison"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        textAreaInput("genes", "Gene(s) — comma or space separated",
                      value = "Kcnj8, Adra1a, Ednra", rows = 3,
                      placeholder = "e.g. Kcnj8, Adra1a, Ednra"),
        actionButton("go", "Plot", class = "btn-primary"),
        tags$hr(),
        uiOutput("status"),
        tags$hr(),
        lapply(c("left", "right"), function(side) {
          ds <- datasets[[side]]
          div(
            style = sprintf("background:%s; border-left:3px solid %s;
                             padding:6px 8px; margin-bottom:8px; font-size:12px;",
                            tint[[side]]$bg, tint[[side]]$edge),
            tags$b(ds$name),
            if (!is.null(ds$note)) tags$span(style = "font-style:italic; color:#5c5c5c;",
                                             paste0(" ", ds$note)),
            tags$br(),
            ds$tissue,
            if (!is.null(ds$platform)) tagList(tags$br(), ds$platform),
            tags$br(),
            sprintf("%s genes measured", format(length(ds$genes), big.mark = ","))
          )
        }),
        tags$small(style = "color:#777;",
                   sprintf("Up to %d genes plotted at a time.", max_genes))
      ),
      mainPanel(
        width = 9,

        fluidRow(
          column(6, column_header(datasets$left,  "left")),
          column(6, column_header(datasets$right, "right"))
        ),

        div(class = "cmp-notebox",
            tags$b("Reading these panels."), tags$br(),
            "The two datasets come from different tissue (hippocampus versus whole ",
            "brain), were generated on different platforms, and were normalised ",
            "independently. Compare the ", tags$i("pattern across cell types"),
            " within each panel. Do not compare bar heights or colour intensities ",
            "between the left and right panels: a difference there can reflect ",
            "dissection, library preparation or sequencing depth rather than biology."
        ),

        compare_row(
          "Bar graph — subtypes (PC split: C_PC / Ts_PC)",
          paste("All samples pooled | aSMC = aSMC + aaSMC | y-axes independent |",
                "error bars are SEM, so their width tracks group size as much as spread."),
          "bar"
        ),
        compare_row(
          "Violin — subtypes (PC split: C_PC / Ts_PC)",
          "Distribution behind the means above. Points hidden; width is density.",
          "violin"
        ),
        compare_row(
          "Reference UMAP and cell-type composition",
          "Embeddings are computed independently, so position is not comparable between panels.",
          "ref",
          ref_height = "420px",
          extra = list(left  = count_table(datasets$left),
                       right = count_table(datasets$right)),
          extra_at = "right"
        ),
        compare_row(
          "Feature plot (UMAP)",
          "Colour scales are set per panel by that panel's maximum.",
          "feature"
        )
      )
    )
  )
}

login_ui <- function(msg = NULL) {
  div(style = "max-width:340px; margin:80px auto; text-align:center;",
      h3("Mural cell gene expression"),
      p("Enter the password to continue."),
      passwordInput("pw", NULL, placeholder = "Password"),
      actionButton("login", "Enter", class = "btn-primary"),
      if (!is.null(msg)) tags$p(style = "color:#b00; margin-top:12px;", msg)
  )
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML(app_css)),
    # Enter submits the login form. Lives at the top level so it survives the
    # login -> app swap; the guard means it does nothing once #login is gone.
    tags$script(HTML(
      "document.addEventListener('keydown', function(e) {
         if (e.key !== 'Enter') return;
         var btn = document.getElementById('login');
         if (btn) { e.preventDefault(); btn.click(); }
       });"
    ))
  ),
  uiOutput("page")
)

# ══════════════════════════════════════════════════════════════════════════════
# Server
# ══════════════════════════════════════════════════════════════════════════════

# Wires up the four outputs for one dataset. Called once per side, so the plot
# code exists in exactly one place. bg = "transparent" lets the column tint show
# through instead of a white rectangle.
register_panel <- function(output, ds, genes_r, authed) {
  side   <- ds$key
  found  <- function() genes_r()$per_dataset[[side]]$found
  n_rows <- function() genes_r()$n_rows_for   # shared, so rows stay aligned

  output[[paste0("ref_", side)]] <- renderPlot({
    req(authed())
    ds$ref_plot
  }, bg = "transparent")

  output[[paste0("bar_", side)]] <- renderPlot({
    req(authed())
    g <- found(); req(length(g) > 0)
    make_bar(ds, g)
  }, height = function() grid_height(n_rows(), ncol_bar, 240), bg = "transparent")

  output[[paste0("violin_", side)]] <- renderPlot({
    req(authed())
    g <- found(); req(length(g) > 0)
    make_violin(ds, g)
  }, height = function() grid_height(n_rows(), ncol_violin, 260), bg = "transparent")

  output[[paste0("feature_", side)]] <- renderPlot({
    req(authed())
    g <- found(); req(length(g) > 0)
    make_feature(ds, g)
  }, height = function() grid_height(n_rows(), ncol_feature, 280), bg = "transparent")
}

server <- function(input, output, session) {

  authed    <- reactiveVal(FALSE)
  login_msg <- reactiveVal(NULL)

  # Depends on login_msg(), so a failed attempt re-renders the login screen
  # without detaching this output from authed() — a later correct password
  # works without a page reload.
  output$page <- renderUI({
    if (authed()) main_ui() else login_ui(login_msg())
  })

  observeEvent(input$login, {
    if (identical(input$pw, app_password)) {
      login_msg(NULL)
      authed(TRUE)
    } else {
      login_msg("Incorrect password. Check for stray spaces and try again.")
    }
  })

  genes_r <- eventReactive(input$go, {
    resolve_genes(input$genes)
  }, ignoreNULL = FALSE)

  output$status <- renderUI({
    req(authed())
    g <- genes_r()
    msgs <- list()

    if (length(g$queries)) {
      msgs <- c(msgs, list(tags$p(tags$b("Plotting: "),
                                  paste(g$queries, collapse = ", "))))
    }
    # Present in one dataset but not the other — the asymmetry readers need to see
    for (d in datasets) {
      miss <- g$per_dataset[[d$key]]$missing
      if (length(miss)) {
        msgs <- c(msgs, list(tags$p(
          style = "color:#b00;",
          tags$b(paste0("Not in ", d$name, ": ")), paste(miss, collapse = ", ")
        )))
      }
    }
    for (d in datasets) {
      undet <- g$per_dataset[[d$key]]$undetected
      if (length(undet)) {
        msgs <- c(msgs, list(tags$p(
          style = "color:#b8860b;",
          tags$b(paste0("Zero in every ", d$name, " cell: ")),
          paste(undet, collapse = ", ")
        )))
      }
    }
    if (length(g$unknown)) {
      msgs <- c(msgs, list(tags$p(style = "color:#b00;",
                                  tags$b("In neither dataset: "),
                                  paste(g$unknown, collapse = ", "))))
    }
    if (g$dropped > 0) {
      msgs <- c(msgs, list(tags$p(
        style = "color:#b00;",
        sprintf("Showing the first %d genes; %d more were left out.",
                max_genes, g$dropped))))
    }
    if (!length(g$queries)) {
      msgs <- c(msgs, list(tags$p(style = "color:#b00;",
                                  "Type at least one gene present in either dataset.")))
    }
    tagList(msgs)
  })

  register_panel(output, datasets$left,  genes_r, authed)
  register_panel(output, datasets$right, genes_r, authed)
}

shinyApp(ui, server)
