# app.R — Mural cell gene expression explorer (password-gated)
#
# Secrets (set in Connect Cloud, never committed):
#   APP_PASSWORD     — password users type on the login screen
#   DATA_REPO_TOKEN  — GitHub PAT with read access to the private data repo
#
# The .rds is NEVER part of the deployed bundle. It is pulled from the private
# repo at startup into tempdir(), read into memory, then deleted from disk.

library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(httr)

# ── Secrets: fail loudly at startup rather than silently at login ────────────

app_password <- Sys.getenv("APP_PASSWORD")
if (!nzchar(app_password)) {
  stop("APP_PASSWORD is not set. Add it as a secret in Connect Cloud.")
}

data_token <- Sys.getenv("DATA_REPO_TOKEN")
if (!nzchar(data_token)) {
  stop("DATA_REPO_TOKEN is not set. Add it as a secret in Connect Cloud.")
}

# ── Fetch the data object from the PRIVATE GitHub repo ───────────────────────

data_repo <- "Qbottle/data"                  # private repo (owner/name)
data_file <- "data/mural_obj_for_suyeon.rds" # path inside that repo

# tempdir(), not the app directory: keeps the object out of git and out of any
# path the Shiny server could serve, and lets the OS clean it up.
rds_path <- file.path(tempdir(), "mural_obj.rds")

fetch_data <- function(dest) {
  api_url  <- sprintf("https://api.github.com/repos/%s/contents/%s", data_repo, data_file)
  partial  <- paste0(dest, ".part")
  on.exit(unlink(partial), add = TRUE)

  resp <- httr::GET(
    api_url,
    httr::add_headers(
      Authorization          = paste("Bearer", data_token),
      Accept                 = "application/vnd.github.raw",
      `X-GitHub-Api-Version` = "2022-11-28"
    ),
    httr::write_disk(partial, overwrite = TRUE),
    httr::timeout(600)
  )
  httr::stop_for_status(resp)

  # A failed request can still leave a small JSON error body on disk, and an
  # interrupted transfer leaves a truncated file. Either would be cached as if
  # it were valid, so check before promoting .part to the real filename.
  if (!file.exists(partial) || file.size(partial) < 1e5) {
    stop("Download of the data object failed or was truncated (",
         "got ", if (file.exists(partial)) file.size(partial) else 0, " bytes).")
  }
  if (!file.rename(partial, dest)) {
    stop("Could not move the downloaded data object into place.")
  }
  invisible(dest)
}

fetch_data(rds_path)

# ── Load once at startup, then remove the file from disk ─────────────────────

mu <- readRDS(rds_path)
unlink(rds_path)

DefaultAssay(mu) <- "RNA"
stopifnot("umap" %in% Reductions(mu))

# ── Class definitions (labels live in mural_final) ───────────────────────────

# 3-class grouping (PC collapsed). Not plotted on the page any more, but kept
# so the collapsed view can be restored without redefining anything.
class_map    <- c(aSMC = "aSMC", aaSMC = "aSMC", C_PC = "PC", Ts_PC = "PC", vSMC = "vSMC")
class_levels <- c("aSMC", "PC", "vSMC")
class_cols   <- c(aSMC = "firebrick", PC = "darkorange2", vSMC = "steelblue")

# Subtype grouping — this is what the page shows.

class_map_detail    <- c(aSMC = "aSMC", aaSMC = "aSMC", C_PC = "C_PC",
                         Ts_PC = "Ts_PC", vSMC = "vSMC")
class_levels_detail <- c("aSMC", "C_PC", "Ts_PC", "vSMC")
class_cols_detail   <- c(aSMC = "firebrick", C_PC = "darkorange2",
                         Ts_PC = "goldenrod3", vSMC = "steelblue")

mu$mural_class <- factor(
  unname(class_map[as.character(mu$mural_final)]), levels = class_levels
)
mu$mural_class_detail <- factor(
  unname(class_map_detail[as.character(mu$mural_final)]), levels = class_levels_detail
)

all_genes   <- sort(rownames(mu))
gene_lookup <- setNames(all_genes, toupper(all_genes))  # case-insensitive matching
max_genes   <- 12L                                      # cap plot size per query

# Reference UMAP — identical for every query, so build it once
p_ref <- DimPlot(mu, group.by = "mural_final", reduction = "umap", label = TRUE,
                 repel = TRUE, label.size = 4, pt.size = 0.5) +
  labs(title = "Mural subtypes (reference)") +
  theme(plot.title = element_text(face = "bold", size = 12))

# ── Helpers ──────────────────────────────────────────────────────────────────

resolve_genes <- function(txt) {
  if (is.null(txt)) txt <- ""
  raw <- trimws(unlist(strsplit(txt, "[,;[:space:]]+")))
  raw <- unique(raw[nzchar(raw)])
  hits <- gene_lookup[toupper(raw)]
  found <- unname(hits[!is.na(hits)])
  list(
    found   = head(found, max_genes),
    missing = raw[is.na(hits)],
    dropped = max(0L, length(found) - max_genes)
  )
}

# rows = ceiling(n / ncol), where ncol never exceeds the number of genes
grid_height <- function(n, ncol_max, row_h) {
  if (n < 1) return(row_h)
  row_h * ceiling(n / min(ncol_max, n))
}

# ── Plot builders ────────────────────────────────────────────────────────────

make_bar <- function(genes, class_col, cols, title, subtitle = NULL) {
  df <- FetchData(mu, vars = c(genes, class_col), layer = "data") %>%
    pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") %>%
    group_by(gene, .data[[class_col]]) %>%
    summarise(
      mean = mean(expr),
      # sd() is NA for a single cell; treat that as no error bar
      sem  = if (n() > 1) sd(expr) / sqrt(n()) else 0,
      .groups = "drop"
    ) %>%
    mutate(gene = factor(gene, levels = genes))
  names(df)[2] <- "grp"

  ggplot(df, aes(grp, mean, fill = grp)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = mean, ymax = mean + sem), width = 0.2, color = "grey30") +
    geom_text(aes(label = round(mean, 2), y = mean + sem), vjust = -0.4,
              size = 2.8, color = "grey20") +
    facet_wrap(~ gene, ncol = min(5, length(genes)), scales = "free_y") +
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Mean expression (log-norm)") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          strip.text    = element_text(face = "bold.italic"),
          legend.position = "none",
          axis.text.x   = element_text(angle = 30, hjust = 1))
}

make_feature <- function(genes) {
  FeaturePlot(mu, features = genes, reduction = "umap", order = TRUE,
              pt.size = 0.5, cols = c("lightgrey", "firebrick"),
              ncol = min(4, length(genes))) &
    theme(plot.title = element_text(size = 11, face = "bold.italic"),
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank())
}

make_violin <- function(genes, class_col, cols) {
  VlnPlot(mu, features = genes, group.by = class_col, pt.size = 0,
          cols = cols, ncol = min(4, length(genes))) &
    theme(plot.title = element_text(size = 11, face = "bold.italic"),
          axis.title.x = element_blank(),
          axis.text.x  = element_text(angle = 0, hjust = 0.5))
}

# ── UI: gene explorer, shown only after a correct password ───────────────────

main_ui <- function() {
  tagList(
    titlePanel("Mural cell gene expression"),
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
        tags$small(sprintf("%d genes available in this dataset. Up to %d plotted at a time.",
                           length(all_genes), max_genes))
      ),
      mainPanel(
        width = 9,

        # 1. Mean expression per subtype — the headline numbers
        h4("Bar graph — subtypes (PC split: C_PC / Ts_PC)"),
        plotOutput("bar_detail"),

        # 2. Same grouping, full distribution behind those means
        h4("Violin — subtypes (PC split: C_PC / Ts_PC)"),
        tags$small(style = "color:#666;",
                   "All samples pooled | aSMC = aSMC + aaSMC"),
        plotOutput("violin_detail"),

        # 3. Where the subtypes sit in UMAP space
        h4("Reference UMAP"),
        plotOutput("ref", height = "420px"),

        # 4. Per-cell expression on that same embedding
        h4("Feature plot (UMAP)"),
        plotOutput("feature")
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

# ── App ──────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  # Enter submits the login form. Lives at the top level so it survives the
  # login -> app swap; the guard means it does nothing once #login is gone.
  tags$head(tags$script(HTML(
    "document.addEventListener('keydown', function(e) {
       if (e.key !== 'Enter') return;
       var btn = document.getElementById('login');
       if (btn) { e.preventDefault(); btn.click(); }
     });"
  ))),
  uiOutput("page")
)

server <- function(input, output, session) {

  authed    <- reactiveVal(FALSE)
  login_msg <- reactiveVal(NULL)

  # Depends on login_msg(), so a failed attempt re-renders the login screen
  # without detaching this output from authed() — a later correct password
  # now works without a page reload.
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

  # ---- everything below is gated on authed() ----

  genes_r <- eventReactive(input$go, {
    resolve_genes(input$genes)
  }, ignoreNULL = FALSE)

  output$status <- renderUI({
    req(authed())
    g <- genes_r()
    msgs <- list()
    if (length(g$found)) {
      msgs <- c(msgs, list(tags$p(tags$b("Found: "), paste(g$found, collapse = ", "))))
    }
    if (length(g$missing)) {
      msgs <- c(msgs, list(tags$p(style = "color:#b00;",
                                  tags$b("Not in dataset: "),
                                  paste(g$missing, collapse = ", "))))
    }
    if (g$dropped > 0) {
      msgs <- c(msgs, list(tags$p(style = "color:#b00;",
                                  sprintf("Showing the first %d genes; %d more were left out.",
                                          max_genes, g$dropped))))
    }
    if (!length(g$found)) {
      msgs <- c(msgs, list(tags$p(style = "color:#b00;", "Type at least one valid gene.")))
    }
    tagList(msgs)
  })

  feat_h <- function() grid_height(length(genes_r()$found), 4, 280)
  bar_h  <- function() grid_height(length(genes_r()$found), 5, 240)
  vln_h  <- function() grid_height(length(genes_r()$found), 4, 260)

  output$ref <- renderPlot({
    req(authed())
    p_ref
  })

  output$feature <- renderPlot({
    req(authed())
    g <- genes_r()$found
    req(length(g) > 0)
    make_feature(g)
  }, height = feat_h)

  output$bar_detail <- renderPlot({
    req(authed())
    g <- genes_r()$found
    req(length(g) > 0)
    make_bar(g, "mural_class_detail", class_cols_detail,
             "Mean expression across mural subtypes",
             "All samples pooled | aSMC = aSMC + aaSMC")
  }, height = bar_h)

  output$violin_detail <- renderPlot({
    req(authed())
    g <- genes_r()$found
    req(length(g) > 0)
    make_violin(g, "mural_class_detail", class_cols_detail[class_levels_detail])
  }, height = vln_h)
}

shinyApp(ui, server)
