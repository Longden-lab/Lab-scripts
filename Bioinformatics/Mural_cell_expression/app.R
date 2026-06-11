# app.R — Mural cell gene expression explorer (password-gated)
# User must enter a password before the gene search is shown.
# The password is read from the APP_PASSWORD environment variable
# (set as a secret in Connect Cloud, NEVER committed to the repo).

library(shiny)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(httr)

# ── Fetch the data object from a PRIVATE GitHub repo at startup ──────────────

data_repo <- "Qbottle/data"                       # <- PRIVATE repo (owner/name)
data_file <- "data/mural_obj_for_suyeon.rds"      # <- path to the file inside it
rds_path  <- "mural_obj_for_suyeon.rds"           # local filename to write

if (!file.exists(rds_path)) {
  token <- Sys.getenv("DATA_REPO_TOKEN")
  if (!nzchar(token)) stop("DATA_REPO_TOKEN is not set. Add it as a secret in Connect Cloud.")
  api_url <- sprintf("https://api.github.com/repos/%s/contents/%s", data_repo, data_file)
  resp <- httr::GET(
    api_url,
    httr::add_headers(
      Authorization          = paste("Bearer", token),
      Accept                 = "application/vnd.github.raw",
      `X-GitHub-Api-Version` = "2022-11-28"
    ),
    httr::write_disk(rds_path, overwrite = TRUE),
    httr::timeout(600)
  )
  httr::stop_for_status(resp)
}

# ── Load data ONCE at startup (not per request) ──────────────────────────────
mural_obj <- readRDS(rds_path)
DefaultAssay(mural_obj) <- "RNA"
stopifnot("umap" %in% Reductions(mural_obj))

# ── Class definitions (labels live in mural_final) ───────────────────────────
class_map    <- c(aSMC="aSMC", aaSMC="aSMC", C_PC="PC", Ts_PC="PC", vSMC="vSMC")
class_levels <- c("aSMC","PC","vSMC")
class_cols   <- c(aSMC="firebrick", PC="darkorange2", vSMC="steelblue")

class_map_detail    <- c(aSMC="aSMC", aaSMC="aSMC", C_PC="C_PC", Ts_PC="Ts_PC", vSMC="vSMC")
class_levels_detail <- c("aSMC","C_PC","Ts_PC","vSMC")
class_cols_detail   <- c(aSMC="firebrick", C_PC="darkorange2", Ts_PC="goldenrod3", vSMC="steelblue")

mu <- mural_obj
mu$mural_class        <- factor(unname(class_map[as.character(mu$mural_final)]),        levels = class_levels)
mu$mural_class_detail <- factor(unname(class_map_detail[as.character(mu$mural_final)]), levels = class_levels_detail)

all_genes <- sort(rownames(mu))

# Reference UMAP — same for every query, so build it once
p_ref <- DimPlot(mu, group.by = "mural_final", reduction = "umap", label = TRUE,
                 repel = TRUE, label.size = 4, pt.size = 0.5) +
  labs(title = "Mural subtypes (reference)") +
  theme(plot.title = element_text(face = "bold", size = 12))

# ── Plot builders (return ggplot objects instead of saving to disk) ──────────
make_bar <- function(genes, class_col, cols, title, subtitle = NULL) {
  facet_ncol <- min(5, length(genes))
  df <- FetchData(mu, vars = c(genes, class_col), layer = "data") %>%
    pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") %>%
    group_by(gene, .data[[class_col]]) %>%
    summarise(mean = mean(expr), sem = sd(expr) / sqrt(n()), .groups = "drop") %>%
    mutate(gene = factor(gene, levels = genes))
  names(df)[2] <- "grp"
  ggplot(df, aes(grp, mean, fill = grp)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = mean, ymax = mean + sem), width = 0.2, color = "grey30") +
    geom_text(aes(label = round(mean, 2), y = mean + sem), vjust = -0.4, size = 2.8, color = "grey20") +
    facet_wrap(~ gene, ncol = facet_ncol, scales = "free_y") +
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Mean expression (log-norm)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          strip.text = element_text(face = "bold.italic"),
          legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
}

make_feature <- function(genes) {
  feat_ncol <- min(4, length(genes))
  FeaturePlot(mu, features = genes, reduction = "umap", order = TRUE,
              pt.size = 0.5, cols = c("lightgrey", "firebrick"), ncol = feat_ncol) &
    theme(plot.title = element_text(size = 11, face = "bold.italic"),
          axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
}

make_violin <- function(genes) {
  vln_ncol <- min(4, length(genes))
  VlnPlot(mu, features = genes, group.by = "mural_class", pt.size = 0,
          cols = class_cols[class_levels], ncol = vln_ncol) &
    theme(plot.title = element_text(size = 11, face = "bold.italic"),
          axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5))
}

# ── The real app UI, shown only after a correct password ─────────────────────
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
        tags$small(sprintf("%d genes available in this dataset.", length(all_genes)))
      ),
      mainPanel(
        width = 9,
        h4("Reference UMAP"),
        plotOutput("ref", height = "420px"),
        h4("Feature plot (UMAP)"),
        plotOutput("feature"),
        h4("Bar graph — 3 classes (aSMC / PC / vSMC)"),
        plotOutput("bar"),
        h4("Bar graph — subtypes (PC split: C_PC / Ts_PC)"),
        plotOutput("bar_detail"),
        h4("Violin — 3 classes"),
        plotOutput("violin")
      )
    )
  )
}

# ── The login screen ─────────────────────────────────────────────────────────
login_ui <- function(msg = NULL) {
  div(style = "max-width:340px; margin:80px auto; text-align:center;",
      h3("Mural cell gene expression"),
      p("Enter the password to continue."),
      passwordInput("pw", NULL, placeholder = "Password"),
      actionButton("login", "Enter", class = "btn-primary"),
      if (!is.null(msg)) tags$p(style = "color:#b00; margin-top:12px;", msg)
  )
}

# ── App ───────────────────────────────────────────────────────────────────────
ui <- fluidPage(uiOutput("page"))

server <- function(input, output, session) {

  authed <- reactiveVal(FALSE)

  # Show login screen or the real app depending on auth state
  output$page <- renderUI({
    if (authed()) main_ui() else login_ui()
  })

  observeEvent(input$login, {
    pw_set <- Sys.getenv("APP_PASSWORD")
    if (nzchar(pw_set) && identical(input$pw, pw_set)) {
      authed(TRUE)
    } else {
      output$page <- renderUI({ login_ui("Incorrect password.") })
    }
  })

  # ---- everything below only matters once the real UI is on screen ----

  genes_r <- eventReactive(input$go, {
    raw <- unlist(strsplit(input$genes, "[,\\s]+"))
    raw <- trimws(raw)
    raw <- raw[nzchar(raw)]
    found   <- raw[raw %in% all_genes]
    missing <- setdiff(raw, all_genes)
    list(found = unique(found), missing = unique(missing))
  }, ignoreNULL = FALSE)

  output$status <- renderUI({
    req(authed())
    g <- genes_r()
    msgs <- list()
    if (length(g$found))
      msgs <- c(msgs, list(tags$p(tags$b("Found: "), paste(g$found, collapse = ", "))))
    if (length(g$missing))
      msgs <- c(msgs, list(tags$p(style = "color:#b00;",
                                  tags$b("Not in dataset: "), paste(g$missing, collapse = ", "))))
    if (!length(g$found))
      msgs <- c(msgs, list(tags$p(style = "color:#b00;", "Type at least one valid gene.")))
    tagList(msgs)
  })

  feat_h <- function() { g <- genes_r()$found; 280 * ceiling(length(g) / min(4, max(1, length(g)))) }
  bar_h  <- function() { g <- genes_r()$found; 240 * ceiling(length(g) / min(5, max(1, length(g)))) }
  vln_h  <- function() { g <- genes_r()$found; 260 * ceiling(length(g) / min(4, max(1, length(g)))) }

  output$ref <- renderPlot({ req(authed()); p_ref })

  output$feature <- renderPlot({
    req(authed()); g <- genes_r()$found; req(length(g) > 0); make_feature(g)
  }, height = feat_h)

  output$bar <- renderPlot({
    req(authed()); g <- genes_r()$found; req(length(g) > 0)
    make_bar(g, "mural_class", class_cols,
             "Mean expression across mural classes (all samples pooled)")
  }, height = bar_h)

  output$bar_detail <- renderPlot({
    req(authed()); g <- genes_r()$found; req(length(g) > 0)
    make_bar(g, "mural_class_detail", class_cols_detail,
             "Mean expression across mural subtypes",
             "All samples pooled | aSMC = aSMC + aaSMC")
  }, height = bar_h)

  output$violin <- renderPlot({
    req(authed()); g <- genes_r()$found; req(length(g) > 0); make_violin(g)
  }, height = vln_h)
}

shinyApp(ui, server)
