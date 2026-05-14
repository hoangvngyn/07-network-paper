# ============================================================
# Network Centrality in Global Value Chains and Worker Wages
# Nguyen, Viakkerev & Willems (2026) — EBC2109 Network Economics
# ============================================================

# (1) Setup
# install.packages(c("igraph", "readxl", "dplyr", "corrplot",
#                    "sandwich", "lmtest", "stargazer", "car"))
rm(list = ls())
library(igraph)
library(readxl)
library(dplyr)
library(corrplot)
library(sandwich)
library(lmtest)
library(stargazer)

setwd("/Users/arryawillems/07-network-paper")

# ---------------------------------------------------------------
# (2) Load data
load("WIOT2014_October16_ROW")
sea <- read_excel("Socio_Economic_Accounts.xlsx", sheet = "DATA")
er  <- read_excel("Exchange_Rates.xlsx", sheet = "EXR", skip = 3)

# ---------------------------------------------------------------
# (3) Build wage variable
# COMP = employee compensation (millions of national currency)
# EMPE = number of employees (thousands)
# log(COMP_usd / EMPE) = log average annual wage per employee in USD
comp <- sea[sea$variable == "COMP", c("country", "code", "2014")]
empe <- sea[sea$variable == "EMPE", c("country", "code", "2014")]

node_data <- merge(comp, empe, by = c("country", "code"))
colnames(node_data)[3:4] <- c("COMP", "EMPE")

er_2014 <- er[, c("Acronym", "_2014")]
colnames(er_2014) <- c("country", "exr_2014")
node_data <- merge(node_data, er_2014, by = "country", all.x = TRUE)

node_data$COMP_usd     <- node_data$COMP * node_data$exr_2014
node_data$wage_usd     <- node_data$COMP_usd / node_data$EMPE
node_data$log_wage_usd <- log(node_data$wage_usd)

# Drop zero/NA EMPE — China's EMPE is missing from WIOD 2016; ~127 industries
# have zero EMPE, which produces undefined log wages
node_data <- node_data[!is.na(node_data$EMPE) & node_data$EMPE > 0, ]
node_data <- node_data[is.finite(node_data$log_wage_usd), ]
cat("Wage observations after cleaning:", nrow(node_data), "\n")

# ---------------------------------------------------------------
# (4) Build the network
# Rows 1:2464 = country-industry nodes; cols 6:2469 = intermediate flows ($M)
adj_matrix <- as.matrix(wiot[1:2464, 6:2469])
node_names <- paste(wiot$Country[1:2464], wiot$IndustryCode[1:2464], sep = "-")
rownames(adj_matrix) <- node_names
colnames(adj_matrix) <- node_names

g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)
cat("Full network — nodes:", vcount(g), "edges:", ecount(g), "\n")

# ---------------------------------------------------------------
# (5) Filter edges <= $1M
# This threshold removes 89% of edges (noise below 90th percentile)
cat("Edge weight deciles:\n")
print(quantile(E(g)$weight, probs = seq(0, 1, 0.1)))
cat("Share of edges <= $1M:", round(mean(E(g)$weight <= 1), 3), "\n")

g_filtered <- delete_edges(g, E(g)[weight <= 1])
cat("Filtered network — nodes:", vcount(g_filtered), "edges:", ecount(g_filtered), "\n")

# ---------------------------------------------------------------
# (6) Centrality measures
#
# FIX 1: use strength() not degree() — the paper defines in/out-degree as
#         input/output *volumes* in $M, not edge counts
#
# FIX 2: betweenness(..., weights = NA) for structural brokerage.
#         Without this, igraph treats trade values as travel distances —
#         large flows become "long" paths, routing through tiny trades instead.
#         Structural (unweighted) betweenness correctly identifies nodes that
#         bridge otherwise disconnected parts of the production network.

str_in  <- strength(g_filtered, mode = "in")   # total value of inputs received ($M)
str_out <- strength(g_filtered, mode = "out")  # total value of inputs supplied ($M)
btw     <- betweenness(g_filtered, directed = TRUE, normalized = TRUE, weights = NA)

cat("\nIn-strength summary ($M):\n");   print(summary(str_in))
cat("Out-strength summary ($M):\n");   print(summary(str_out))
cat("Betweenness (normalized):\n");    print(summary(btw))

# ---------------------------------------------------------------
# (7) Network-level statistics (for the descriptives section)
cat("\n--- Network-Level Statistics ---\n")
cat("Density:              ", round(edge_density(g_filtered, loops = FALSE), 5), "\n")
cat("Global clustering:    ", round(transitivity(g_filtered, type = "global"), 4), "\n")
cat("Mean distance:        ", round(mean_distance(g_filtered, directed = TRUE), 3), "\n")
cat("Diameter (unweighted):", diameter(g_filtered, directed = TRUE, weights = NA), "\n")
cat("Degree assortativity: ", round(assortativity_degree(g_filtered, directed = TRUE), 4), "\n")

dist_tab <- distance_table(g_filtered, directed = TRUE)
barplot(dist_tab$res,
        names.arg = seq_along(dist_tab$res),
        xlab = "Path Length", ylab = "Node Pairs",
        main = "Distribution of Shortest Path Lengths (Filtered Network)")

# ---------------------------------------------------------------
# (8) GVC participation indices
# Computed directly from the full WIOT intermediate block (before edge filtering)
#
# Backward GVC_i = Σ_j≠country(i) Z[j,i] / GO_i
#   Share of gross output sourced from foreign intermediate inputs
#
# Forward GVC_i  = Σ_j≠country(i) Z[i,j] / GO_i
#   Share of gross output supplied as intermediate inputs to foreign industries
#
# Countries appear in contiguous blocks of 56 industries (44 × 56 = 2,464)
# so domestic vs. foreign flows can be isolated by block arithmetic.

Z          <- as.matrix(wiot[1:2464, 6:2469])
go_vals    <- as.numeric(wiot[wiot$IndustryCode == "GO", 6:2469])
ctry_vec   <- wiot$Country[1:2464]

gvc_back <- numeric(2464)
gvc_fwd  <- numeric(2464)

for (c in seq_len(44)) {
  idx     <- ((c - 1) * 56 + 1):(c * 56)   # this country's 56 industry indices
  foreign <- setdiff(seq_len(2464), idx)    # all other indices

  # Backward: foreign inputs flowing INTO this country-industry block
  gvc_back[idx] <- colSums(Z[foreign, idx]) / go_vals[idx]

  # Forward: this country-industry block's outputs flowing TO foreign industries
  gvc_fwd[idx]  <- rowSums(Z[idx, foreign]) / go_vals[idx]
}

# Cap at [0, 1] — a small number of nodes exceed 1 due to re-imports
gvc_back <- pmin(pmax(gvc_back, 0), 1)
gvc_fwd  <- pmin(pmax(gvc_fwd,  0), 1)

gvc_df <- data.frame(node_id  = node_names,
                     gvc_back = gvc_back,
                     gvc_fwd  = gvc_fwd)

# ---------------------------------------------------------------
# (9) Merge all node-level data
node_data$node_id <- paste(node_data$country, node_data$code, sep = "-")

centrality_df <- data.frame(
  node_id = names(str_in),
  str_in  = as.numeric(str_in),
  str_out = as.numeric(str_out),
  btw     = as.numeric(btw)
)

node_data <- merge(node_data, centrality_df, by = "node_id", all.x = TRUE)
node_data <- merge(node_data, gvc_df,        by = "node_id", all.x = TRUE)
cat("Final merged sample:", nrow(node_data), "observations\n")

# ---------------------------------------------------------------
# (10) Descriptive statistics
vars_desc <- c("log_wage_usd", "str_in", "str_out", "btw", "gvc_back", "gvc_fwd")
cat("\nDescriptive statistics:\n")
print(summary(node_data[, vars_desc]))

stargazer(as.data.frame(node_data[, vars_desc]),
          type = "text",
          title = "Table 1 — Descriptive Statistics",
          covariate.labels = c("Log Wage (USD)", "In-strength ($M)",
                               "Out-strength ($M)", "Betweenness centrality",
                               "Backward GVC", "Forward GVC"),
          out = "descriptives.txt")

# ---------------------------------------------------------------
# (11) Correlation matrix
cor_mat <- cor(node_data[, vars_desc], use = "pairwise.complete.obs")
cat("\nCorrelation matrix:\n")
print(round(cor_mat, 3))

corrplot(cor_mat,
         method = "color", type = "upper",
         addCoef.col = "black", tl.col = "black", tl.srt = 45,
         title = "Correlation Matrix — Node-Level Variables",
         mar = c(0, 0, 2, 0))

# ---------------------------------------------------------------
# (12) Regression analysis
# Three specifications to isolate the betweenness effect step by step:
#   M1 — Baseline: degree centrality only (trading volume)
#   M2 — Main:     add betweenness (brokerage on top of volume)
#   M3 — Full:     add GVC participation controls (separate network position
#                  from upstream/downstream integration)
#
# Heteroskedasticity-robust SEs (HC1) throughout — standard for cross-section
# wage regressions where residual variance scales with industry size.

reg_data <- node_data[complete.cases(node_data[, vars_desc]), ]
cat("\nRegression sample:", nrow(reg_data), "observations\n")

# Scale strength to billions USD for readable coefficients
reg_data$str_in_bn  <- reg_data$str_in  / 1000
reg_data$str_out_bn <- reg_data$str_out / 1000

m1 <- lm(log_wage_usd ~ str_in_bn + str_out_bn,
         data = reg_data)
m2 <- lm(log_wage_usd ~ str_in_bn + str_out_bn + btw,
         data = reg_data)
m3 <- lm(log_wage_usd ~ str_in_bn + str_out_bn + btw + gvc_back + gvc_fwd,
         data = reg_data)

# Robust standard errors
se1 <- sqrt(diag(vcovHC(m1, type = "HC1")))
se2 <- sqrt(diag(vcovHC(m2, type = "HC1")))
se3 <- sqrt(diag(vcovHC(m3, type = "HC1")))

cat("\n--- Model 1: Baseline (volume only) ---\n")
print(coeftest(m1, vcov = vcovHC(m1, type = "HC1")))
cat("\n--- Model 2: + Betweenness ---\n")
print(coeftest(m2, vcov = vcovHC(m2, type = "HC1")))
cat("\n--- Model 3: + GVC Controls ---\n")
print(coeftest(m3, vcov = vcovHC(m3, type = "HC1")))

# Multicollinearity check — VIF > 10 signals a problem
if (requireNamespace("car", quietly = TRUE)) {
  cat("\nVIF — Model 3 (concern if > 10):\n")
  print(car::vif(m3))
}

# Regression table — outputs to console and regression_table.txt
stargazer(m1, m2, m3,
          se             = list(se1, se2, se3),
          type           = "text",
          title          = "Table 2 — Network Centrality and Log Wages (Cross-Section 2014)",
          dep.var.labels = "Log Avg. Wage (USD)",
          column.labels  = c("Baseline", "Main", "Full"),
          covariate.labels = c("In-strength (bn USD)",
                               "Out-strength (bn USD)",
                               "Betweenness centrality",
                               "Backward GVC participation",
                               "Forward GVC participation"),
          add.lines = list(c("GVC controls", "No", "No", "Yes")),
          notes     = paste("Robust SEs (HC1) in parentheses.",
                            "Unit of observation: country-industry pair, 2014.",
                            "ROW and China (missing EMPE) excluded."),
          out = "regression_table.txt")

# ---------------------------------------------------------------
# (13) Network visualisation — top 30 broker nodes
top_nodes <- order(btw, decreasing = TRUE)[1:30]
g_sub     <- induced_subgraph(g_filtered, top_nodes)

ctry_sub     <- gsub("-.*", "", V(g_sub)$name)
unique_ctry  <- unique(ctry_sub)
pal          <- rainbow(length(unique_ctry))
node_colors  <- pal[match(ctry_sub, unique_ctry)]

btw_sub     <- btw[top_nodes]
node_sizes  <- btw_sub / max(btw_sub) * 20 + 5
edge_widths <- log(E(g_sub)$weight + 1) / max(log(E(g_sub)$weight + 1)) * 3

plot(g_sub,
     layout             = layout_with_fr(g_sub),
     vertex.size        = node_sizes,
     vertex.color       = node_colors,
     vertex.label       = V(g_sub)$name,
     vertex.label.cex   = 0.5,
     vertex.label.color = "black",
     edge.width         = edge_widths,
     edge.arrow.size    = 0.15,
     edge.curved        = 0.3,
     main               = "Top 30 Broker Nodes — Global Value Chain Network (2014)")

legend("bottomleft",
       legend = unique_ctry, fill = pal,
       title = "Country", cex = 0.7, bty = "n")
legend("bottomright",
       legend = c("Node size = betweenness centrality",
                  "Edge width = log trade value"),
       bty = "n", cex = 0.6)
