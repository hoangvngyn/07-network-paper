# ============================================================
# Network Centrality in Global Value Chains and Worker Wages
# Nguyen, Viakkerev & Willems (2026) — EBC2109 Network Economics
# ============================================================

# (1) Setup
# install.packages(c("igraph", "readxl", "dplyr", "corrplot",
#                    "sandwich", "lmtest", "stargazer", "car", "ggplot2", "scales"))
rm(list = ls())
library(igraph)
library(readxl)
library(dplyr)
library(corrplot)
library(sandwich)
library(lmtest)
library(stargazer)
library(ggplot2)
library(scales)

setwd("/Users/arryawillems/07-network-paper")
dir.create("plots", showWarnings = FALSE)

# ---------------------------------------------------------------
# (2) Load data
load("WIOT2014_October16_ROW")
sea <- read_excel("Socio_Economic_Accounts.xlsx", sheet = "DATA")
er  <- read_excel("Exchange_Rates.xlsx", sheet = "EXR", skip = 3)

# ---------------------------------------------------------------
# (3) Build wage variable
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

node_data <- node_data[!is.na(node_data$EMPE) & node_data$EMPE > 0, ]
node_data <- node_data[is.finite(node_data$log_wage_usd), ]
cat("Wage observations after cleaning:", nrow(node_data), "\n")

# ---------------------------------------------------------------
# (4) Build the network
adj_matrix <- as.matrix(wiot[1:2464, 6:2469])
node_names <- paste(wiot$Country[1:2464], wiot$IndustryCode[1:2464], sep = "-")
rownames(adj_matrix) <- node_names
colnames(adj_matrix) <- node_names

g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)
cat("Full network — nodes:", vcount(g), "edges:", ecount(g), "\n")

# ---------------------------------------------------------------
# (5) Filter edges <= $1M
cat("Edge weight deciles:\n")
print(quantile(E(g)$weight, probs = seq(0, 1, 0.1)))
cat("Share of edges <= $1M:", round(mean(E(g)$weight <= 1), 3), "\n")

g_filtered <- delete_edges(g, E(g)[weight <= 1])
cat("Filtered network — nodes:", vcount(g_filtered), "edges:", ecount(g_filtered), "\n")

# ---------------------------------------------------------------
# (6) Centrality measures
str_in  <- strength(g_filtered, mode = "in")
str_out <- strength(g_filtered, mode = "out")
btw     <- betweenness(g_filtered, directed = TRUE, normalized = TRUE, weights = NA)

# ---------------------------------------------------------------
# (7) Network-level statistics
cat("\n--- Network-Level Statistics ---\n")
net_density  <- edge_density(g_filtered, loops = FALSE)
net_clust    <- transitivity(g_filtered, type = "global")
net_apl      <- mean_distance(g_filtered, directed = TRUE)
net_diam     <- diameter(g_filtered, directed = TRUE, weights = NA)
net_assort   <- assortativity_degree(g_filtered, directed = TRUE)

cat("Density:              ", round(net_density, 5), "\n")
cat("Global clustering:    ", round(net_clust, 4), "\n")
cat("Mean distance:        ", round(net_apl, 3), "\n")
cat("Diameter (unweighted):", net_diam, "\n")
cat("Degree assortativity: ", round(net_assort, 4), "\n")

dist_tab <- distance_table(g_filtered, directed = TRUE)

# FIGURE 1 — Path length distribution
png("plots/01_path_length_distribution.png", width = 900, height = 580, res = 150)
par(mar = c(4.5, 4.5, 3, 1.5), bg = "white")
bp <- barplot(dist_tab$res / 1e6,
              names.arg = seq_along(dist_tab$res),
              col = "#3B82F6", border = "white",
              xlab = "Shortest Path Length (hops)",
              ylab = "Node Pairs (millions)",
              main = "Distribution of Shortest Path Lengths",
              cex.main = 1.1, cex.lab = 0.95, cex.axis = 0.85,
              ylim = c(0, max(dist_tab$res / 1e6) * 1.15))
text(bp, dist_tab$res / 1e6 + max(dist_tab$res / 1e6) * 0.03,
     labels = paste0(round(dist_tab$res / 1e6, 1), "M"),
     cex = 0.75, col = "#1e3a5f")
mtext(paste0("Mean = ", round(net_apl, 2), " hops  |  Diameter = ", net_diam,
             "  |  Density = ", round(net_density * 100, 1), "%"),
      side = 3, line = 0.1, cex = 0.75, col = "gray40")
dev.off()

# ---------------------------------------------------------------
# (8) GVC participation indices
Z        <- as.matrix(wiot[1:2464, 6:2469])
go_vals  <- as.numeric(wiot[wiot$IndustryCode == "GO", 6:2469])

gvc_back <- numeric(2464)
gvc_fwd  <- numeric(2464)

for (c in seq_len(44)) {
  idx     <- ((c - 1) * 56 + 1):(c * 56)
  foreign <- setdiff(seq_len(2464), idx)
  gvc_back[idx] <- colSums(Z[foreign, idx]) / go_vals[idx]
  gvc_fwd[idx]  <- rowSums(Z[idx, foreign]) / go_vals[idx]
}

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

# FIGURE 2 — Correlation heatmap
png("plots/02_correlation_matrix.png", width = 800, height = 780, res = 150)
par(oma = c(0, 0, 2, 0), bg = "white")
corrplot(cor_mat,
         method       = "color",
         type         = "upper",
         addCoef.col  = "black",
         number.cex   = 0.82,
         tl.col       = "black",
         tl.srt       = 40,
         tl.cex       = 0.85,
         col          = colorRampPalette(c("#2563EB", "white", "#DC2626"))(200),
         cl.cex       = 0.75,
         mar          = c(0, 0, 1, 0))
mtext("Correlation Matrix — Node-Level Variables",
      side = 3, outer = TRUE, cex = 1.0, font = 2, line = 0.5)
dev.off()

# ---------------------------------------------------------------
# (12) Regression analysis
reg_data <- node_data[complete.cases(node_data[, vars_desc]), ]
cat("\nRegression sample:", nrow(reg_data), "observations\n")

reg_data$str_in_bn  <- reg_data$str_in  / 1000
reg_data$str_out_bn <- reg_data$str_out / 1000

m1 <- lm(log_wage_usd ~ str_in_bn + str_out_bn,
         data = reg_data)
m2 <- lm(log_wage_usd ~ str_in_bn + str_out_bn + btw,
         data = reg_data)
m3 <- lm(log_wage_usd ~ str_in_bn + str_out_bn + btw + gvc_back + gvc_fwd,
         data = reg_data)

se1 <- sqrt(diag(vcovHC(m1, type = "HC1")))
se2 <- sqrt(diag(vcovHC(m2, type = "HC1")))
se3 <- sqrt(diag(vcovHC(m3, type = "HC1")))

cat("\n--- Model 1: Baseline ---\n")
print(coeftest(m1, vcov = vcovHC(m1, type = "HC1")))
cat("\n--- Model 2: + Betweenness ---\n")
print(coeftest(m2, vcov = vcovHC(m2, type = "HC1")))
cat("\n--- Model 3: + GVC Controls ---\n")
print(coeftest(m3, vcov = vcovHC(m3, type = "HC1")))

if (requireNamespace("car", quietly = TRUE)) {
  cat("\nVIF — Model 3:\n")
  print(car::vif(m3))
}

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
                            "Unit: country-industry pair, 2014.",
                            "ROW and China (missing EMPE) excluded."),
          out = "regression_table.txt")

# ---------------------------------------------------------------
# (13) Figures for README

# --- Country region groupings for wage distribution plot
high_income  <- c("AUS","AUT","BEL","CAN","CHE","CYP","DEU","DNK","ESP","EST",
                  "FIN","FRA","GBR","GRC","HRV","IRL","ITA","JPN","KOR","LTU",
                  "LUX","LVA","MLT","NLD","NOR","PRT","SVK","SVN","SWE","TWN","USA")
upper_middle <- c("BGR","BRA","CHN","MEX","ROU","RUS","TUR")
lower_middle <- c("IDN","IND")

reg_data$region <- ifelse(reg_data$country %in% high_income,  "High income",
                   ifelse(reg_data$country %in% upper_middle, "Upper-middle income",
                                                              "Lower-middle income"))
reg_data$region <- factor(reg_data$region,
                          levels = c("High income", "Upper-middle income", "Lower-middle income"))

# FIGURE 3 — Wage distribution by income group
pal3 <- c("High income" = "#2563EB", "Upper-middle income" = "#F59E0B",
          "Lower-middle income" = "#10B981")

p_wage <- ggplot(reg_data, aes(x = log_wage_usd, fill = region)) +
  geom_histogram(binwidth = 0.25, colour = "white", linewidth = 0.2, alpha = 0.9,
                 position = "identity") +
  scale_fill_manual(values = pal3, name = NULL) +
  scale_x_continuous(breaks = seq(-2, 7, 1)) +
  labs(title    = "Distribution of Log Wages by Income Group",
       subtitle = "Country-industry level, 2014  (n = 2,171)",
       x        = "Log Average Wage (USD)",
       y        = "Count") +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(colour = "gray40", size = 9),
        legend.position = "top",
        panel.grid.minor = element_blank())

ggsave("plots/03_wage_distribution.png", p_wage,
       width = 8, height = 5, dpi = 150, bg = "white")

# FIGURE 4 — Betweenness vs log wage scatter
p_btw <- ggplot(reg_data, aes(x = btw, y = log_wage_usd, colour = region)) +
  geom_point(alpha = 0.35, size = 1.0) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              colour = "#1e3a5f", fill = "#93C5FD", linewidth = 1) +
  scale_colour_manual(values = pal3, name = NULL) +
  scale_x_continuous(labels = label_scientific()) +
  labs(title    = "Betweenness Centrality vs. Log Wage",
       subtitle = "OLS trend line with 95% CI  |  Structural brokerage predicts higher wages",
       x        = "Betweenness Centrality (normalized)",
       y        = "Log Average Wage (USD)") +
  theme_minimal(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(colour = "gray40", size = 9),
        legend.position = "top",
        panel.grid.minor = element_blank())

ggsave("plots/04_betweenness_vs_wage.png", p_btw,
       width = 8, height = 5, dpi = 150, bg = "white")

# FIGURE 5 — GVC participation vs log wage (two panels)
gvc_long <- rbind(
  data.frame(log_wage_usd = reg_data$log_wage_usd,
             gvc          = reg_data$gvc_back,
             type         = "Backward GVC participation",
             region       = reg_data$region),
  data.frame(log_wage_usd = reg_data$log_wage_usd,
             gvc          = reg_data$gvc_fwd,
             type         = "Forward GVC participation",
             region       = reg_data$region)
)

p_gvc <- ggplot(gvc_long, aes(x = gvc, y = log_wage_usd, colour = region)) +
  geom_point(alpha = 0.25, size = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              colour = "#1e3a5f", fill = "#93C5FD", linewidth = 1) +
  scale_colour_manual(values = pal3, name = NULL) +
  facet_wrap(~ type) +
  labs(title    = "GVC Participation vs. Log Wage",
       subtitle = "Both backward and forward integration predict higher wages",
       x        = "GVC Participation Index (0–1)",
       y        = "Log Average Wage (USD)") +
  theme_minimal(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(colour = "gray40", size = 9),
        legend.position  = "top",
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"))

ggsave("plots/05_gvc_vs_wage.png", p_gvc,
       width = 10, height = 5, dpi = 150, bg = "white")

# FIGURE 6 — Standardised coefficient forest plot (Model 3)
reg_std <- as.data.frame(scale(
  reg_data[, c("log_wage_usd","str_in_bn","str_out_bn","btw","gvc_back","gvc_fwd")]
))

m3_std   <- lm(log_wage_usd ~ str_in_bn + str_out_bn + btw + gvc_back + gvc_fwd,
               data = reg_std)
m3_coef  <- coeftest(m3_std, vcov = vcovHC(m3_std, type = "HC1"))

coef_plot_df <- data.frame(
  term     = factor(c("In-strength", "Out-strength", "Betweenness",
                      "Backward GVC", "Forward GVC"),
                    levels = rev(c("In-strength", "Out-strength", "Betweenness",
                                   "Backward GVC", "Forward GVC"))),
  estimate = m3_coef[-1, 1],
  se       = m3_coef[-1, 2],
  pval     = m3_coef[-1, 4]
)
coef_plot_df$ci_lo  <- coef_plot_df$estimate - 1.96 * coef_plot_df$se
coef_plot_df$ci_hi  <- coef_plot_df$estimate + 1.96 * coef_plot_df$se
coef_plot_df$sig    <- ifelse(coef_plot_df$pval < 0.01,  "p < 0.01",
                       ifelse(coef_plot_df$pval < 0.05,  "p < 0.05",
                       ifelse(coef_plot_df$pval < 0.1,   "p < 0.10",
                                                          "n.s.")))
coef_plot_df$sig    <- factor(coef_plot_df$sig,
                              levels = c("p < 0.01","p < 0.05","p < 0.10","n.s."))

pal_sig <- c("p < 0.01" = "#1D4ED8", "p < 0.05" = "#2563EB",
             "p < 0.10" = "#60A5FA", "n.s."     = "#CBD5E1")

p_coef <- ggplot(coef_plot_df, aes(x = estimate, y = term, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray60", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.9) +
  geom_point(size = 3.5) +
  scale_colour_manual(values = pal_sig, name = "Significance") +
  labs(title    = "Standardised Regression Coefficients — Full Model (M3)",
       subtitle = "Coefficients show effect of 1 SD change in predictor on log wage (SD units)  |  Error bars: 95% CI",
       x        = "Standardised Coefficient (β)",
       y        = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(colour = "gray40", size = 9),
        legend.position  = "right",
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_text(size = 10.5))

ggsave("plots/06_regression_coefficients.png", p_coef,
       width = 9, height = 5, dpi = 150, bg = "white")

# FIGURE 7 — R-squared progression across the three models
r2_df <- data.frame(
  model   = factor(c("M1: Volume\n(in/out-strength)",
                     "M2: + Betweenness",
                     "M3: + GVC\ncontrols"),
                   levels = c("M1: Volume\n(in/out-strength)",
                              "M2: + Betweenness",
                              "M3: + GVC\ncontrols")),
  r2      = c(summary(m1)$r.squared,
              summary(m2)$r.squared,
              summary(m3)$r.squared),
  r2_adj  = c(summary(m1)$adj.r.squared,
              summary(m2)$adj.r.squared,
              summary(m3)$adj.r.squared)
)

p_r2 <- ggplot(r2_df, aes(x = model, y = r2)) +
  geom_col(fill = "#3B82F6", width = 0.5, alpha = 0.9) +
  geom_text(aes(label = paste0("R² = ", round(r2, 3))),
            vjust = -0.6, fontface = "bold", size = 3.5, colour = "#1e3a5f") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     limits = c(0, max(r2_df$r2) * 1.25)) +
  labs(title    = "Explained Variance by Model Specification",
       subtitle = "Adding betweenness centrality (M2) and GVC controls (M3) each improve fit",
       x        = NULL,
       y        = "R-squared") +
  theme_minimal(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(colour = "gray40", size = 9),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank())

ggsave("plots/07_r2_progression.png", p_r2,
       width = 7, height = 4.5, dpi = 150, bg = "white")

# ---------------------------------------------------------------
# (14) Network visualisation — top 30 broker nodes
top_nodes <- order(btw, decreasing = TRUE)[1:30]
g_sub     <- induced_subgraph(g_filtered, top_nodes)

ctry_sub    <- gsub("-.*", "", V(g_sub)$name)
unique_ctry <- unique(ctry_sub)
pal_net     <- rainbow(length(unique_ctry), s = 0.7, v = 0.85)
node_colors <- pal_net[match(ctry_sub, unique_ctry)]

btw_sub     <- btw[top_nodes]
node_sizes  <- btw_sub / max(btw_sub) * 22 + 6
edge_widths <- log(E(g_sub)$weight + 1) / max(log(E(g_sub)$weight + 1)) * 3

# FIGURE 8 — Network visualisation
png("plots/08_network_top30_brokers.png", width = 1400, height = 1050, res = 150)
par(mar = c(1, 1, 2, 1), bg = "white")
set.seed(42)
plot(g_sub,
     layout             = layout_with_fr(g_sub),
     vertex.size        = node_sizes,
     vertex.color       = node_colors,
     vertex.frame.color = "white",
     vertex.label       = V(g_sub)$name,
     vertex.label.cex   = 0.55,
     vertex.label.color = "black",
     vertex.label.dist  = 0,
     edge.width         = edge_widths,
     edge.arrow.size    = 0.2,
     edge.color         = rgb(0.4, 0.4, 0.4, 0.5),
     edge.curved        = 0.25,
     main               = "Top 30 Broker Nodes in the Global Value Chain Network (2014)")
legend("bottomleft",
       legend = unique_ctry, fill = pal_net,
       title = "Country", cex = 0.65, bty = "n", ncol = 2)
legend("bottomright",
       legend = c("Node size = betweenness centrality",
                  "Edge width = log trade value ($M)"),
       bty = "n", cex = 0.65)
dev.off()

cat("\nAll plots saved to plots/\n")
cat("Regression tables saved to descriptives.txt and regression_table.txt\n")
