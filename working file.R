# ------------------------------------------------------------- 
# (1) Setup - install + load packages and clean the environment
# install.packages(c("igraph", "readxl", "dplyr", "corrplot"))
rm(list = ls()) 
library(igraph)
library(readxl)
library(dplyr)
library(corrplot)

# -------------------------------------------------------------
# (2) Set working directory and then load data — WIOT and SEA
setwd("/Users/hoangvietnguyen/Documents/06 projects /07 network project")
load("/Users/hoangvietnguyen/Documents/06 projects /07 network project/WIOTS_in_R/WIOT2014_October16_ROW.RData")
sea <- read_excel("Socio_Economic_Accounts.xlsx", sheet = "DATA")
er  <- read_excel("Exchange_Rates.xlsx", sheet = "EXR", skip = 3)

# -------------------------------------------------------------
# (3) Build wage variable
# COMP = employee compensation (millions of national currency)
# EMPE = number of employees (thousands)
# wage = COMP / EMPE = average annual compensation per employee
# log(wage) used because wage distributions are right-skewed (Mincer, 1974)
comp <- sea[sea$variable == "COMP", c("country", "code", "2014")]
empe <- sea[sea$variable == "EMPE", c("country", "code", "2014")]

node_data <- merge(comp, empe, by = c("country", "code"))
colnames(node_data)[3:4] <- c("COMP", "EMPE")

node_data$wage     <- node_data$COMP / node_data$EMPE
node_data$log_wage <- log(node_data$wage)

# Convert COMP to USD for cross-country comparability
# Without conversion, comparing Indonesian rupiah to Swiss franc is meaningless
er_2014 <- er[, c("Acronym", "_2014")]
colnames(er_2014) <- c("country", "exr_2014")
node_data <- merge(node_data, er_2014, by = "country", all.x = TRUE)
node_data$COMP_usd     <- node_data$COMP * node_data$exr_2014
node_data$wage_usd     <- node_data$COMP_usd / node_data$EMPE
node_data$log_wage_usd <- log(node_data$wage_usd)

summary(node_data)

# -------------------------------------------------------------
# (4) Build the network
# Nodes = country-industry pairs (e.g. AUS-A01, DEU-C28)
# Edges = intermediate input flows (directed, weighted)
# Rows 1:2464 = country-industry pairs; cols 6:2469 = intermediate flows only
adj_matrix <- as.matrix(wiot[1:2464, 6:2469])
node_names <- paste(wiot$Country[1:2464], wiot$IndustryCode[1:2464], sep = "-")
rownames(adj_matrix) <- node_names
colnames(adj_matrix) <- node_names

g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)
summary(g)

# Filter edges below $1 million — 89% of edges are noise below this threshold
# $1M sits at the 90th percentile of the edge weight distribution
quantile(E(g)$weight, probs = seq(0, 1, 0.1))
mean(E(g)$weight <= 1)
g_filtered <- delete_edges(g, E(g)[weight <= 1])
summary(g_filtered)

# -------------------------------------------------------------
# (5) Network statistics
## Density: fraction of all possible edges that exist (sparse = meaningful variation)
edge_density(g_filtered, loops = FALSE)

## Degree: in = inputs received (buyer), out = inputs sent (supplier)
deg_in  <- degree(g_filtered, mode = "in")
deg_out <- degree(g_filtered, mode = "out")
summary(deg_in)
summary(deg_out)

## Betweenness: how often a node sits on shortest paths — our broker measure
btw <- betweenness(g_filtered, directed = TRUE, normalized = TRUE)
summary(btw)

## Clustering: probability that two neighbors are also connected to each other
clust <- transitivity(g_filtered, type = "global")
clust

## Average path length: average steps between any two nodes (network level)
apl <- mean_distance(g_filtered, directed = TRUE)
apl

## Diameter: longest shortest path (worst case reachability, unweighted)
diam <- diameter(g_filtered, directed = TRUE, weights = NA)
diam

## Distance table: distribution of path lengths across all node pairs
dist_table <- distance_table(g_filtered, directed = TRUE)
dist_table$res
barplot(dist_table$res,
        names.arg = 1:length(dist_table$res),
        xlab = "Path Length",
        ylab = "Number of Node Pairs",
        main = "Distribution of Path Lengths")

## Assortativity: do hubs connect to hubs? Negative = hub-and-spoke structure
assort <- assortativity_degree(g_filtered, directed = TRUE)
assort

## Eigenvector centrality: importance weighted by neighbor importance
eigen_cent <- eigen_centrality(g_filtered, directed = TRUE)$vector
summary(eigen_cent)

## Closeness: how quickly a node can reach all others (not used in regression)
close_cent <- closeness(g_filtered, mode = "all", normalized = TRUE)
hist(close_cent)

# -------------------------------------------------------------
# (6) Merge centrality with wage data
node_data$node_id <- paste(node_data$country, node_data$code, sep = "-")

centrality_df <- data.frame(
  node_id = names(deg_in),
  deg_in  = as.numeric(deg_in),
  deg_out = as.numeric(deg_out),
  btw     = as.numeric(btw)
)

node_data <- merge(node_data, centrality_df, by = "node_id", all.x = TRUE)
dim(node_data)

# -------------------------------------------------------------
# (7) Descriptive statistics table
summary(node_data[, c("log_wage_usd", "deg_in", "deg_out", "btw")])

# -------------------------------------------------------------
# (8) Correlation matrix
cor_matrix <- cor(node_data[, c("log_wage_usd", "deg_in", "deg_out", "btw")],
                  use = "pairwise.complete.obs")
round(cor_matrix, 3)

corrplot(cor_matrix,
         method = "color",
         type = "upper",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         title = "Correlation Matrix — Node-Level Variables",
         mar = c(0, 0, 1, 0))

# -------------------------------------------------------------
# (9) Network plot — top 30 broker nodes
top_nodes <- order(btw, decreasing = TRUE)[1:30]
g_sub <- induced_subgraph(g_filtered, top_nodes)

countries        <- gsub("-.*", "", V(g_sub)$name)
unique_countries <- unique(countries)
colors           <- rainbow(length(unique_countries))
node_colors      <- colors[match(countries, unique_countries)]

btw_sub    <- btw[top_nodes]
node_sizes <- btw_sub / max(btw_sub) * 20 + 5

edge_weights <- log(E(g_sub)$weight + 1)
edge_widths  <- edge_weights / max(edge_weights) * 3

plot(g_sub,
     layout = layout_with_fr(g_sub),
     vertex.size = node_sizes,
     vertex.color = node_colors,
     vertex.label = V(g_sub)$name,
     vertex.label.cex = 0.5,
     vertex.label.color = "black",
     edge.width = edge_widths,
     edge.arrow.size = 0.15,
     edge.curved = 0.3,
     main = "Top 30 Broker Nodes — Global Value Chain Network (2014)")

legend("bottomleft", legend = unique_countries, fill = colors,
       title = "Country", cex = 0.7, bty = "n")
legend("bottomright",
       legend = c("Node size = betweenness centrality",
                  "Edge width = trade flow value"),
       bty = "n", cex = 0.6)