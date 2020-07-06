library("DEP")
library("dplyr")
library("readxl")

#get the data
data <- read_excel("C:/Users/thilo/Documents/data.xlsx")

#check for douplicates
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)

data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

#Get a SummarizedExperiment object
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
experimental_design <- read_excel("C:/Users/thilo/Documents/data_se.xlsx")
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

#for plotting:
data_filt <- filter_missval(data_se, thr = 0)
data_filt2 <- filter_missval(data_se, thr = 1)

data_norm <- normalize_vsn(data_filt)

data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

# chapter 4.9 is working now:
data_diff <- test_diff(data_imp, type = "control", control = "non.infected")
data_diff_all_contrasts <- test_diff(data_imp, type = "all")

dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))

#Plot

plot_frequency(data_se)

plot_numbers(data_filt)

plot_coverage(data_filt)

plot_normalization(data_filt, data_norm)

plot_missval(data_filt)

plot_detect(data_filt)

plot_imputation(data_norm, data_imp)

#plot chapter 4.9

plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)

plot_volcano(dep, contrast = "infected_vs_non.infected", label_size = 2, add_names = TRUE)

plot_single(dep, proteins = c("Fabp1", "Hbb-b1"))

plot_single(dep, proteins = "Fabp1", type = "centered")

#

data_results <- get_results(dep)

df_wide <- get_df_wide(dep)

df_long <- get_df_long(dep)

save(data_se, data_norm, data_imp, data_diff, dep, file = "data.RData")

#load again:
load("data.RData")

