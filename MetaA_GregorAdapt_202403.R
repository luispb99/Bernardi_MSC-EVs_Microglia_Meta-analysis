#############################
# Package Install/Load
#############################
## DOI: 10.1016/j.brs.2021.05.014
library(metafor)
library(dplyr)
library(metaviz)
library(ggplot2)

#############################
# Data
#############################
## Data.frame ##
### Columns description ###
#### - study: list of different studies
#### - brain_region: list of brain regions in each study
#### - n_c: number of samples in reference group for each study
#### - m_c: mean of samples in reference group for each study
#### - sd_c: standard deviation in reference group for each study
#### - n_t: number of samples in test group for each study
#### - m_t: mean of samples in test group for each study
#### - sd_t: standard deviation in test group for each study
#### - n_s: number of subjects in study
colnames(data_original) <- tolower(colnames(data_original))

data_split <- split(x = data_original, f = data_original$table.id)

resList <- lapply(data_split, function(data) {
     #############################
     # Calculate Effect Sizes
     #############################
     data_escalc <- escalc(n1i = n_t, m1i = m_t, sd1i = sd_t, n2i = n_c, m2i = m_c, sd2i = sd_c, measure = "SMD",
                           data = data, vtype = "UB", append = F) 
     # I use vtype = UB in my publications, it's more conservative
     
     data <- cbind(data, data_escalc)
     
     dup <- data$study[duplicated(data$study)]
     
     data_dup <- data[data$study%in%dup,]
     data_unique <- data[!data$study%in%dup,]
     row.names(data_unique) <- data_unique$study
     
     if(nrow(data_dup)>0) {
          n_samples <- split(x = data_dup, f = data_dup$study)
          n_samples <- lapply(n_samples, function(x) c(sum(x$n_c), sum(x$n_t)))
          n_samples <- do.call("rbind", n_samples)
          colnames(n_samples) <- c("n_c", "n_t")
          n_samples <- rbind(n_samples, data_unique[,c("n_c","n_t")])
          colnames(n_samples) <- c("n_samples_c", "n_samples_t")
     } else {
          n_samples <- data_unique[,c("n_c","n_t")]
          colnames(n_samples) <- c("n_samples_c", "n_samples_t")
     }
     
     #############################
     # Combined SMD (adjusted standardized mean difference by fixed-effect model)
     #############################
     studs <- unique(data_dup$study)
     res <- lapply(studs, function(x, dt) {
          ourdata <- dt[dt$study==x,]
          z <- rma.uni(yi = yi, vi = vi, measure = "SMD",
                       data = ourdata, method = "FE")
          c(z[1],x)
     }, dt = data_dup)
     data_comb <- data_dup[!duplicated(data_dup$study),]
     row.names(data_comb) <- data_comb$study
     for(s in studs) {
          data_comb[s, "n_c"] <- max(data_comb[s, "n_c"])
          data_comb[s, "n_t"] <- max(data_comb[s, "n_t"])
     }
     data_comb$yi <- as.numeric(lapply(res, '[[',1))
     data_comb$vi <- ((data_comb$n_c + data_comb$n_t)/(data_comb$n_c * data_comb$n_t) + 
                           data_comb$yi^2)/(2*(data_comb$n_c + data_comb$n_t))
     # here you have to be consistent with the vtype above 
     # (you can find the formula corresponding to your outcome measure in the escalc function script, 
     # for SMD and vtype UB  it is: (n1+n2)/(n1*n2) +y^2 / (2*(n1+n2)) )
     
     #############################
     # Random-effect model
     #############################
     data_final <- rbind(data_unique, data_comb)
     data_final$ci_low <- data_final$yi - (stats::qnorm(1 - (1 - 0.95)/2) * as.numeric(sqrt(data_final$vi)))
     data_final$ci_high <- data_final$yi + (stats::qnorm(1 - (1 - 0.95)/2) * as.numeric(sqrt(data_final$vi)))
     data_final <- cbind(data_final, n_samples[data_final$study,])
     row.names(data_final) <- data_final$id
     data_final <- data_final[as.character(unique(data$id)),]
     
     res <- rma.uni(yi, vi, data = data_final, method="DL", 
                    slab = paste(data_final$stud,
                                 sep=""))
     return(res)
})

table_indiv <- lapply(resList, function(x) {
     d <- x$data
     w <- weights(x)
     d$weights <- w[d$study]
     return(d)
})
table_indiv <- do.call("rbind", table_indiv)
row.names(table_indiv) <- NULL

table_synthesis <- lapply(resList, function(x) {
     w <- weights(x)
     fstats <- fitstats(x)[,1]
     names(fstats) <- gsub(pattern = ":", replacement = "", x = names(fstats))
     c("SMD" = round(x$beta, 2), "CI_lower" = round(x$ci.lb, 2), 
       "CI_upper" = round(x$ci.ub, 2), "pvalue" = round(x$pval,3), "zvalue" = round(x$zval,3),
       "I2" = round(x$I2, 2), "tau2" = round(x$tau2, 2), round(fstats, 3),
       # "weights" = paste(paste(names(w), round(w,2), sep = " = "), collapse = "; "),
       "sum_sample_size_c" = sum(x$data$n_samples_c),
       "sum_sample_size_t" = sum(x$data$n_samples_t),
       "sum_indiv_size_c" = sum(x$data$n_c),
       "sum_indiv_size_t" = sum(x$data$n_t))
     
})
table_synthesis <- do.call("rbind", table_synthesis)
table_synthesis <- data.frame("ID" = row.names(table_synthesis), table_synthesis,
                              stringsAsFactors = FALSE)

table_egger <- lapply(names(resList), function(nm, l) {
     x <- l[[nm]]
     if(nrow(x$data)>2) {
          egger <- regtest(x)
          data.frame("ID" = nm, "Std_Eff" = c("intercept", "slope_bias"), "coefficient" = egger$fit$beta,
                     "std_error" = egger$fit$se, "zvalue" = egger$fit$zval, 
                     "pvalue" = egger$fit$pval, "CI_low" = egger$fit$ci.lb,
                     "CI_high" = egger$fit$ci.ub)
     } else {
          data.frame("ID" = nm, "Std_Eff" = NA, "coefficient" = NA,
                     "std_error" = NA, "zvalue" = NA, 
                     "pvalue" = NA, "CI_low" = NA,
                     "CI_high" = NA)
     }
}, l = resList)
table_egger <- do.call("rbind", table_egger)

table_jackknife <- lapply(names(resList), function(nm, l) {
     x <- l[[nm]]
     jackk <- data.frame(leave1out(x))
     data.frame("ID" = nm, "study" = row.names(jackk), jackk)
}, l = resList)
table_jackknife <- do.call("rbind", table_jackknife)

forestPlots <- lapply(resList, function(x) {
     sum_table <- x$data[,c("study", "n_c", "n_t", "n_samples_c", "n_samples_t")]
     sum_table$weights <- round(weights(x), 3)
     sum_table <- rbind(sum_table, 
                        "All studies" = c("study" = "All studies",
                                      "n_c" = sum(sum_table$n_c),
                                      "n_t" = sum(sum_table$n_t),
                                      "n_samples_c" = sum(sum_table$n_samples_c),
                                      "n_samples_t" = sum(sum_table$n_samples_t),
                                      "weights" = round(sum(sum_table$weights))))
     sum_table$"n_c/n_t" <- paste(sum_table$n_c, sum_table$n_t, sep = "/")
     sum_table$"n_samples_c/n_samples_t" <- paste(sum_table$n_samples_c, sum_table$n_samples_t, sep = "/")
     sum_table <- sum_table[,c(1,7,8,6)]
     x_lim_vec <- c((min(x$data$ci_low) - 1), (max(x$data$ci_high) + 1))
     vfplot <- viz_forest(x = x, type = "standard", variant = "classic",
                study_labels = x$data$study,
                summary_label = "Summary RE Model", xlab = "SMD",
                study_table = sum_table[1:(nrow(sum_table)-1),],
                summary_table = sum_table[nrow(sum_table),,drop=FALSE],
                table_headers = c("Study", "Individuals (HC/AD)",
                                  "Samples (HC/AD)", "Weights (%)"), 
                text_size = 4,
                annotate_CI = TRUE, col = "gray",
                table_layout = matrix(c(2, 2, 1, 1, 3), nrow = 1),
                x_limit = x_lim_vec)
     return(vfplot)
})

# funnelPlots <- lapply(resList, function(x) {
#      viz_funnel(x = x, egger = TRUE, trim_and_fill = TRUE)
# })

#############################
# Visualization
#############################
# ndir <- paste("results_giovanna3", sep = "/")
# if(!dir.exists(ndir)) {
#      dir.create(ndir)   
# }
ndir <- getwd()
lapply(names(resList), function(table, resfit, fplot) {
     resfit <- resfit[[table]]
     fname <- paste(ndir, "/", table, ".pdf", sep = "")
     title_plot <- paste(table, " - pvalue=", round(resfit$pval,3), sep = "")
     pdf(file = fname, width = 12, height = 6)
     plot(forestPlots[[table]])
     funnel(resfit, main=title_plot)
     dev.off()
}, resfit = resList, fplot = forestPlots)
fname <- paste(ndir, "table_indiv.txt", sep = "/")
write.table(x = table_indiv, file = fname, sep = "\t", row.names = FALSE)
fname <- paste(ndir, "table_synthesis.txt", sep = "/")
write.table(x = table_synthesis, file = fname, sep = "\t", row.names = FALSE)
fname <- paste(ndir, "table_egger.txt", sep = "/")
write.table(x = table_egger, file = fname, sep = "\t", row.names = FALSE)
fname <- paste(ndir, "table_jackknife.txt", sep = "/")
write.table(x = table_jackknife, file = fname, sep = "\t", row.names = FALSE)
