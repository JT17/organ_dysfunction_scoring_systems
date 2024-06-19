rm(list = ls())
library(MASS)
library(effects)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(car)
library(kableExtra)
library(knitr)
library(table1)
library(gridExtra)
library(glmnet)
library(gtsummary)
setwd("/Users/JonathanTiao/Documents/CummingsResearch/")

#Helper Functions
fit_ordinal_regression <- function(response_col, predictor_cols, data) {
  formula_str <- as.formula(paste(response_col, "~", paste(predictor_cols, collapse = "+")))
  model <- MASS::polr(formula = formula_str, data = data, Hess=T)
  return(model)
}

fit_mv_ordinal_regression <- function(response_col, var_col, covariates, data){
  formula_str <- as.formula(paste(response_col, "~",var_col, " + ", paste(covariates,collapse="+" )))
  model <- polr(formula=formula_str, data=data, Hess=T)
  return (model)
}

calculate_ordinal_pval <- function(input_model){
  ctable <- coef(summary(input_model))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  ctable <- cbind(ctable, "p value" = p)
  return (ctable)
}

calculate_ordinal_ci <- function(input_model){
  ci <- exp(confint(input_model))
  or <- exp(coef(input_model))
  return(cbind(OR=or, rbind(ci)))
}

create_p_ci_df <- function(pvals_list, ci_list){
  p_list <- lapply(pvals_list,function(mat){
    as.data.frame(mat[1, ,drop=FALSE]) %>%
      tibble::rownames_to_column(var="RowName")
  })
  p_df <- bind_rows(p_list)
  ci_list <- lapply(ci_list, function(mat){
    as.data.frame(mat[1, , drop=FALSE]) %>% 
      tibble::rownames_to_column(var="RowName")
  })
  ci_df <- bind_rows(ci_list)
  merged_res <- merge(p_df, ci_df)
  return(merged_res)
}

create_combined_kable_table <- function(uv_input_table, 
                                        uv_input_cols, 
                                        mv_input_table, 
                                        mv_input_cols,
                                        title){
  output <- uv_input_table[,uv_input_cols]
  
  output$mv_p_val <- mv_input_table$`p value`
  output$mv_OR <- mv_input_table$OR
  output$mv_2.5 <- mv_input_table$`2.5 %`
  output$mv_97.5 <- mv_input_table$`97.5 %`
  is.num <- sapply(output, is.numeric)
  
  output[is.num] <- lapply(output[is.num], round,2)
  output$`p value`[output$`p value` < "0.001"] <- "<0.001"
  output$mv_p_val[output$mv_p_val < "0.001"] <- "<0.001" 
  output$formatted_p <- cell_spec(output$`p value`, 
                                  color = ifelse(output$mv_p_val < 0.05, "red", "black"))
  output$formatted_or <- cell_spec(output$OR, 
                                   color = ifelse(output$mv_p_val < 0.05, "red", "black"))
  calculated_ci <- paste("(", output$`2.5 %`, ",", output$`97.5 %`, ")")
  calculated_mv_ci <- paste("(", output$mv_2.5, ",", output$mv_97.5, ")")
  output$formatted_ci <- cell_spec(calculated_ci, color = ifelse(output$mv_p_val < 0.05, "red", "black"))
  output$formatted_mv_p <- cell_spec(output$mv_p_val,
                                     color=ifelse(output$mv_p_val < 0.05, "red", "black"))
  output$formatted_mv_OR <- cell_spec(output$mv_OR,
                                     color=ifelse(output$mv_p_val < 0.05, "red", "black"))
  output$formatted_mv_ci <- cell_spec(calculated_mv_ci, color = ifelse(output$mv_p_val < 0.05, "red", "black"))
  formatted_output <- output[,c("RowName", "formatted_or","formatted_ci", "formatted_p", "formatted_mv_OR", "formatted_mv_ci", "formatted_mv_p")]
  colnames(formatted_output) <- c("Soluble Mediator", "OR",  "95% CI", "p","Adjusted OR", "Adjusted 95% CI", "Adjusted p")
  formatted_output <- formatted_output %>% 
    mutate(`Soluble Mediator` = str_replace(`Soluble Mediator`, "log", ""))
  to_return <- kable(formatted_output,
                     row.names = FALSE,
                     format = "html", 
                     escape=FALSE) %>% 
    kable_styling(bootstrap_options = c("striped", "hover", full_width=F)) %>%
    footnote(general = "", general_title=title)
  return (to_return)
}

##end helper functions

sepsis_dt <-fread("/Users/JonathanTiao/Downloads/RESERVE-U-1 Cleaned Master Dataset_for_analyses_nonpreg_adult_n_301_biomarkers_trclust_medclust_csv.csv")
sepsis_dt$qsofa_factor <- factor(sepsis_dt$qsofa_score)
sepsis_dt$uva_factor <- factor(sepsis_dt$UVA_score)
sepsis_dt$mews_factor <- factor(sepsis_dt$MEWS_score)
adjustment_vars <- c("age", "sex", "hivrdtresult", "malariardtresult", "illnessduration")

#initialize empty lists
biomarker_names <- list(length=14)

#define strings, fill names list
logifn <- c('logifn')
biomarker_names[1] <- logifn
label(sepsis_dt$logifn) <- "log(ifn)"
logil6 <- c('logil6')
biomarker_names[2] <- logil6
logil8 <- c('logil8')
label(sepsis_dt$logil6) <- "log(IL-6)"
biomarker_names[3] <- logil8
label(sepsis_dt$logil8) <- "log(il8)"
logil10 <- c('logil10')
biomarker_names[4] <- logil10
label(sepsis_dt$logil10) <- "log(il10)"
logip10 <- c('logip10')
biomarker_names[5] <- logip10
label(sepsis_dt$logip10) <- "log(ip10)"
logmip1a <- c('logmip1a')
label(sepsis_dt$logmip1a) <- "log(mip1a)"
biomarker_names[6] <- logmip1a
logmip1b <- c('logmip1b')
label(sepsis_dt$logmip1b) <- "log(mip1b)"
biomarker_names[7] <- logmip1b
label(sepsis_dt$logtnf) <- "log(tnf)"
logtnf <- c('logtnf')
biomarker_names[8] <- logtnf
logang1 <- c('logang1')
biomarker_names[9] <- logang1
label(sepsis_dt$logang1) <- "log(ang1)"
logang2 <- c('logang2')
label(sepsis_dt$logang2) <- "log(ang-2)"
biomarker_names[10] <- logang2
logsil2ra <- c('logsil2ra')
label(sepsis_dt$logsil2ra) <- "log(sIL-2Ra)"
biomarker_names[11] <- logsil2ra
logstnfr1 <- c('logstnfr1')
label(sepsis_dt$logstnfr1) <- "log(stnfr1)"
biomarker_names[12] <- logstnfr1
logmif <- c('logmif')
label(sepsis_dt$logmif) <- "log(mif)"
biomarker_names[13] <- logmif
logpai1 <- c('logpai1')
biomarker_names[14] <- logpai1
label(sepsis_dt$logpai1) <- "log(pai1)"
label(sepsis_dt$qsofa_score) <- "qSOFA Score"
label(sepsis_dt$UVA_score) <- "UVA Score"
label(sepsis_dt$MEWS_score) <- "MEWS Score"

#drop any rows with missing values
sepsis_dt <- sepsis_dt[!(is.na(sepsis_dt$logifn))]
sepsis_dt <- sepsis_dt[!(is.na(sepsis_dt$logmif))]

#create univariate qsofa ordinal regressions
biomarker_ordinal_regs <- list(length=14)
biomarker_qsofa_pvals <- list(length=14)
biomarker_qsofa_cis <- list(length=14)

mv_biomarkers_qsofa <- list(length=14)
mv_biomarkers_pvals <- list(length=14)
mv_biomarkers_ci <- list(length=14)
for(i in 1:length(biomarker_names)){
  #univariate qsofa models
  biomarker_ordinal_regs[[i]] <- fit_ordinal_regression('qsofa_factor', biomarker_names[[i]], sepsis_dt)
  biomarker_qsofa_pvals[[i]] <- calculate_ordinal_pval(biomarker_ordinal_regs[[i]])
  biomarker_qsofa_cis[[i]] <- calculate_ordinal_ci(biomarker_ordinal_regs[[i]])
  
  #multivariate qsofa models
  mv_biomarkers_qsofa[[i]] <- fit_mv_ordinal_regression('qsofa_factor', 
                                                        var_col = biomarker_names[[i]], 
                                                        covariates = adjustment_vars, 
                                                        data=sepsis_dt)
  mv_biomarkers_pvals[[i]] <- calculate_ordinal_pval(mv_biomarkers_qsofa[[i]])
  mv_biomarkers_ci[[i]] <- calculate_ordinal_ci(mv_biomarkers_qsofa[[i]])
}
uv_ord_qsofa_results <- create_p_ci_df(pvals_list = biomarker_qsofa_pvals, ci_list = biomarker_qsofa_cis)
mv_ord_qsofa_results <- create_p_ci_df(pvals_list = mv_biomarkers_pvals, ci_list = mv_biomarkers_ci)

combined_qsofa <- create_combined_kable_table(uv_ord_qsofa_results, 
                                              c("RowName", "OR", "2.5 %", "97.5 %", "p value" ), 
                                              mv_ord_qsofa_results, 
                                              c("p value", "OR"),
                                              "Association of qSOFA Scores with Immune Mediators")

# uva ordinal regression
uv_ord_biomarkers_uva <- list(length=14)
uv_ord_biomarkers_uva_pvals <- list(length=14)
uv_ord_biomarkers_uva_ci <- list(length=14)

mv_ord_biomarkers_uva <- list(length=14)
mv_ord_biomarkers_uva_pvals <- list(length=14)
mv_ord_biomarkers_uva_ci <- list(length=14)

#remove hiv from UVA
uva_adjustment_vars <- c("age", "sex", "malariardtresult", "illnessduration")

for (i in 1:length(biomarker_names)){
  #univariate uva models
  uv_ord_biomarkers_uva[[i]] <- fit_ordinal_regression('uva_factor', biomarker_names[[i]], data=sepsis_dt)
  uv_ord_biomarkers_uva_pvals[[i]] <- calculate_ordinal_pval(uv_ord_biomarkers_uva[[i]])
  uv_ord_biomarkers_uva_ci[[i]] <- calculate_ordinal_ci(uv_ord_biomarkers_uva[[i]])
  
  #multivariate uva models
  mv_ord_biomarkers_uva[[i]] <- fit_mv_ordinal_regression('uva_factor', 
                                                          var_col = biomarker_names[[i]], 
                                                          covariates = uva_adjustment_vars, 
                                                          data=sepsis_dt)
  mv_ord_biomarkers_uva_pvals[[i]] <- calculate_ordinal_pval(mv_ord_biomarkers_uva[[i]])
  mv_ord_biomarkers_uva_ci[[i]] <- calculate_ordinal_ci(mv_ord_biomarkers_uva[[i]])
}

#create tables to display uva results
uv_ord_uva_results <- create_p_ci_df(pvals_list = uv_ord_biomarkers_uva_pvals, ci_list = uv_ord_biomarkers_uva_ci)
mv_ord_uva_results <- create_p_ci_df(pvals_list = mv_ord_biomarkers_uva_pvals, ci_list = mv_ord_biomarkers_uva_ci)

combined_uva <- create_combined_kable_table(uv_ord_uva_results, 
                                            c("RowName", "OR", "2.5 %", "97.5 %", "p value" ), 
                                            mv_ord_uva_results, 
                                            c("p value"),
                                            "Association of UVA Scores with Immune Mediators")

# mews ordinal regression
uv_ord_biomarkers_mews <- list(length=14)
uv_ord_biomarkers_mews_pvals <- list(length=14)
uv_ord_biomarkers_mews_ci <- list(length=14)

mv_ord_biomarkers_mews <- list(length=14)
mv_ord_biomarkers_mews_pvals <- list(length=14)
mv_ord_biomarkers_mews_ci <- list(length=14)

for (i in 1:length(biomarker_names)){
  uv_ord_biomarkers_mews[[i]] <- fit_ordinal_regression('mews_factor', biomarker_names[[i]], data=sepsis_dt)
  uv_ord_biomarkers_mews_pvals[[i]] <- calculate_ordinal_pval(uv_ord_biomarkers_mews[[i]])
  uv_ord_biomarkers_mews_ci[[i]] <- calculate_ordinal_ci(uv_ord_biomarkers_mews[[i]])
  
  mv_ord_biomarkers_mews[[i]] <- fit_mv_ordinal_regression('mews_factor', 
                                                           var_col = biomarker_names[[i]], 
                                                           covariates = adjustment_vars, 
                                                           data=sepsis_dt)
  mv_ord_biomarkers_mews_pvals[[i]] <- calculate_ordinal_pval(mv_ord_biomarkers_mews[[i]])
  mv_ord_biomarkers_mews_ci[[i]] <- calculate_ordinal_ci(mv_ord_biomarkers_mews[[i]])
}
uv_ord_mews_results <- create_p_ci_df(pvals_list = uv_ord_biomarkers_mews_pvals, ci_list = uv_ord_biomarkers_mews_ci)
mv_ord_mews_results <- create_p_ci_df(pvals_list = mv_ord_biomarkers_mews_pvals, ci_list = mv_ord_biomarkers_mews_ci)
combined_mews <- create_combined_kable_table(uv_ord_mews_results, 
                                             c("RowName", "OR", "2.5 %", "97.5 %", "p value" ), 
                                             mv_ord_mews_results, 
                                             c("p value"),
                                             "Association of MEWS Scores with Immune Mediators")

#lasso modeling
unlisted_bm_nm <- unlist(biomarker_names)
sepsis_lasso_input <- c(unlisted_bm_nm, "qsofa_binary", "uva_binary", "mews_binary")
sepsis_dt$sex <- factor(sepsis_dt$sex, levels=c(0,1), labels=c("Female", "Male"))
sepsis_dt$hiv_factor <- factor(sepsis_dt$hivrdtresult, levels=c(0,1), labels=c("Absent", "Present"))
sepsis_dt$malaria_factor <- factor(sepsis_dt$malariardtresult, levels=c(0,1), labels=c("Absent", "Present"))
sepsis_dt$qsofa_binary <- as.factor(as.numeric(sepsis_dt$qsofa_score >= 2))
sepsis_dt$uva_binary <- as.factor(as.numeric(sepsis_dt$UVA_score >= 4))
sepsis_dt$mews_binary <- as.factor(as.numeric(sepsis_dt$MEWS_score >= 4))

#select for just biomarkers
sepsis_lasso_dt <- sepsis_dt[,..sepsis_lasso_input]

#drop any missing biomarkers
sepsis_lasso_dt <- sepsis_lasso_dt[complete.cases(sepsis_lasso_dt),]
biomarker_input <- data.matrix(sepsis_lasso_dt[,..unlisted_bm_nm])

#LASSO logistic regression using qsofa >=2 as outcome
#first step - identify optimal lambda
qsofa_lasso <- cv.glmnet(biomarker_input, 
                         sepsis_lasso_dt$qsofa_binary, 
                         type.measure = "mse", 
                         alpha = 1, 
                         nfolds=10, 
                         family = "binomial")

qsofa_best_lambda <- qsofa_lasso$lambda.min

#re-train to identify coefficients associated with optimal lambda
best_qsofa_model <- glmnet(biomarker_input, 
                           sepsis_lasso_dt$qsofa_binary,
                           type.measure = "mse",
                           alpha=1,
                           nfolds=10,
                           family="binomial",
                           lambda = qsofa_best_lambda)
coef(best_qsofa_model)

#MEWS LASSO
mews_lasso <- cv.glmnet(biomarker_input, 
                        sepsis_lasso_dt$mews_binary, 
                        type.measure = "mse", 
                        alpha = 1, 
                        nfolds=10, 
                        family = "binomial")

mews_best_lambda <- mews_lasso$lambda.min

#re-train to identify coefficients associated with optimal lambda
best_mews_model <- glmnet(biomarker_input, sepsis_lasso_dt$mews_binary,
                          type.measure = "mse",
                          alpha=1,
                          nfolds=10,
                          family="binomial",
                          lambda = mews_best_lambda)
coef(best_mews_model)

#UVA LASSO
uva_lasso <- cv.glmnet(biomarker_input, 
                       sepsis_lasso_dt$uva_binary, 
                       type.measure = "mse", 
                       alpha = 1, 
                       nfolds=10, 
                       family = "binomial")

uva_best_lambda <- uva_lasso$lambda.min

#re-train to identify coefficients associated with optimal lambda
best_uva_model <- glmnet(biomarker_input, sepsis_lasso_dt$uva_binary,
                         type.measure = "mse",
                         alpha=1,
                         nfolds=10,
                         family="binomial",
                         lambda = uva_best_lambda)
coef(best_uva_model)


combined_coefs <- cbind(as.matrix(coef(best_qsofa_model)), as.matrix(coef(best_mews_model)), as.matrix(coef(best_uva_model)))
combined_coefs_df <- as.data.frame(combined_coefs)
colnames(combined_coefs_df) <- c("qSOFA Lasso Coefficients", "MEWS Lasso Coefficients", "UVA Lasso Coefficients")
combined_coefs_df <- combined_coefs_df[-1,]
combined_coefs_df <- round(combined_coefs_df, 2)
combined_coefs_df <- combined_coefs_df %>% mutate_all(~ifelse(.==0, "---", .))

lasso_coefs_kable <- kable(combined_coefs_df, "html", escape=FALSE, 
                           caption="Lasso Coefficients of Logistic Regression Comparing Biomarkers and qSOFA/MEWS/UVA") %>% 
  kable_styling()
lasso_coefs_kable

