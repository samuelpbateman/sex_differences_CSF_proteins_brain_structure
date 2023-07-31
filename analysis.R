# install GWASTools of BiocManager 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GWASTools")
# install qvalue of BiocManager 
BiocManager::install("qvalue")
# install sva of BiocManager 
BiocManager::install("sva")

# install meff for multiple testing corrections ====
meff = function(x, eigen = FALSE, correl = FALSE) {
  # if the input x is eigenvalues
  if(eigen) {
    if(!class(x) %in% c("numeric", "integer")) {
      stop("Eigenvalues are not numeric or integer.")
    }
    k = length(x)
    evs = x
    # if the input x is a correlation matrix       
  } else if(correl) {
    # number of variables
    k = ncol(x)
    # dimension checks
    if(!isSymmetric(x)) {
      stop("Correlation matrix is not symmetric.")
    } 
    # convert the correlation matrix to positive-definite
    require(Matrix)
    x = as.matrix(nearPD(x)$mat)
    # get eigenvalues and absolute eigenvalues of correlation matrix
    evs = eigen(x)$values
    # otherwise the input x is assumed to be a dataframe with variables in columns and observations in rows
  } else {
    # number of variables
    k = ncol(x)
    # get correlation matrix
    x = cor(x, use = "pairwise.complete.obs")
    # convert the correlation matrix to positive-definite.
    require(Matrix)
    x = as.matrix(nearPD(x)$mat)
    # get eigenvalues and absolute eigenvalues of R matrix
    evs = eigen(x)$values
  }
  # effective number of tests (Nyholt, 2004)
  eff = 1 + (k - 1) * (1 - var(evs) / k)
  return(eff)
}
meff(d[37:439])
  # 320.6846 for 3SD winsorisation
  # 320.3376 for 1.5SD winsorisation

# load packages and set seed ====
library(skimr) # for data inspection
library(tidyverse) # optional; ease of life 
library(summarytools) # for data inspection 
library(car) # for qqPlots and other stats
library(ggsignif) # to add significance labels onto ggplots
library(modelsummary) # 
library(psych) # for general stats
library(GWASTools) # for qqPlots of p values 
library(lmtest) # for robust SEs
library(sandwich) # for robust SEs
library(lavaan) # library for SEM
library(effectsize) # for effect sizes from models
  options(es.use_symbols = TRUE) # use mathematical symbols
set.seed(1)

# import data and create df for only complete cases ====
d <- as.data.frame(read_csv('finaldata_winsorized.csv', na = c("NA", "", -4))) # 817 rows, 442 vars
d <- d %>% mutate(...1 = NULL, ...2 = NULL, ...4 = NULL) # 439 cols now 
# names(d) # 37:439 are proteins (357 - 439 are from RBM data)
  # drop_na(d) %>% dim() # this significantly reduces sample size; 223 rows
d_noNA_full <- drop_na(d) # 75 cases with full data (also removes MRI NAs)
d_noNA <- d[complete.cases(d[,c(37:439)]),] # 99 cases with full protein data but not necessarily MRI

# Running checks on variables ====
glimpse(d[,1:38]) # appropriately convert variables 
  d$PTGENDER <- factor(d$PTGENDER)
  d$DX <- factor(d$DX)
  d$SITE <- factor(d$SITE)
  d$APOE4 <- factor(d$APOE4)
  d$SMOK <- factor(d$SMOK, levels = c(0,1))
  
# Winsorise proteins at ±1.5SD
d[37:439] <- DescTools::Winsorize(d[39:441], maxval=1.5, minval=-1.5, na.rm=F)


# Subsetting data for demographics table to reflect regression analyses with na.exclude() 
dim(na.exclude(d)) # n = 223, 439 vars, less than the 286 from RQ1 analyses
length(na.exclude(d$NEGR1.VVVNFAPTIQEIK)) # 286 
dim(na.exclude(d[37:439])) # 286, good, this is used to create demographics table 
na.exclude(d[37:439]) %>% select("RID":"ORIGPROT") %>% dfSummary()
dfSummary(d[na.exclude(d[39:439])])

d %>% select("RID", "A1AT.AVLTIDEK":"von.Willebrand.Factor..vWF...ug.mL.") %>% 
  na.exclude() %>% dim() # 286 cases, good
demTable <- d %>% select("RID", "A1AT.AVLTIDEK":"von.Willebrand.Factor..vWF...ug.mL.") %>% 
  na.exclude() 
d2 <- d[1:36]
demTable2 <- merge(demTable, d2, by="RID", all=F)
dim(demTable2) # 286 cases, 439 vars

# get stats for demographics table
dfSummary(demTable2[405:436])


# check whether brain structures have 0s.
  d %>% filter(Ventricles < 1) %>% glimpse() # none
  d %>% filter(Ventricles_corrected < 0.001) %>% glimpse() # none
  
# investigate sites
skim(d$SITE) # 58 unique 
unique(d$SITE) # doesn't look right - 914 is one unique value? 
hist(d$SITE)


# calculate brain phenotypes by diagnosis
demTable2 %>% select(DX, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% filter(DX == "Dementia") %>% dfSummary()
demTable2 %>% select(DX, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% filter(DX == "MCI") %>% dfSummary()
demTable2 %>% select(DX, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% filter(DX == "CN") %>% dfSummary()

# calcualte brain phenotypes by sex
demTable2 %>% select(DX, PTGENDER, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% filter(PTGENDER == "Female") %>% dfSummary()
demTable2 %>% select(DX, PTGENDER, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% filter(PTGENDER == "Male") %>% dfSummary()

# calcualte brain phenotypes by sex and diagnosis
demTable2 %>% select(DX, PTGENDER, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% 
  filter(PTGENDER == "Male" & DX == "CN") %>% dfSummary()
demTable2 %>% select(DX, PTGENDER, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% 
  filter(PTGENDER == "Male" & DX == "MCI") %>% dfSummary()
demTable2 %>% select(DX, PTGENDER, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% 
  filter(PTGENDER == "Male" & DX == "Dementia") %>% dfSummary()
demTable2 %>% select(DX, PTGENDER, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% 
  filter(PTGENDER == "Female" & DX == "CN") %>% dfSummary()
demTable2 %>% select(DX, PTGENDER, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% 
  filter(PTGENDER == "Female" & DX == "MCI") %>% dfSummary()
demTable2 %>% select(DX, PTGENDER, ICV, WholeBrain, Ventricles, Hippocampus, Fusiform) %>% 
  filter(PTGENDER == "Female" & DX == "Dementia") %>% dfSummary()


# Sex difference in brain phenotypes and proteins ==== 
  ## Total brain volume ====
    # run multiple regression for all proteins against corrected whole brain and extract only sex interaction stats
    BV_proteins <- lapply(d[37:439], function(p) 
      round(summary(lm(WholeBrain_corrected ~ p * PTGENDER + AGE + DX, na.action=na.exclude, data=d))$coefficients[7,], 7))
    # convert output to a table
    BV_proteins_table <- do.call(rbind, BV_proteins) 
    write.csv(BV_proteins_table, "BV_proteins_table.csv", row.names = T) # all vars
    # subset where p < .05
      BV_proteins_df <- as.data.frame(BV_proteins_table) # convert to df
      (BV_proteins_sig <- BV_proteins_df %>% filter(.[,4] < 0.05)) # 
    # subset where p < FDR correction
      BV_proteins_fdr <- p.adjust(BV_proteins_table[,4], "fdr", n=length(BV_proteins_table[,4]))
      BV_proteins_fdr <- p.adjust(BV_proteins_table[,4], "fdr", n=320.6846)
      
      BV_proteins_fdr_df <- as.data.frame(BV_proteins_fdr) # convert to df
      (BV_proteins_fdr_sig <- BV_proteins_fdr_df %>% filter(BV_proteins_fdr < 0.05)) # 
    # subset where p < Bonferroni correction
      BV_proteins_bon <- p.adjust(BV_proteins_table[,4], "bonferroni", n=length(BV_proteins_table[,4])) # all are 1?? 
      BV_proteins_bon_df <- as.data.frame(BV_proteins_bon) # convert to df
      (BV_proteins_bon_sig <- BV_proteins_bon_df %>% filter(BV_proteins_bon < 0.05)) #
      # manual calcualtion to check
        BV_proteins_bon_sig_manual <- BV_proteins_df %>% filter(.[,4] < 0.00012) # 
        BV_proteins_bon_sig_manual <- BV_proteins_df %>% filter(.[,4] < (0.05/320.6846)) # 
    # make qqplot of p values
    par(mfrow=c(2,2))
    GWASTools::qqPlot(BV_proteins_table[,4], ci = T, main = "Total Brain") # underinflation present
    
    # nobs(model) = 283 included cases in example model
    
    # calculate robust standard errors for model coefficients
      # create new list of regression outputs (no rounded summary)
      BV_proteins_robustSE <- lapply(d[37:439], function(p) 
        lm(WholeBrain_corrected ~ p * PTGENDER + AGE + DX, na.action=na.exclude, data=d))
      # run all results through coeftest() 
      BV_proteins_robustSE <- lapply(BV_proteins_robustSE, function(x) 
        coeftest(x, vcov = vcovHC(x, type = "HC0"))) # run models 
      # extract only interaction coefficients
      BV_proteins_robustSE_int_coef <- lapply(BV_proteins_robustSE, function(x) 
        x[grepl("p:PTGENDER", rownames(x)),]) # extract only interaction terms
      BV_proteins_robustSE_int_coef_table <- do.call(rbind, BV_proteins_robustSE_int_coef) # convert to table
      BV_proteins_robustSE_int_coef_df <- as.data.frame(BV_proteins_robustSE_int_coef_table) # convert to df
      write.csv(BV_proteins_robustSE_int_coef_table, "BV_proteins_table_rSE_1.5.csv", row.names = T) # all vars
      # make qqplot of robust SE p values
      GWASTools::qqPlot(BV_proteins_robustSE_int_coef_table[,4], ci = T, main = "Total Brain Robust SE")
      # Bonferroni correction
        View(BV_proteins_robustSE_int_coef_df[order(BV_proteins_robustSE_int_coef_df[,4]), ])
      # FDR correction
        BV_proteins_fdr_rSE <- p.adjust(BV_proteins_robustSE_int_coef_table[,4], "fdr", 
                                        n=length(BV_proteins_robustSE_int_coef_table[,4]))
        BV_proteins_fdr_rSE_df <- as.data.frame(BV_proteins_fdr_rSE) # convert to df
        (BV_proteins_fdr_rSe_sig <- BV_proteins_fdr_rSE_df %>% filter(BV_proteins_fdr_rSE_df < 0.05))
        # checking automated results are same as manual
          m1 <- lm(WholeBrain_corrected ~ A1AT.AVLTIDEK * PTGENDER + AGE + DX, na.action=na.exclude, data=d)
          summary(m1)
          coeftest(m1, vcov = vcovHC(m1, type = "HC0")) # success 
        
      
  ## Ventricles ==== 
    VE_proteins <- lapply(d[37:439], function(p) 
      round(summary(lm(Ventricles_corrected ~ p * PTGENDER + AGE + DX, na.action=na.exclude, data=d))$coefficients[7,], 7))
    # convert output to a table
    VE_proteins_table <- do.call(rbind, VE_proteins) 
    write.csv(VE_proteins_table, "VE_proteins_table.csv", row.names = T) # all vars
      # subset where p < .05
        VE_proteins_df <- as.data.frame(VE_proteins_table) # convert to df
        (VE_proteins_sig <- VE_proteins_df %>% filter(.[,4] < 0.05))
      # subset where p < FDR correction
        VE_proteins_fdr <- p.adjust(VE_proteins_table[,4], "fdr", n=length(VE_proteins_table[,4]))
        VE_proteins_fdr_df <- as.data.frame(VE_proteins_fdr) # convert to df
        (VE_proteins_fdr_sig <- VE_proteins_fdr_df %>% filter(VE_proteins_fdr < 0.05)) 
      # subset where p < Bonferroni correction
        VE_proteins_bon <- p.adjust(VE_proteins_table[,4], "bonferroni", n=length(VE_proteins_table[,4]))
        VE_proteins_bon_df <- as.data.frame(VE_proteins_bon) # convert to df
        (VE_proteins_bon_sig <- VE_proteins_bon_df %>% filter(VE_proteins_bon < 0.05))
        # manual calcualtion to check
          VE_proteins_bon_sig_manual <- VE_proteins_df %>% filter(.[,4] < 0.00012)
        # SCG1.HLEEPGETQNAFLNER nearing significance at FDR level. No BF significances. 
    # make qqplot of p values
    GWASTools::qqPlot(VE_proteins_table[,4], ci = T, main = "Ventricles")
    
    # calculate robust standard errors for model coefficients
      # create new list of regression outputs (no rounded summary)
      VE_proteins_robustSE <- lapply(d[37:439], function(p) 
        lm(Ventricles_corrected ~ p * PTGENDER + AGE + DX, na.action=na.exclude, data=d))
      # run all results through coeftest() 
      VE_proteins_robustSE <- lapply(VE_proteins_robustSE, function(x) 
        coeftest(x, vcov = vcovHC(x, type = "HC0")))
      # extract only interaction coefficients
      VE_proteins_robustSE_int_coef <- lapply(VE_proteins_robustSE, function(x) 
        x[grepl("p:PTGENDER", rownames(x)),]) # extract only interaction terms
      VE_proteins_robustSE_int_coef_table <- do.call(rbind, VE_proteins_robustSE_int_coef) # convert to table
      VE_proteins_robustSE_int_coef_df <- as.data.frame(VE_proteins_robustSE_int_coef_table) # convert to df
      write.csv(VE_proteins_robustSE_int_coef_df, "VE_proteins_table_rSE_1.5.csv", row.names = T) # all vars
      # make qqplot of robust SE p values
      GWASTools::qqPlot(VE_proteins_robustSE_int_coef_table[,4], ci = T, main = "Ventricles Robust SE")
        # Bonferroni correction
        View(VE_proteins_robustSE_int_coef_df[order(VE_proteins_robustSE_int_coef_df[,4]), ])
        # FDR correction
        VE_proteins_fdr_rSE <- p.adjust(VE_proteins_robustSE_int_coef_table[,4], "fdr", 
                                        n=length(VE_proteins_robustSE_int_coef_table[,4]))
        VE_proteins_fdr_rSE_df <- as.data.frame(VE_proteins_fdr_rSE) # convert to df
        (VE_proteins_fdr_rSE_sig <- VE_proteins_fdr_rSE_df %>% filter(VE_proteins_fdr_rSE_df < 0.05))
        # Insulin.like.Growth.Factor.Binding.Prote..ng.mL. is v close to significance at FDR level

          
  ## Hippocampus ==== 
    HI_proteins <- lapply(d[37:439], function(p) 
      round(summary(lm(Hippocampus_corrected ~ p * PTGENDER + AGE + DX, na.action=na.exclude, data=d))$coefficients[7,], 7))
    # convert output to a table
    HI_proteins_table <- do.call(rbind, HI_proteins) 
    write.csv(HI_proteins_table, "HI_proteins_table.csv", row.names = T) # all vars
      # subset where p < .05
        HI_proteins_df <- as.data.frame(HI_proteins_table) # convert to df
        (HI_proteins_sig <- HI_proteins_df %>% filter(.[,4] < 0.05))
      # subset where p < FDR correction
        HI_proteins_fdr <- p.adjust(HI_proteins_table[,4], "fdr", n=length(HI_proteins_table[,4]))
        HI_proteins_fdr_df <- as.data.frame(HI_proteins_fdr) # convert to df
        (HI_proteins_fdr_sig <- HI_proteins_fdr_df %>% filter(HI_proteins_fdr < 0.05)) 
      # subset where p < Bonferroni correction
        HI_proteins_bon <- p.adjust(HI_proteins_table[,4], "bonferroni", n=length(HI_proteins_table[,4]))
        HI_proteins_bon_df <- as.data.frame(HI_proteins_bon) # convert to df
        (HI_proteins_bon_sig <- HI_proteins_bon_df %>% filter(HI_proteins_bon < 0.05))
        # manual calcualtion to check
          HI_proteins_bon_sig_manual <- HI_proteins_df %>% filter(.[,4] < 0.00012)
          # none sig at BF or FDR level
    # make qqplot of p values
    GWASTools::qqPlot(HI_proteins_table[,4], ci = T, main = "Hippocampus")
    
    # calculate robust standard errors for model coefficients
      # create new list of regression outputs (no rounded summary)
      HI_proteins_robustSE <- lapply(d[37:439], function(p) 
        lm(Hippocampus_corrected ~ p * PTGENDER + AGE + DX, na.action=na.exclude, data=d))
      # run all results through coeftest() 
      HI_proteins_robustSE <- lapply(HI_proteins_robustSE, function(x) 
        coeftest(x, vcov = vcovHC(x, type = "HC0"))) 
      # extract only interaction coefficients
      HI_proteins_robustSE_int_coef <- lapply(HI_proteins_robustSE, function(x) 
        x[grepl("p:PTGENDER", rownames(x)),]) # extract only interaction terms
      HI_proteins_robustSE_int_coef_table <- do.call(rbind, HI_proteins_robustSE_int_coef) # convert to table
      HI_proteins_robustSE_int_coef_df <- as.data.frame(HI_proteins_robustSE_int_coef_table) # convert to df 
      write.csv(HI_proteins_robustSE_int_coef_df, "HI_proteins_table_rSE_1.5.csv", row.names = T) # all vars
      # make qqplot of robust SE p values
      GWASTools::qqPlot(HI_proteins_robustSE_int_coef_table[,4], ci = T, main = "Hippocampus Robust SE")
        # Bonferroni correction
        View(HI_proteins_robustSE_int_coef_df[order(HI_proteins_robustSE_int_coef_df[,4]), ])
          # none sig but ENOG.LGAEVYHTLK nearly 
        # FDR correction
        HI_proteins_fdr_rSE <- p.adjust(HI_proteins_robustSE_int_coef_table[,4], "fdr", 
                                        n=length(HI_proteins_robustSE_int_coef_table[,4]))
        HI_proteins_fdr_rSE_df <- as.data.frame(HI_proteins_fdr_rSE) # convert to df
        (HI_proteins_fdr_rSE_sig <- HI_proteins_fdr_rSE_df %>% filter(HI_proteins_fdr_rSE_df < 0.05))
        # T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. sig at FDR level 
    
    
  ## Fusiform Gyrus ====
    FG_proteins <- lapply(d[37:439], function(p) 
      round(summary(lm(Fusiform_corrected ~ p * PTGENDER + AGE + DX, na.action=na.exclude, data=d))$coefficients[7,], 7))
    # convert output to a table
    FG_proteins_table <- do.call(rbind, FG_proteins) 
    write.csv(FG_proteins_table, "FG_proteins_table.csv", row.names = T) # all vars
      # subset where p < .05
        FG_proteins_df <- as.data.frame(FG_proteins_table) # convert to df
        (FG_proteins_sig <- FG_proteins_df %>% filter(.[,4] < 0.05)) 
      # subset where p < FDR correction
        FG_proteins_fdr <- p.adjust(FG_proteins_table[,4], "fdr", n=length(FG_proteins_table[,4]))
        FG_proteins_fdr_df <- as.data.frame(FG_proteins_fdr) # convert to df
        (FG_proteins_fdr_sig <- FG_proteins_fdr_df %>% filter(FG_proteins_fdr < 0.05)) 
      # subset where p < Bonferroni correction
        FG_proteins_bon <- p.adjust(FG_proteins_table[,4], "bonferroni", n=length(FG_proteins_table[,4]))
        FG_proteins_bon_df <- as.data.frame(FG_proteins_bon) # convert to df
        (FG_proteins_bon_sig <- FG_proteins_bon_df %>% filter(FG_proteins_bon < 0.05)) 
        # manual calcualtion to check
          FG_proteins_bon_sig_manual <- FG_proteins_df %>% filter(.[,4] < 0.00012) 
    # make qqplot of log10-p values
    GWASTools::qqPlot(FG_proteins_table[,4], ci = T, main = "Fusiform Gyrus")
    
    # calculate robust standard errors for model coefficients
      # create new list of regression outputs (no rounded summary)
      FG_proteins_robustSE <- lapply(d[37:439], function(p) 
        lm(Fusiform_corrected ~ p * PTGENDER + AGE + DX, na.action=na.exclude, data=d))
      # run all results through coeftest() 
      FG_proteins_robustSE <- lapply(FG_proteins_robustSE, function(x) 
        coeftest(x, vcov = vcovHC(x, type = "HC0")))
      # extract only interaction coefficients
      FG_proteins_robustSE_int_coef <- lapply(FG_proteins_robustSE, function(x) 
        x[grepl("p:PTGENDER", rownames(x)),]) # extract only interaction terms
      FG_proteins_robustSE_int_coef_table <- do.call(rbind, FG_proteins_robustSE_int_coef) # convert to table
      FG_proteins_robustSE_int_coef_df <- as.data.frame(FG_proteins_robustSE_int_coef_table) # convert to df
      write.csv(FG_proteins_robustSE_int_coef_df, "FG_proteins_table_rSE_1.5.csv", row.names = T) # all vars
      # make qqplot of robust SE p values
      GWASTools::qqPlot(FG_proteins_robustSE_int_coef_table[,4], ci = T, main = "Fusiform Gyrus Robust SE")
        # Bonferroni correction 
        View(FG_proteins_robustSE_int_coef_df[order(FG_proteins_robustSE_int_coef_df[,4]), ])
          # ENOG.LGAEVYHTLK and SCG1.HLEEPGETQNAFLNER are significant 
        # FDR correction
        FG_proteins_fdr_rSE <- p.adjust(FG_proteins_robustSE_int_coef_table[,4], "fdr", 
                                        n=length(FG_proteins_robustSE_int_coef_table[,4]))
        FG_proteins_fdr_rSE_df <- as.data.frame(FG_proteins_fdr_rSE) # convert to df
        (FG_proteins_fdr_rSE_sig <- FG_proteins_fdr_rSE_df %>% filter(FG_proteins_fdr_rSE_df < 0.05))

          
  ### Regression diagnostics for 5 sig proteins ====
  na.exclude(d) %>% length() # n = 439
  
  # CNTNAP2
  WB_CNTNAP2 <- lm(WholeBrain_corrected ~ CNTP2.VDNAPDQQNSHPDLAQEEIR * PTGENDER + AGE + DX, na.action=na.exclude, data=d)
  confint(WB_CNTNAP2)
  par(mfrow=c(2,2))
  plot(WB_CNTNAP2)
    # resid vs fitted: roughly horizontal, fine
    # qq: fine
    # scale-location: roughly horizontal, fine
    # resid vs lev: no outliers of cook's d
  # multicollinearity 
  library(car)
  vif(WB_CNTNAP2) # protein main effect and interaction > 2 indicating possible mild correlation between predictor vars
  ## plot effect
    d %>% ggplot(., aes(x = CNTP2.VDNAPDQQNSHPDLAQEEIR, y = WholeBrain_corrected, colour=PTGENDER)) +
          geom_smooth(method = lm) + 
          geom_point(shape=d$PTGENDER, size = 2) + 
          labs(x = "Scaled and Centred CNTNAP2 Concentration", 
               y = "ICV-corrected Whole Brain Volume",
               color = "Sex") +
          theme_classic(base_size = 14)
          # plot shows significant skew from a few v low winsorized cases. Stronger p cor in F than M
  # effect size 
  eta_squared(WB_CNTNAP2) # 0.01=2 (0.00, 1)
    # sensitivity analysis to remove extreme values and re-plot
    d$CNTNAP2_2 <- DescTools::Winsorize(d$CNTP2.VDNAPDQQNSHPDLAQEEIR, maxval=1.5, minval=-1.5, na.rm=F)
    WB_CNTNAP2_2 <- lm(WholeBrain_corrected ~ CNTNAP2_2 * PTGENDER + AGE + DX, na.action=na.exclude, data=d)
    summary(WB_CNTNAP2_2) # t = -1.820  p = 0.070
    par(mfrow=c(2,2))
    plot(WB_CNTNAP2_2)
      # all fine still 
    # multicollinearity 
    library(car)
    vif(WB_CNTNAP2_2) # protein main effect and interaction reduced under 2 now 
    ## re-plot effect 
      d %>% ggplot(., aes(x = CNTNAP2_2, y = WholeBrain_corrected, colour=PTGENDER)) + 
            geom_smooth(method = lm) + 
            geom_point(shape=d$PTGENDER, size = 2) + 
            labs(x = "Scaled and Centred CNTNAP2 Concentration", 
                 y = "ICV-corrected Whole Brain Volume",
                 color = "Sex") +
            theme_classic(base_size = 14)
            # plot shows same pattern as before, stronger p correlation in F than M, perhaps slightly smaller difference 
      
  # CCL5 
  HV_CCL5 <- lm(Hippocampus_corrected ~ T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. * PTGENDER + AGE + DX, na.action=na.exclude, data=d)
  plot(HV_CCL5)
    # resid vs fitted: roughly horizontal, fine
    # qq: fine but a few cases towards top end, 60
    # scale-location: roughly horizontal, fine
    # resid vs lev: no outliers of c'd
  # multicollinearity
  vif(HV_CCL5) # all < 2, fine 
  ## plot effect
    d %>% ggplot(., aes(x = T.Cell.Specific.Protein.RANTES..RANTES...ng.mL., y = Hippocampus_corrected, colour=PTGENDER)) +
          geom_smooth(method = lm) + 
          geom_point(shape = d$PTGENDER, size=2) + 
          labs(x = "Scaled and Centred CCL5 Concentration", 
               y = "ICV-corrected Hippocampal Volume",
               color= "Sex") + 
          theme_classic(base_size = 14)
    # sensitivity analysis to remove extreme values and re-plot
    d$CCL5_2 <- DescTools::Winsorize(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL., maxval=1.5, minval=-1.5, na.rm=F)
    HV_CCL5_2 <- lm(Hippocampus_corrected ~ CCL5_2 * PTGENDER + AGE + DX, na.action=na.exclude, data=d)
      summary(HV_CCL5_2) # t = -2.426  p = 0.01602 *  
      par(mfrow=c(2,2))
      plot(HV_CCL5_2)
      # all okay but case 60 skews qqplot at high end 
    # multicollinearity 
      library(car)
      vif(HV_CCL5_2) # all < 2 
    # re-plot effect 
    d %>% ggplot(., aes(x = CCL5_2, y = Hippocampus_corrected, colour=PTGENDER)) +
      geom_smooth(method = lm) + 
      geom_point(shape = d$PTGENDER, size= 2) + 
      labs(x = "Scaled and Centred CCL5 Concentration", 
           y = "ICV-corrected Hippocampal  Volume",
           color = "Sex") + 
      theme_classic(base_size = 14)
    # relationships are much closer but F still show a greater 
    
  # CHGB
    FG_CHGB <- lm(Fusiform_corrected ~ SCG1.HLEEPGETQNAFLNER * PTGENDER + AGE + DX, na.action=na.exclude, data=d)
    plot(FG_CHGB)
    # resid vs fitted: roughly horizontal, fine
    # qq: fine 
    # scale-location: roughly horizontal, fine
    # resid vs lev: no outliers of c'd
  # multicollinearity
    vif(FG_CHGB) # all < 2
  ## plot effect
    d %>% ggplot(., aes(x = SCG1.HLEEPGETQNAFLNER, y = Fusiform_corrected, colour=PTGENDER)) +
          geom_smooth(method = lm) + 
          geom_point(shape = d$PTGENDER, size= 2) + 
          labs(x = "Scaled and Centred CHGB Concentration", 
               y = "ICV-corrected Fusiform Gyrus Volume",
               color = "Sex") + 
          theme_classic(base_size = 14)
          # plot shows a few low winsorized cases from F. Stronger p cor in M than F.
    # sensitivity analysis to remove extreme values and re-plot
    d$CHGB_2 <- DescTools::Winsorize(d$SCG1.HLEEPGETQNAFLNER, maxval=1.5, minval=-1.5, na.rm=F)
    FG_CHGB_2 <- lm(Fusiform_corrected ~ CHGB_2 * PTGENDER + AGE + DX, na.action=na.exclude, data=d)
      summary(FG_CHGB_2) # t = 2.760, p = 0.00626
      par(mfrow=c(2,2))
      plot(FG_CHGB_2)
        # all fine still 
    # multicollinearity 
    library(car)
    vif(FG_CHGB_2) # all < 2 
    # re-plot effect 
      d %>% ggplot(., aes(x = CHGB_2, y = Fusiform_corrected, colour=PTGENDER)) +
        geom_smooth(method = lm) + 
        geom_point(shape = d$PTGENDER, size= 2) + 
        labs(x = "Scaled and Centred CHGB Concentration", 
             y = "ICV-corrected Fusiform Gyrus Volume",
             color = "Sex") + 
        theme_classic(base_size = 14)
        # relationships are much closer but F still show a greater 
      
  # CLSTN3
    FG_CLSTN3 <- lm(Fusiform_corrected ~ CSTN3.ESLLLDTTSLQQR * PTGENDER + AGE + DX, na.action=na.exclude, data=d)
    plot(FG_CLSTN3)
    # resid vs fitted: roughly horizontal, fine
    # qq: fine
    # scale-location: roughly horizontal, slight curve 
    # resid vs lev: no outliers of c'd
  # multicollinearity
    vif(FG_CLSTN3) # all < 2 
  ## plot effect
    d %>% ggplot(., aes(x = CSTN3.ESLLLDTTSLQQR, y = Fusiform_corrected, colour=PTGENDER)) +
          geom_smooth(method = lm) + 
          geom_point(shape = d$PTGENDER, size = 2) + 
          labs(color = "Sex",
               x = "Scaled and Centred CLSTN3 Concentration", 
               y = "ICV-corrected Fusiform Gyrus Volume") + 
          theme_classic(base_size = 14)
    
    
  # NEGR1
    FG_NEGR1 <- lm(Fusiform_corrected ~ NEGR1.SSIIFAGGDK * PTGENDER + AGE + DX, na.action=na.exclude, data=d)
    plot(FG_NEGR1)
    # resid vs fitted: roughly horizontal, fine
    # qq: fine
    # scale-location: roughly horizontal, slight curve at end 
    # resid vs lev: no outliers of c'd
  # multicollinearity
    vif(FG_NEGR1) # none > 2 
  ## plot effect
    d %>% ggplot(., aes(x = NEGR1.SSIIFAGGDK, y = Fusiform_corrected, colour=PTGENDER)) +
      geom_smooth(method = lm) + 
      geom_point(shape = d$PTGENDER, size = 2) + 
      labs(color = "Sex",
           x = "Scaled and Centred NEGR1 Concentration", 
           y = "ICV-corrected Fusiform Gyrus Volume") + 
      theme_classic(base_size = 14) \
  # effect size 
    eta_squared(FG_NEGR1) # 0.04 (LCI: 0.01 - UCI: 1.00)
    ## plot by DX
    d %>% ggplot(., aes(x = NEGR1.SSIIFAGGDK, y= Fusiform_corrected, colour = PTGENDER)) + 
      geom_smooth(method = lm) + 
      geom_point(shape=d$PTGENDER, size = 2) + 
      labs(x = "Scaled and Centred NEGR1 Concentration", 
           y = "ICV-corrected Fusiform Gyrus Volume",
           color = "Sex") +
      facet_wrap(~DX, nrow=3) + 
      theme_classic(base_size = 14)
    # different by diagnosis
    
    
  # code for removing specific cases 
  d2 <- d[!(d$...1 == 221), ]

  ### Sex-disaggregated analyses for phenotypes ====
  # Do the 3 proteins of interest sig predict their respective phenotype in females
  females <- subset(d, d$PTGENDER == "Female")
  f <- list()
  f$CCL5 <- lm(Hippocampus_corrected ~ T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. + AGE + DX, na.action=na.exclude, data=females)
  f$CLSTN3 <- lm(Fusiform_corrected ~ CSTN3.ESLLLDTTSLQQR + AGE + DX, na.action=na.exclude, data=females)
  f$NEGR1 <- lm(Fusiform_corrected ~ NEGR1.SSIIFAGGDK + AGE + DX, na.action=na.exclude, data=females)
  
  f_results <- do.call(cbind, lapply(f, function(z) 
    summary(z)$coefficients[2,]))
  f_results <- t(f_results)
  write.csv(f_results, "sex_disaggregated_results_females.csv") # save 
  
  # Do the 3 proteins of interest sig predict their respective phenotypes in males
  males <- subset(d, d$PTGENDER == "Male")
  m <- list()
  m$CCL5 <- lm(Hippocampus_corrected ~ T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. + AGE + DX, na.action=na.exclude, data=males)
  m$CLSTN3 <- lm(Fusiform_corrected ~ CSTN3.ESLLLDTTSLQQR + AGE + DX, na.action=na.exclude, data=males)
  m$NEGR1 <- lm(Fusiform_corrected ~ NEGR1.SSIIFAGGDK + AGE + DX, na.action=na.exclude, data=males)
  
  m_results <- do.call(cbind, lapply(m, function(z) 
    summary(z)$coefficients[2,]))
  m_results <- t(m_results)
  write.csv(m_results, "sex_disaggregated_results_males.csv") # save 
  # 
  
  #### Protein descriptive stats ====
  # CNTNP2: CNTP2.VDNAPDQQNSHPDLAQEEIR
  describeBy(d$CNTP2.VDNAPDQQNSHPDLAQEEIR, group = d$PTGENDER)
    hist(d$CNTP2.VDNAPDQQNSHPDLAQEEIR) # extreme values on low end 
    hist(d$CNTP2.VDNAPDQQNSHPDLAQEEIR[d$PTGENDER == "Female"]) # similar to above
    hist(d$CNTP2.VDNAPDQQNSHPDLAQEEIR[d$PTGENDER == "Male"]) # similar to above 
    describeBy(d$CNTP2.VDNAPDQQNSHPDLAQEEIR, group = d$DX)
    hist(d$CNTP2.VDNAPDQQNSHPDLAQEEIR[d$DX == "CN"]) # outliers on low end
    hist(d$CNTP2.VDNAPDQQNSHPDLAQEEIR[d$DX == "MCI"]) # similar 
    hist(d$CNTP2.VDNAPDQQNSHPDLAQEEIR[d$DX == "Dementia"]) # similar 
    
  # CCL5: T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.
    describeBy(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL., group = d$PTGENDER)
    hist(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) # extreme values at high end 
    hist(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.[d$PTGENDER == "Female"]) # similar
    hist(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.[d$PTGENDER == "Male"]) # similar
    describeBy(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL., group = d$DX)
    hist(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.[d$DX == "CN"]) # non-normal, right skewed and outliers on high end
    hist(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.[d$DX == "MCI"]) # similar 
    hist(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.[d$DX == "Dementia"]) # similar 
  
  # CHGB: SCG1.HLEEPGETQNAFLNER
    describeBy(d$SCG1.HLEEPGETQNAFLNER, group = d$PTGENDER)
    hist(d$SCG1.HLEEPGETQNAFLNER) # extreme values at low end 
    hist(d$SCG1.HLEEPGETQNAFLNER[d$PTGENDER == "Female"]) # similar
    hist(d$SCG1.HLEEPGETQNAFLNER[d$PTGENDER == "Male"]) # normal
    describeBy(d$SCG1.HLEEPGETQNAFLNER, group = d$DX)
    hist(d$SCG1.HLEEPGETQNAFLNER[d$DX == "CN"]) # outliers on low end
    hist(d$SCG1.HLEEPGETQNAFLNER[d$DX == "MCI"]) # normal 
    hist(d$SCG1.HLEEPGETQNAFLNER[d$DX == "Dementia"]) # slight right skew

  # CLSTN3: CSTN3.ESLLLDTTSLQQR
    describeBy(d$CSTN3.ESLLLDTTSLQQR, group = d$PTGENDER)
    hist(d$CSTN3.ESLLLDTTSLQQR) # normal
    hist(d$CSTN3.ESLLLDTTSLQQR[d$PTGENDER == "Female"]) # normal
    hist(d$CSTN3.ESLLLDTTSLQQR[d$PTGENDER == "Male"]) # normal
    describeBy(d$CSTN3.ESLLLDTTSLQQR, group = d$DX)
    hist(d$CSTN3.ESLLLDTTSLQQR[d$DX == "CN"]) # slight left skew
    hist(d$CSTN3.ESLLLDTTSLQQR[d$DX == "MCI"]) # normal 
    hist(d$CSTN3.ESLLLDTTSLQQR[d$DX == "Dementia"]) # odd looking, not enough data points?
      length(na.exclude(d$CSTN3.ESLLLDTTSLQQR[d$DX == "Dementia"])) # 66, seems enough people
    
  # NEGR1: NEGR1.SSIIFAGGDK
    describeBy(d$NEGR1.SSIIFAGGDK, group = d$PTGENDER)
    hist(d$NEGR1.SSIIFAGGDK) # normal
    hist(d$NEGR1.SSIIFAGGDK[d$PTGENDER == "Female"]) # normal
    hist(d$NEGR1.SSIIFAGGDK[d$PTGENDER == "Male"]) # normal
    describeBy(d$NEGR1.SSIIFAGGDK, group = d$DX)
    hist(d$NEGR1.SSIIFAGGDK[d$DX == "CN"]) # normal but outlier on low end
    hist(d$NEGR1.SSIIFAGGDK[d$DX == "MCI"]) # normalish, slight left skew 
    hist(d$NEGR1.SSIIFAGGDK[d$DX == "Dementia"]) # slightly left skewed
    
# Environmental/clinical associations with proteins ====
# checking sex differences in variables 
plot(d$APOE4) # unequal group sizes but still a lot in each 
  plot(d$APOE4[d$PTGENDER == "Male"])
  plot(d$APOE4[d$PTGENDER == "Female"]) # look quite similar between genders 
plot(d$SMOK) # less have smoked than never smoked
  plot(d$SMOK[d$PTGENDER == "Male"])
  plot(d$SMOK[d$PTGENDER == "Female"]) # likely sex difference but hard to tell as bar chart 
hist(d$BMI) # normal but some extreme values on high end
  describeBy(d$BMI, group=d$PTGENDER) # mild sex difference
hist(d$PTEDUCAT) # strongly left skewed 
  hist(d$PTEDUCAT[d$PTGENDER == "Male"])
  hist(d$PTEDUCAT[d$PTGENDER == "Female"])
  describeBy(d$PTEDUCAT, group = d$PTGENDER) # mild sex difference
hist(d$DIGITSCOR) # normal
  describeBy(d$DIGITSCOR, group=d$PTGENDER) # mild sex difference
hist(d$VSBPSYS) # roughly normal
  describeBy(d$VSBPSYS, group=d$PTGENDER) # no likely sex difference
hist(d$VSBPDIA) # roughly normal
  describeBy(d$VSBPDIA, group=d$PTGENDER) # no likely sex difference
hist(d$VSPULSE) # roughly normal
  describeBy(d$VSPULSE, group=d$PTGENDER) # moderate sex difference 
    
# CCL5 associations
  CCL5 <- list()
  CCL5$education <- lm(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. ~ PTEDUCAT + AGE + DX, na.action=na.exclude, data=d)
  CCL5$APOE4 <- lm(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. ~ APOE4 + AGE + DX, na.action=na.exclude, data=d)
  CCL5$BMI <- lm(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. ~ BMI + AGE + DX, na.action=na.exclude, data=d)
  CCL5$smoking <- lm(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. ~ SMOK + AGE + DX, na.action=na.exclude, data=d)
  CCL5$pulse <- lm(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. ~ VSPULSE + AGE + DX, na.action=na.exclude, data=d)
  CCL5$diastolic_bp <- lm(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. ~ VSBPDIA + AGE + DX, na.action=na.exclude, data=d)
  CCL5$systolic_bp <- lm(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. ~ VSBPSYS + AGE + DX, na.action=na.exclude, data=d)
  CCL5$RAVLT <- lm(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. ~ RAVLT.immediate + AGE + DX, na.action=na.exclude, data=d)
  CCL5$digit_score <- lm(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL. ~ DIGITSCOR + AGE + DX, na.action=na.exclude, data=d)
  
  CCL5_results <- do.call(cbind, lapply(CCL5, function(z) 
    summary(z)$coefficients[2,])) # aggregate model summaries  
  write.csv(CCL5_results, "CCL5_results.csv") # save 
  
  lapply(CCL5, function(x)
    plot(x)) # show diagnostic plots for all proteins. Notes taken elsewhere 
  
  meff(CCL5_results[,1:9]) # 5.194087
  # 0.05 / 5.194087 = critical p < .0096
  # correcting for anything > 2 would make observed p values non-significant as they are all > 0.025
  
  hist(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.)
  qqPlot(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.)
  hist(sqrt(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.))
  qqPlot(sqrt(d$T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.))
  
  # sqrt transformation
  CCL5sqrt <- list()
  CCL5sqrt$education <- lm(sqrt(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) ~ PTEDUCAT + AGE + DX, na.action=na.exclude, data=d)
  CCL5sqrt$APOE4 <- lm(sqrt(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) ~ APOE4 + AGE + DX, na.action=na.exclude, data=d)
  CCL5sqrt$BMI <- lm(sqrt(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) ~ BMI + AGE + DX, na.action=na.exclude, data=d)
  CCL5sqrt$smoking <- lm(sqrt(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) ~ SMOK + AGE + DX, na.action=na.exclude, data=d)
  CCL5sqrt$heart_rate <- lm(sqrt(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) ~ VSPULSE + AGE + DX, na.action=na.exclude, data=d)
  CCL5sqrt$diastolic_bp <- lm(sqrt(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) ~ VSBPDIA + AGE + DX, na.action=na.exclude, data=d)
  CCL5sqrt$systolic_bp <- lm(sqrt(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) ~ VSBPSYS + AGE + DX, na.action=na.exclude, data=d)
  CCL5sqrt$RAVLT <- lm(sqrt(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) ~ RAVLT.immediate + AGE + DX, na.action=na.exclude, data=d)
  CCL5sqrt$digit_score <- lm(sqrt(T.Cell.Specific.Protein.RANTES..RANTES...ng.mL.) ~ DIGITSCOR + AGE + DX, na.action=na.exclude, data=d)
    
  CCL5sqrt_results <- t(do.call(cbind, lapply(CCL5sqrt, function(z)
    summary(z)$coefficients[2,])))
  write.csv(CCL5sqrt_results, "CCL5sqrt_results.csv")
  
  lapply(CCL5sqrt, function(x)
    plot(x))
  
  
# CLSTN3 associations 
  CLSTN3 <- list()
  CLSTN3$education <- lm(CSTN3.ESLLLDTTSLQQR ~ PTEDUCAT + AGE + DX, na.action=na.exclude, data=d)
  CLSTN3$APOE4 <- lm(CSTN3.ESLLLDTTSLQQR ~ APOE4 + AGE + DX, na.action=na.exclude, data=d)
  CLSTN3$BMI <- lm(CSTN3.ESLLLDTTSLQQR ~ BMI + AGE + DX, na.action=na.exclude, data=d)
  CLSTN3$smoking <- lm(CSTN3.ESLLLDTTSLQQR ~ SMOK + AGE + DX, na.action=na.exclude, data=d)
  CLSTN3$pulse <- lm(CSTN3.ESLLLDTTSLQQR ~ VSPULSE + AGE + DX, na.action=na.exclude, data=d)
  CLSTN3$diastolic_bp <- lm(CSTN3.ESLLLDTTSLQQR ~ VSBPDIA + AGE + DX, na.action=na.exclude, data=d)
  CLSTN3$systolic_bp <- lm(CSTN3.ESLLLDTTSLQQR ~ VSBPSYS + AGE + DX, na.action=na.exclude, data=d)
  CLSTN3$RAVLT <- lm(CSTN3.ESLLLDTTSLQQR ~ RAVLT.immediate + AGE + DX, na.action=na.exclude, data=d)
  CLSTN3$digit_score <- lm(CSTN3.ESLLLDTTSLQQR ~ DIGITSCOR + AGE + DX, na.action=na.exclude, data=d)
    
  CLSTN3_results <- do.call(cbind, lapply(CLSTN3, function(z) 
    summary(z)$coefficients[2,]))
  write.csv(CLSTN3_results, "CLSTN3_results.csv")

  lapply(CLSTN3, function(x)
    plot(x)) # show diagnostic plots for all proteins. Notes taken elsewhere 
  
  meff(CLSTN3_results[,1:9]) # 3.932619
  # 0.05 / 3.932619 = critical p < 0.013
  
  
# NEGR1 associations
  NEGR1 <- list()
  NEGR1$education <- lm(NEGR1.SSIIFAGGDK ~ PTEDUCAT + AGE + DX, na.action=na.exclude, data=d)
  NEGR1$APOE4 <- lm(NEGR1.SSIIFAGGDK ~ APOE4 + AGE + DX, na.action=na.exclude, data=d)
  NEGR1$BMI <- lm(NEGR1.SSIIFAGGDK ~ BMI + AGE + DX, na.action=na.exclude, data=d)
  NEGR1$smoking <- lm(NEGR1.SSIIFAGGDK ~ SMOK + AGE + DX, na.action=na.exclude, data=d)
  NEGR1$pulse <- lm(NEGR1.SSIIFAGGDK ~ VSPULSE + AGE + DX, na.action=na.exclude, data=d)
  NEGR1$diastolic_bp <- lm(NEGR1.SSIIFAGGDK ~ VSBPDIA + AGE + DX, na.action=na.exclude, data=d)
  NEGR1$systolic_bp <- lm(NEGR1.SSIIFAGGDK ~ VSBPSYS + AGE + DX, na.action=na.exclude, data=d)
  NEGR1$RAVLT <- lm(NEGR1.SSIIFAGGDK ~ RAVLT.immediate + AGE + DX, na.action=na.exclude, data=d)
  NEGR1$digit_score <- lm(NEGR1.SSIIFAGGDK ~ DIGITSCOR + AGE + DX, na.action=na.exclude, data=d)
    
  NEGR1_results <- t(do.call(cbind, lapply(NEGR1, function(z) 
    summary(z)$coefficients[2,])))
  write.csv(NEGR1_results, "NEGR1_results.csv")
  
  lapply(NEGR1, function(x)
    plot(x)) # show diagnostic plots for all proteins. Notes taken elsewhere 
    

# Example structural equation modelling ====
# Example of code used if SEM were to be performed
m1 <- ' 
  phenotype ~ c * protein
  environmentalVar ~ a * protein
	phenotype ~ b * environmentalVar
  ab := a * b
  c prime := c + (a * b)
  '
fittedmodel <- sem(m1, data=d, test="bootstrap") # bootstrapped
r <- parameterEstimates(fittedmodel) # CIs
r$pvalue[8] # get p value for c

# Example of code used if mediation were to be performed 
library(rockchalk)
  m1.5 <- lm(phenotype ~ protein * scale(environmentalVar), data=d)
  summary(m1.5)

  simpleSlopes <- plotSlopes(m1.5, plotx="protein", modx="environmentalVar", modxVals="std.dev.")
  testSlopes(simpleSlopes) # statistical test
    
  