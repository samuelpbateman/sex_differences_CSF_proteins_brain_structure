getRversion() # 4.2.1
### installing ADNIMERGE and loading libraries
install.packages("./ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
library(ADNIMERGE)
library(tidyverse)
library(skimr)
library(psych)
library(summarytools)

# load csv 
finalData <- read_csv("finaldata2.csv", na = c("", "NA", -4), ) # 439 vars, 818 obs
# 

### adnimerge wrangling ======
data(adnimerge) # 16,294 cases, 115 vars
adnimerge <- adnimerge %>% filter(COLPROT == "ADNI1") # 5013 cases, 115 vars
adnimerge <- adnimerge %>% filter(VISCODE == "bl") # baseline only, 819 cases, 115 vars

# checking and deleting empty variables
adnimerge %>% select() %>% skim() # this line used to check completion for each var before deletion
glimpse(adnimerge)

adnimergeFinal <- adnimerge %>% mutate(EcogPtMem.bl = NULL, EcogPtLang.bl = NULL, EcogPtVisspat.bl = NULL, 
  EcogPtPlan.bl = NULL, EcogPtOrgan.bl = NULL, EcogPtDivatt.bi = NULL, EcogPtTotal.bl = NULL, 
  EcogSPMem.bl = NULL, EcogSPLang.bl = NULL, EcogSPVisspat.bl = NULL, EcogSPPlan.bl = NULL, 
  EcogSPOrgan.bl = NULL, EcogSPDivatt.bl = NULL, EcogSPTotal.bl = NULL, M = NULL, Month = NULL,
  Years.bl = NULL, Month.bl = NULL, MOCA.bl = NULL, FBB.bl = NULL, AV45.bl = NULL, EcogPtDivatt.bl = NULL,
  EcogPtMem = NULL, EcogPtLang = NULL, EcogPtVisspat = NULL, EcogPtPlan = NULL, EcogPtOrgan = NULL,
  EcogPtDivatt = NULL, EcogPtTotal = NULL, EcogSPMem = NULL, EcogSPLang = NULL, EcogSPVisspat = NULL,
  EcogSPPlan = NULL, EcogSPOrgan = NULL, EcogSPDivatt = NULL, EcogSPTotal = NULL, MOCA = NULL,
  AV45 = NULL, FBB = NULL, RAVLT.immediate = adnimerge$RAVLT.immediate
  ) # removes all empty variables, 77 vars remaining

adnimergeFinal <- adnimergeFinal %>% mutate(RID = factor(RID), COLPROT = NULL, ORIGPROT = NULL,
  PTID = NULL, VISCODE = NULL, EXAMDATE = NULL, .bl = NULL, AGE = as.double(AGE), 
  PTGENDER = factor(PTGENDER), PTEDUCAT = factor(PTEDUCAT), PTMARRY = factor(PTMARRY), 
  APOE4 = factor(APOE4), FDG = NULL, PIB = NULL, ABETA = NULL,
  TAU = as.double(TAU), CDRSB = NULL, ADAS11 = NULL)

# refine variables 
adnimergeFinal <- adnimergeFinal %>% select(RID, SITE, AGE, PTGENDER, PTEDUCAT, PTETHCAT, PTRACCAT, 
  PTMARRY, APOE4, Ventricles, Hippocampus, WholeBrain, Entorhinal, Fusiform, MidTemp, ICV, DX, DIGITSCOR, 
  DIGITSCOR.bl, RAVLT.immediate)

# change variable types
adnimergeFinal <- adnimergeFinal %>% mutate(Ventricles = as.double(Ventricles), 
                                            Hippocampus = as.double(Hippocampus),
                                            WholeBrain = as.double(WholeBrain),
                                            Entorhinal = as.double(Entorhinal),
                                            Fusiform = as.double(Fusiform),
                                            MidTemp = as.double(MidTemp),
                                            ICV = as.double(ICV),
                                            DIGITSCOR = as.integer(DIGITSCOR),
                                            RAVLT.immediate = as.integer(RAVLT.immediate))

# correct for ICV
adnimergeFinal <- adnimergeFinal %>% mutate(WholeBrain_corrected = WholeBrain / ICV)
adnimergeFinal <- adnimergeFinal %>% mutate(Hippocampus_corrected = Hippocampus / ICV)
adnimergeFinal <- adnimergeFinal %>% mutate(Fusiform_corrected = Fusiform / ICV)
adnimergeFinal <- adnimergeFinal %>% mutate(Ventricles_corrected = Ventricles / ICV)

### MRM proteins wrangling =====
data(csfmrm) # 303 cases, 322 cols
csfMRM <- data.frame(csfmrm)
csfMRM %>% select(everything()) %>% distinct(RID) %>% glimpse() # 288 unique of 304 total cases
csfMRM <- csfMRM[!duplicated(csfMRM$RID), ] # remove duplicates in RID, 288 cases remaining
csfMRM %>% select(ORIGPROT) %>% skim() # all are ADNI1 so no filter required
# no empty variables by visual inspection

# centering MRM protein values
cent_csfMRM <- csfMRM %>% select(-RID, -VISCODE, -ORIGPROT) %>% scale() # centering, 320 vars
cent_csfMRM <- data.frame(cent_csfMRM) # convert to dataframe
id_csfMRM <- csfMRM %>% select(RID, VISCODE, ORIGPROT) # extract vars to append, 3 vars
cent_csfMRM <- add_column(cent_csfMRM, id_csfMRM, .before = T) # combine to 323 vars



### RBM proteins wrangling ======
csfRBM <- read_csv("CSFRBM.csv") # 327 cases, 163 vars
csfRBM <- rename(csfRBM, "RID" = "rid") # renaming RID to match case of other dataframes
csfRBM %>% select(everything()) %>% distinct(RID) %>% glimpse() # 311 unique of 327 total cases
csfRBM <- csfRBM[!duplicated(csfRBM$RID), ] # remove duplicates in RID, 311 cases remaining

# checking and deleting empty variables
csfRBM %>% select(everything(where(is.numeric))) %>% skim()

# centering RBM protein values
cent_csfRBM <- csfRBM %>% select(where(is.numeric), -RID, -sampid, -visit_code, -id) %>% 
  scale() # centering and scaling all numeric variables, 83 vars
cent_csfRBM <- data.frame(cent_csfRBM)
id_csfRBM <- csfRBM %>% select(RID, sampid, visit_code, id) # extract vars to append, 4 vars
cent_csfRBM <- add_column(cent_csfRBM, id_csfRBM, .before = T) # combine to 87 vars
 


### vitals wrangling =====
vitals %>% select(Phase, VISCODE, VSHEIGHT, VSWEIGHT) %>% 
  View() # for height: bl == -4 == NA

vitals <- read_csv("VITALS.csv") # "VITALS.csv", 15,049 cases, 22 vars
vitals_ADNI1 <- subset(vitals, Phase == "ADNI1") # only ADNI1, reduce to 5044 cases, 22 vars
vitalsBl <- subset(vitals_ADNI1, VISCODE != "bl") # only baseline, reduce to 4225 cases, 22 vars
vitalsSc <- subset(vitalsBl, VISCODE == "sc")
vitalsOfInterest <- vitalsSc %>% select(RID, ID, VSHEIGHT, VSHTUNIT, VSWEIGHT, VSWTUNIT, VISCODE, 
                                        VSBPSYS, VSBPDIA, VSPULSE) # height and ID, 4425 cases, 6 vars
vitalsFinal <- vitalsOfInterest %>% arrange(RID) # one row per participant
vitalsFinal <- vitalsFinal %>% select(RID, VSHEIGHT, VSHTUNIT, VSWEIGHT, VSWTUNIT, VSPULSE, VSBPDIA, 
                                      VSBPSYS)
vitalsFinal$RID <- factor(vitalsFinal$RID)
vitalsFinal$VSHTUNIT <- factor(vitalsFinal$VSHTUNIT)
vitalsFinal$VSWTUNIT <- factor(vitalsFinal$VSWTUNIT)

vitalsFinal$VSHEIGHT[which(vitalsFinal$VSHTUNIT == 1)] <- vitalsFinal$VSHEIGHT[which(vitalsFinal$VSHTUNIT == 1)] * 2.54
# convert pounds into kilograms
vitalsFinal$VSWEIGHT[which(vitalsFinal$VSWTUNIT == 1)] <- vitalsFinal$VSWEIGHT[which(vitalsFinal$VSWTUNIT == 1)] * 0.45359237
# creating BMI variable
vitalsFinal$BMI <- vitalsFinal$VSWEIGHT / ((vitalsFinal$VSHEIGHT/100)^2)

### strokesum wrangling ====
# get white matter hyperintensities
strokeSum <- read_csv("STROKESUM.csv") # STROKESUM.csv, 3469 rows, 10 vars
skim(strokeSum)
strokeSumOfInterest <- subset(strokeSum, VISCODE == "sc") # removes later measurements, 838 rows
WMHyper <- strokeSumOfInterest %>% select(RID, VISCODE, WHITMATHYP, MAGSTRENGTH) # 838 rows, 4 vars
WMHyper <- WMHyper %>% arrange(RID) # one row per participant
WMHyper <- WMHyper %>% select(RID, WHITMATHYP, MAGSTRENGTH)
WMHyper$RID <- factor(WMHyper$RID)

### medical history wrangling (with imported data) ====
med <- read_csv("MEDHIST.csv") # 
skim(med)
smoking <- med %>% select(RID, SITEID, MH16SMOK) # 3083 rows
glimpse(finalData) # imported data has 818 rows, 439 vars 


# merge smoking variable 
finalData_smoking <- merge(finalData, smoking, by="RID", all=F) # 1293 rows, 441 vars
finalData_smoking %>% distinct(RID) %>% pull() %>% length() # 817 unique RIDs
# remove non-unique RIDs
finalData_smoking_unique <- finalData_smoking %>% select(everything()) %>% distinct(RID, .keep_all=T)
# remove SITEID and change name of smoking variable 
finalData_smoking_unique <- finalData_smoking_unique %>% mutate(SITEID = NULL, SMOK = MH16SMOK)
# delete MH16SMOK
finalData_smoking_unique <- finalData_smoking_unique %>% mutate(MH16SMOK = NULL) 

# check:
glimpse(finalData_smoking_unique) # 817 rows, 440 vars

# relocate MH16SMOK to 
names(finalData_smoking_unique) # see variable names 
fD_relocate <- finalData_smoking_unique %>% relocate(SMOK, .after = BMI) # relocate SMOK
glimpse(fD_relocate) # 817 rows, 440 vars, SMOK moved to demographic area

write.csv(fD_relocate, "finaldata2.csv")

### merging dataframes =====

# adnimerge + vitals
ADNI_Vitals <- merge(adnimergeFinal, vitalsFinal, by="RID", all=T) # 822 rows, 24 vars
# + smoking
ADNI_Vitals_2 <- merge(ADNI_Vitals, smoking, by="RID", all=T)

# adni + white matter hyp
ADNI_Final <- merge(ADNI_Vitals_2, WMHyper, by="RID", all=T) # 843 rows, 26 vars

# RBM + CSF 
proteins_Final <- merge(cent_csfMRM, cent_csfRBM, by="RID", all=T) # 311 rows, 409 vars

# adni + proteins
finalData <- merge(ADNI_Final, proteins_Final, by="RID", all=T) # 844 rows, 434 vars

finalData %>% distinct(RID) %>% pull() %>% 
  length() # 843 unique participants of 844 cases, remaining case likely labels

write.csv(finalData, "finaldata2.csv")

### data exploration =====

(dataSummary <- dfSummary(finalData[,1:36])) # inspects all non-protein variables
# VSHEIGHT, WHITMATHYP not normally distributed 

# delete unnecessary variables
  finalData <- finalData %>% mutate(DIGITSCOR.bl = NULL, sampid = NULL, visit_code = NULL, id = NULL, SITEID = NULL)
  finalData <- finalData %>% mutate(...1 = NULL, ...2 = NULL, ...3 = NULL)

# delete 3 cases which have almost all NAs
  finalData %>% select(PTGENDER, RID) %>% filter(!complete.cases(PTGENDER)) # RIDs 460, 542, 662
  finalData2 <- subset(finalData, RID != 460 & RID != 542 & RID != 662)

# inspecting BMI
  finalData %>% filter(BMI < 10) %>% glimpse() # participant 276 is an outlier, 844 obs
  finalData2 <- finalData2 %>% filter(BMI > 10) # 822 obs
  
  finalData2 %>% filter(BMI > 100) %>% glimpse() # participant 280
  finalData3 <- finalData2 %>% filter(BMI < 99) # 821 obs
    # this is now saved to finalData.csv

## saving 
  write.csv(new_df, "finaldata3.csv")


# Winsorize proteins 
  finalData <- read_csv("finaldata2.csv", na = c("", "NA", -4), ) # 439 vars, 818 obs
  df <- finalData
  df[39:441] <- DescTools::Winsorize(df[39:441], maxval=3, minval=-3, na.rm=F)
  # check if it worked
    describe(finalData$A1AT.AVLTIDEK)
    describe(df$A1AT.AVLTIDEK)
      # success, saved to finalData2.csv

# removing BMI outliers
  finalData %>% filter(BMI > 50) %>% glimpse()
  finalData <- finalData %>% filter(BMI < 50)
  skim(finalData$BMI) # success, saved to finalData.csv 

# removing unnecessary cols
  finalData <- finalData[3:437]

# removing 3 NAs for almost all vars
  finalData %>% filter(!complete.cases(AGE)) %>% glimpse()
  finalData2 <- finalData %>% filter(complete.cases(AGE))
  dim(finalData) # 821 
dim(finalData2) # 818, success, saved to finalData.csv