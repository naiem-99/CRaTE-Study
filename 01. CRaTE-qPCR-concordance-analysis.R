#------------------- Set Directory------------------------------------------------
setwd("D:/ICDDRB_Feb/CRaTE_Study/1.qPCR-Manuscript")
getwd()
#-----------------load packages----------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(stringr)
library(tibble)
library(vroom)
library(data.table)
library(irr)
rm(list = ls());gc()
#--------------------------------Data Loading---------------------------------------
# 1. Bioperfectus Data
CRaTE_Bio <- read.csv("D:/ICDDRB_Feb/CRaTE_Study/1.qPCR-Manuscript/tidy/qpcrbio.csv", check.names = FALSE) %>%
  select(patient_id, fam_result_obs, rox_result_obs, vic_result_obs) %>%
  rename_with(~ paste0("Bioperfectus_", .), -patient_id) %>%
  mutate(bioperfectus_result = case_when(
    Bioperfectus_fam_result_obs == "positive" &
      Bioperfectus_rox_result_obs == "positive" &
      Bioperfectus_vic_result_obs == "negative" ~ 1,
    TRUE ~ 0))
#---------------------------------------------------------------------------
# 2. Culture Data
CRaTE_Culture_D <- read.csv("D:/ICDDRB_Feb/CRaTE_Study/1.qPCR-Manuscript/tidy/culture.csv", check.names = FALSE) %>%
  filter(direct_enriched == "direct") %>%
  select(patient_id, inaba_agglut, ogawa_agglut) %>%
  rename_with(~ paste0("Culture_", .), -patient_id) %>%
  mutate(
    # Combined result: returns 1 if either is yes, 0 if neither, and 0 if missing (NA)
    culture_result = case_when(Culture_inaba_agglut == "yes" | Culture_ogawa_agglut == "yes" ~ 1,
      is.na(Culture_inaba_agglut) & is.na(Culture_ogawa_agglut) ~ 0,
      TRUE ~ 0))
#-----------------------------------------------------------------------
# 3. PCR Stool Data (Separated Columns)
CRaTE_PCR_D_Stool <- read.csv("D:/ICDDRB_Feb/CRaTE_Study/1.qPCR-Manuscript/tidy/pcr_1816.csv", check.names = FALSE) %>%
  filter(`metadata-dna_extracted_from` == "stool",
         `metadata-direct_enriched` == "direct",
         `metadata-pcr_time_point` == "week1") %>%
  select(patient_id, ctxa_result, rfbO1_result) %>%
  rename_with(~ paste0("PCR_", .), -patient_id) %>% 
  mutate(pcr_ctxa_result = if_else(PCR_ctxa_result == "positive", 1, 0),
         pcr_rfbO1_result = if_else(PCR_rfbO1_result == "positive", 1, 0))
#-------------------------------------------------------------------------------
# 4. LSHTM Data (Separated Columns)
CRaTE_LSHTM <- read.csv("D:/ICDDRB_Feb/CRaTE_Study/1.qPCR-Manuscript/tidy/qpcrlshtm.csv", check.names = FALSE) %>%
  select(patient_id, fam_result_obs, rox_result_obs) %>%
  rename_with(~ paste0("LSHTM_", .), -patient_id) %>%
  mutate(lshtm_fam_result = if_else(LSHTM_fam_result_obs == "positive", 1, 0),
         lshtm_rox_result = if_else(LSHTM_rox_result_obs == "positive", 1, 0))
#-------------Merge the Data------------------------------------------------------------------
CraTE_All <- CRaTE_Bio %>%
  inner_join(CRaTE_Culture_D, by = "patient_id") %>%
  inner_join(CRaTE_PCR_D_Stool, by = "patient_id") %>%
  inner_join(CRaTE_LSHTM, by = "patient_id") %>% select(patient_id,bioperfectus_result,culture_result,pcr_ctxa_result,pcr_rfbO1_result,lshtm_rox_result,lshtm_fam_result)

write.csv(CraTE_All,"CRaTE_1000_Bioperfectus_compiled_v2.csv",row.names = F)

#------------------------------Kappa correlation-----------------------------------------------------------
library(tidyverse)
library(irr)  # kappa2()
library(corrplot)  
#-----------------------------------------------------------------------------------------------------------
d <- read.csv("CRaTE_1000_Bioperfectus_compiled_v2.csv", check.names = FALSE) %>%filter(!is.na(patient_id) & patient_id != "")

## ---- helper: build a kappa matrix from a set of 0/1 columns ----
kappa_matrix <- function(data, cols, labels) {
  n <- length(cols); M <- diag(n); dimnames(M) <- list(labels, labels)
  for (i in 1:(n-1)) for (j in (i+1):n) {
    k <- kappa2(data[, c(cols[i], cols[j])])$value
    M[i, j] <- M[j, i] <- k
  }
  M
}

## palette matching the report (blue -> white -> pink/red)
rep_col <- colorRampPalette(c("#3b6d99", "#f4eef0", "#c98a95"))(200)

## ================= PANEL A : ctxA =================
## Endpoint PCR ctxA  vs  LSHTM ctxA (= lshtm_fam)
A <- kappa_matrix(d,cols   = c("pcr_ctxa_result", "lshtm_fam_result"),labels = c("Endpoint PCR", "LSHTM qPCR"))
print(round(A, 2))



corrplot(A, method = "color", type = "upper", col = rep_col,
         col.lim = c(-1, 1), addCoef.col = "black", number.cex = 1,
         tl.col = "black", tl.srt = 45, cl.pos = "r", 
         cl.ratio = 0.2, cl.align.text = "l", cl.offset = 0.5,
         diag = TRUE, title = "A  ctxA", mar = c(0, 0, 2, 1))

## ================= PANEL B : O1 rfb =================
#pairwise Cohen's kappa in both panels)
## Endpoint PCR O1 rfb | Bioperfectus O1 rfbM | LSHTM O1 rfb (= lshtm_rox)
B <- kappa_matrix(d,cols   = c("pcr_rfbO1_result", "bioperfectus_result", "lshtm_rox_result"),
                    labels = c("Endpoint PCR", "Bioperfectus qPCR", "LSHTM qPCR"))
print(round(B, 2))

corrplot(B, method = "color", type = "upper", col = rep_col,
         col.lim = c(-1, 1), addCoef.col = "black", number.cex = 1,
         tl.col = "black", tl.srt = 45, cl.pos = "r", 
         cl.ratio = 0.2, cl.align.text = "l", cl.offset = 0.5,
         diag = TRUE, title = "B  O1 rfb", mar = c(0, 0, 2, 1))

## ---- optional: positive / negative concordance for every pair ----
concordance <- function(x, y) {
  a <- sum(x==1 & y==1); b <- sum(x==1 & y==0)
  c <- sum(x==0 & y==1); dd <- sum(x==0 & y==0)
  tibble(both_pos=a, A_only=b, B_only=c, both_neg=dd,
         pos_conc = round(100*a/(a+b+c), 1),
         neg_conc = round(100*dd/(b+c+dd), 1))
}
concordance(d$pcr_ctxa_result,  d$lshtm_fam_result)   # ctxA: EPCR vs LSHTM
concordance(d$pcr_rfbO1_result, d$bioperfectus_result) # O1: EPCR vs Bioperfectus
concordance(d$pcr_rfbO1_result, d$lshtm_rox_result)    # O1: EPCR vs LSHTM
concordance(d$bioperfectus_result, d$lshtm_rox_result) # O1: Bioperfectus vs LSHTM

## ============ SAVE BOTH HEATMAPS IN ONE FILE ============
png("Figure16_kappa_heatmaps.png", width = 11, height = 5, units = "in", res = 900)
par(mfrow = c(1, 2))              # two panels side by side

corrplot(A, method = "color", type = "upper", col = rep_col,
         col.lim = c(-1, 1), addCoef.col = "black", number.cex = 1,
         tl.col = "black", tl.srt = 45, cl.pos = "r", 
         cl.ratio = 0.2, cl.align.text = "l", cl.offset = 0.5,
         diag = TRUE, title = "A  ctxA", mar = c(0, 0, 2, 1))


corrplot(B, method = "color", type = "upper", col = rep_col,
         col.lim = c(-1, 1), addCoef.col = "black", number.cex = 1,
         tl.col = "black", tl.srt = 45, cl.pos = "r", 
         cl.ratio = 0.2, cl.align.text = "l", cl.offset = 0.5,
         diag = TRUE, title = "B  O1 rfb", mar = c(0, 0, 2, 1))

dev.off()

## ============ SAVE ALL CONCORDANCE RESULTS IN ONE TABLE ============
results <- bind_rows(
  concordance(d$pcr_ctxa_result,   d$lshtm_fam_result)   %>% mutate(pair = "ctxA: Endpoint vs LSHTM"),
  concordance(d$pcr_rfbO1_result,  d$bioperfectus_result)%>% mutate(pair = "O1: Endpoint vs Bioperfectus"),
  concordance(d$pcr_rfbO1_result,  d$lshtm_rox_result)   %>% mutate(pair = "O1: Endpoint vs LSHTM"),
  concordance(d$bioperfectus_result, d$lshtm_rox_result) %>% mutate(pair = "O1: Bioperfectus vs LSHTM")
) %>% select(pair, everything())

print(results)
write.csv(results, "kappa_concordance_results.csv", row.names = FALSE)


#------------------------------------------------------------------------------------------------
