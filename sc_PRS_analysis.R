#scRNAseq informed polygenic risk scoring
#Mo Dehestani, 27.05.2024


.libPaths( c( "/data/Common_Folder/R/Single_cell_packages/", .libPaths()) )

library(data.table)
library(dplyr)

#load up annotated gwas genes
gwasgenes <- fread("/data/dehestani/scPRS_analysis/Annotated_PD_GWAS/gwasgenes")

#load up your DEG list 
DEGs <- read.csv("path/to/DEGs") %>% filter(p_val_adj < 0.05) 
#DEGs <- fread("/data/nasser/Manuscript/DEGs/Figure1_c_iODC_DEGs.csv")  %>% filter(p_val_adj < 0.05) 


#if there is no Gene column set it here  
DEGs$Gene <- rownames(DEGs)

#to remove SNCA and LRRK2 from DEG lists - optional 
DEGs <- subset(DEGs, !(Gene %in% c("LRRK2", "SNCA")))
 

#merge DEGs with gwas genes 
merged  <- merge(gwasgenes, DEGs, by.x = "Gene.refGene", by.y = "Gene") 
 

#define a SNP col as: chrom-pos-A1-A2 
merged$SNP <- paste(merged$Chr, merged$Start, merged$Ref, merged$Alt, sep = ":")
 

#load up PD GWAS summ stat 
GWAS <- fread("/data/dehestani/scPRS_analysis/Sum_stats/Chang2017_GWAS.tab")
#GWAS$chr_bp <- paste0("chr",GWAS$CHR,":",GWAS$BP)
#GWAS$SNP <- paste0(GWAS$CHR,":",GWAS$BP, ":", GWAS$A1, ":", GWAS$A2)
#colnames(GWAS)[colnames(GWAS) == "P.META"] <- "P"
#colnames(GWAS)[colnames(GWAS) == "SE.META"] <- "SE"
#colnames(GWAS)[colnames(GWAS) == "OR.META"] <- "OR"
#write.table(GWAS, file="/data/dehestani/scPRS_analysis/Sum_stats/Chang2017_GWAS.tab", quote=F, row.names = F, sep="\t")



#merge SNPs 
merged_SNPs <- merge(GWAS, merged, by = "SNP")
 


#create a PRSice2 input summ stat and  save it 
subset_merged_SNPs <- merged_SNPs[, c("SNP", "A1", "A2", "OR","SE","P")]
write.table(unique(subset_merged_SNPs), "/data/dehestani/scPRS_analysis/PRsice2_input_sumstats/ODC_test.txt", row.names = FALSE, quote = F, sep = "\t")
 

#### Calculate PRS using PRSice2
#run /data/dehestani/scPRS_analysis/Scripts/PRsice2_code.sh in bash terminal
 
#### CLINICAL OUTCOME PREDICTION 

library("readxl")
library(lm.beta)

PRS_all <- read_xlsx("/data/dehestani/scPRS_analysis/PRsice2_output/bestPRS")

jing <- merge(PRS_all,Covariates_file,by="IID")


#Predict clincal outcomes
my_data <- fread("/data/dehestani/scPRS_analysis/Clinical_outcomes/clinical_outcomes.csv")

my_dataV2 <- as.data.frame(my_data)
my_dataV3 <- my_dataV2[(my_dataV2$`last diagnosis` == "PD"),]
colnames(my_dataV3)[1] <- "IID" 
PRS_all_clinic <- merge(PRS_all,my_dataV3, by="IID")
Covariates_file <- read.delim("/data/dehestani/scPRS_analysis/Clinical_outcomes/Covariates", stringsAsFactors=F)
PRS_all_clinicV2 <- merge(PRS_all_clinic,Covariates_file, by="IID")
PRS_all_clinic <- PRS_all_clinicV2



#UPDRS-III prediction 

PRS_all_clinic_UPDRS <- PRS_all_clinic[PRS_all_clinic$UPDRS_III_total_score != "unknown",]
PRS_all_clinic_UPDRS$UPDRS_III_total_score <- as.numeric(PRS_all_clinic_UPDRS$UPDRS_III_total_score)

summary(lm(UPDRS_III_total_score~sex+AGE+PC1+PC2+PC3+PC4+ODC+OPC+CEP+NL1+NL2+POPC+RGC+NPC, data=PRS_all_clinic_UPDRS))
lm.beta(lm(UPDRS_III_total_score~sex+AGE+PC1+PC2+PC3+PC4+ODC+OPC+CEP+NL1+NL2+POPC+RGC+NPC, data=PRS_all_clinic_UPDRS))
 

#MoCA prediction 

PRS_all_clinic_MOCA <- PRS_all_clinic[PRS_all_clinic$MoCA_total_score != "unknown",]
PRS_all_clinic_MOCA$MoCA_total_score <- as.numeric(PRS_all_clinic_MOCA$MoCA_total_score)

summary(lm(MoCA_total_score~sex+AGE+PC1+PC2+PC3+PC4+ODC+OPC+CEP+NL1+NL2+POPC+RGC+NPC, data=PRS_all_clinic_MOCA))
lm.beta(lm(MoCA_total_score~sex+AGE+PC1+PC2+PC3+PC4+ODC+OPC+CEP+NL1+NL2+POPC+RGC+NPC, data=PRS_all_clinic_MOCA))
 

#BDI-II prediction 

PRS_all_clinic_MOCA <- PRS_all_clinic[PRS_all_clinic$`BDI-II_total score` != "unknown",]
PRS_all_clinic_MOCA$`BDI-II_total score` <- as.numeric(PRS_all_clinic_MOCA$`BDI-II_total score`)

summary(lm(`BDI-II_total score`~sex+AGE+PC1+PC2+PC3+PC4+ODC+OPC+CEP+NL1+NL2+POPC+RGC+NPC, data=PRS_all_clinic_MOCA))
lm.beta(lm(`BDI-II_total score`~sex+AGE+PC1+PC2+PC3+PC4+ODC+OPC+CEP+NL1+NL2+POPC+RGC+NPC, data=PRS_all_clinic_MOCA))

 

### K fold cross validation to figure out overfitting in linear regression

# Load necessary packages
library(caret)
library(e1071)

# Define your data and formula
data <- PRS_all_clinic_MOCA
formula <- `BDI-II_total score` ~ sex + AGE + PC1 + PC2 + PC3 + PC4 + ODC + OPC + CEP + NL1 + NL2 + POPC + RGC + NPC

# Set up train control for k-fold cross-validation
set.seed(123) # For reproducibility
train_control <- trainControl(method = "cv", number = 10) # 10-fold cross-validation

# Train the model using cross-validation
model <- train(formula, data = data, method = "lm", trControl = train_control)

# Print the results
print(model)
summary(model$finalModel)


#####  logistic regression for case/control prediction 
library(pROC)

PRS_cov <- merge(PRS_all, Covariates_file, by = "FID")
PRS_cov$Pheno <- as.factor(PRS_cov$Pheno)
PRS_cov$SEX <- as.factor(PRS_cov$SEX)
PRS_cov <- na.omit(PRS_cov)

#run glm for your PRS model
model <- glm(Pheno ~ PRSmodel + AGE + SEX + PC1 + PC2 + PC3 + PC4, 
             data = PRS_cov, 
             family = binomial(link='logit'))
summary(model)

#check for colinearity by Variance Inflation Factor 
#VIF value around 1 means no colinearity
vif(model)


# Get predicted probabilities
predicted_probabilities <- predict(model, type = "response")


# Calculate AUC
roc_curve <- roc(PRS_cov$Pheno, predicted_probabilities)
auc <- auc(roc_curve)
print(paste("AUC:", auc))

# Plot ROC curve
plot(roc_curve, main="ROC Curve", col="blue")




# Set up cross-validation for logistic regression 
train_control <- trainControl(method = "cv", number = 10)

# Define the formula for the logistic regression model
formula <- Pheno ~ PRSmodel + AGE + SEX + PC1 + PC2 + PC3 + PC4

# Train the logistic regression model with cross-validation
model <- train(formula, data = PRS_cov, method = "glm", family = binomial(link = 'logit'), trControl = train_control)

# Print the results
print(model)




# Set up cross-validation for logistic regression 
train_control <- trainControl(method = "cv", number = 10)

# Define the formula for the logistic regression model
formula <- Pheno ~ ODC + AGE + SEX + PC1 + PC2 + PC3 + PC4

# Train the logistic regression model with cross-validation
model <- train(formula, data = PRS_cov, method = "glm", family = binomial(link = 'logit'), trControl = train_control)

# Print the results
print(model)




###do regularization to account for overfitting

# Define the matrix of predictors and the response variable
X <- as.matrix(PRS_cov[, c("ODC", "AGE", "SEX", "PC1", "PC2", "PC3", "PC4")])
PRS_cov$Pheno <- factor(PRS_cov$Pheno, levels = c("1", "2"), labels = c("Class1", "Class2"))
y <- PRS_cov$Pheno

# Set up cross-validation
train_control <- trainControl(method = "cv", 
                              number = 10, 
                              classProbs = TRUE, 
                              summaryFunction = twoClassSummary)

# Train the logistic regression model with Lasso regularization (alpha = 1)
lasso_model <- train(X, y, 
                     method = "glmnet", 
                     trControl = train_control, 
                     tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 1, by = 0.01)), 
                     family = "binomial")

# Train the logistic regression model with Ridge regularization (alpha = 0)
ridge_model <- train(X, y, 
                     method = "glmnet", 
                     trControl = train_control, 
                     tuneGrid = expand.grid(alpha = 0, lambda = seq(0.001, 1, by = 0.01)), 
                     family = "binomial")

# Print the results
print(lasso_model)
print(ridge_model)

# Best tuning parameters
print(lasso_model$bestTune)
print(ridge_model$bestTune)


# Lasso Model
lasso_preds <- predict(lasso_model, newdata = PRS_cov)
confusionMatrix(lasso_preds, PRS_cov$Pheno)

# Ridge Model
ridge_preds <- predict(ridge_model, newdata = PRS_cov)
confusionMatrix(ridge_preds, PRS_cov$Pheno)


# Predict probabilities
lasso_probs <- predict(lasso_model, newdata = PRS_cov, type = "prob")
ridge_probs <- predict(ridge_model, newdata = PRS_cov, type = "prob")

# ROC curve for lasso model
roc_lasso <- roc(PRS_cov$Pheno, lasso_probs$Class1)
plot(roc_lasso, col = "blue", main = "ROC Curves for Lasso and Ridge Models")

# ROC curve for ridge model
roc_ridge <- roc(PRS_cov$Pheno, ridge_probs$Class1)
plot(roc_ridge, col = "red", add = TRUE)

# Add legend
legend("bottomright", legend = c("Lasso", "Ridge"), col = c("blue", "red"), lwd = 2)





pd <- readRDS("/data/nasser/Manuscript/processedobject/ODC35_woClus8_subclust3_res0.15_NK")
table(a@meta.data$BroadCellType) 



