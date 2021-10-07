dbds1 <- read.table(file='moslemi_dbds1_questionnaire_decode_subset_sf12_smoke_drink_sf12_bmi.tsv', sep='\t', header = T, stringsAsFactors = F, fill  = T)
dbds2 <- read.table(file='moslemi_dbds2_questionnaire_decode_subset_sf12_smoke_drink_bmi.tsv', sep='\t', header = T, stringsAsFactors = F, fill  = T)

dbds1 <- dbds1[!is.na(dbds1$q32),]
dbds1 <- dbds1[!is.na(dbds1$q36),]

dbds1.bmi <- data.frame('cpr_enc'=dbds1$cpr_enc,'age'= as.integer(floor(difftime(as.Date(dbds1$donation_date),as.Date(dbds1$birthdate),units = 'days')/365.25)),'bmi'=dbds1$q32/((dbds1$q36/100)^2))

dbds2 <- dbds2[!is.na(dbds2$Blandet_6),]
dbds2 <- dbds2[!is.na(dbds2$Blandet_7),]

dbds2.bmi <- data.frame('cpr_enc'=dbds2$cpr_enc,'age'= as.integer(floor(difftime(as.Date(dbds2$date),as.Date(dbds2$birthdate),units = 'days')/365.25)),'bmi'=dbds2$Blandet_7/((dbds2$Blandet_6/100)^2))



#hist(dbds1.bmi$bmi)
#hist(dbds2.bmi$bmi)

#nrow(dbds1.bmi)
#nrow(dbds2.bmi)

#length(unique(dbds1.bmi$cpr_enc, dbds2.bmi$cpr_enc))

vel <- read.table(file='predictions.tsv', sep='\t', header = T, stringsAsFactors = F, fill  = T)
#head(vel)
#unique(vel$serologi)
#unique(vel$algorithm)

par(mfrow=c(1,2))
#-- Taqman
#nrow(vel[vel$algorithm < 0,])
vel.pos <- vel[which(vel$serologi > 0),]$donor
vel.neg <- vel[which(vel$serologi < 0),]$donor

#nrow(dbds1.bmi[dbds1.bmi$cpr_enc %in% vel.neg,])
#nrow(dbds2.bmi[dbds2.bmi$cpr_enc %in% vel.neg,])

#dbds1.bmi[dbds1.bmi$cpr_enc %in% vel.neg,]
#dbds2.bmi[dbds2.bmi$cpr_enc %in% vel.neg,]

#mean(dbds1.bmi[dbds1.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]$bmi)
#mean(dbds2.bmi[dbds2.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]$bmi
#boxplot(dbds1.bmi[dbds1.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]$bmi,dbds2.bmi[dbds2.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]$bmi)

dbds1.bmi <- dbds1.bmi[ !dbds1.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]
dbds12.bmi <- rbind(dbds1.bmi,dbds2.bmi)

nrow(dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.pos,])
nrow(dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.neg,])

ages <- dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.neg,]$age
#################### END ####################

boxplot(dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.pos & dbds12.bmi$age %in% ages,]$bmi,dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.neg,]$bmi, names = c("VEL.POS", "VEL.NEG"), main="Taqman")
legend(x="topright", legend = c("POS: 6601","NEG: 15", "P = 0.8261"))
t.test(dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.pos & dbds12.bmi$age %in% ages,]$bmi,dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.neg,]$bmi)

# Algorithm
vel.pos <- vel[which(vel$algorithm == 1),]$donor
vel.neg <- vel[which(vel$algorithm == -1),]$donor

ages <- dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.neg,]$age

#nrow(dbds1.bmi[dbds1.bmi$cpr_enc %in% vel.neg,])
#nrow(dbds2.bmi[dbds2.bmi$cpr_enc %in% vel.neg,])

#dbds1.bmi[dbds1.bmi$cpr_enc %in% vel.neg,]
#dbds2.bmi[dbds2.bmi$cpr_enc %in% vel.neg,]

#mean(dbds1.bmi[dbds1.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]$bmi)
#mean(dbds2.bmi[dbds2.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]$bmi)

#boxplot(dbds1.bmi[dbds1.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]$bmi,dbds2.bmi[dbds2.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]$bmi)

#dbds1.bmi <- dbds1.bmi[ !dbds1.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc),]
#dbds12.bmi <- rbind(dbds1.bmi,dbds2.bmi)

nrow(dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.pos & dbds12.bmi$age %in% ages,])
nrow(dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.neg,])

boxplot(dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.pos & dbds12.bmi$age %in% ages,]$bmi,dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.neg,]$bmi, names = c("VEL.POS", "VEL.NEG"), main="Imputed")
legend(x="topright", legend = c("POS: 40739","NEG: 37", "P = 0.02944"))
t.test(dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.pos & dbds12.bmi$age %in% ages,]$bmi,dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.neg,]$bmi)

##############################
diagnosis <- read.table(file='/data/projects/rbc_gwas/moslemi_registry.tsv', sep=' ', quote="\"", header = T, stringsAsFactors = F, fill  = T)
diagnosis.selection <- diagnosis[diagnosis$cpr_enc %in% vel.neg,]
length(unique(diagnosis$cpr_enc))
unique(diagnosis$diagnosis)


### BMI - smoke

read.table(file='/data/projects/rbc_gwas/DBDS1_questions_explained.txt', sep=' ', quote="\"", header = T, stringsAsFactors = F, fill  = T)
read.table('/data/projects/rbc_gwas/dbds2_variables.csv', sep=' ', quote="\"", header = T, stringsAsFactors = F, fill  = T)

#dbds3.variables<-read.table('/data/projects/rbc_gwas/dbds3_variables.tsv', sep='\t', header = T, stringsAsFactors = F)

dbds3<-read.table('/data/projects/rbc_gwas/dbds3_rbc_gwas.tsv', sep='\t', header = T, stringsAsFactors = F)

### DBDS3 smoke
smokers<- dbds3[which(!(dbds3$Smoking_1 == 'NaN') & dbds3$Smoking_1 %in% c("EC4","EC1","EC3","EC2")),]$cpr_enc
smokers<-unique(c(smokers,dbds3[which(!(dbds3$Smoking_2 == 'NaN') & dbds3$Smoking_2 %in% c("ED4","ED2","ED1","ED3")),]$cpr_enc))

dbds3.smoke<-data.frame('cpr_enc'=dbds3$cpr_enc[which(!(dbds3$Smoking_1 == 'NaN') & !(dbds3$Smoking_2 == 'NaN'))],'smoker'=FALSE,stringsAsFactors = F)
dbds3.smoke[which(dbds3.smoke$cpr_enc %in% smokers),]$smoker=T

### DBDS3 bmi
dbds3.bmi<-dbds3[,c(1,4,5)]
dbds3.bmi<-dbds3.bmi[which(!(dbds3.bmi$Weight == 'NaN' | dbds3.bmi$Height == 'NaN') ),]
dbds3.bmi[which(is.na(as.numeric(dbds3.bmi$Weight))),]$Weight = '150'
dbds3.bmi$Weight <- as.numeric(dbds3.bmi$Weight)
dbds3.bmi$Height <- as.numeric(dbds3.bmi$Height)
dbds3.bmi<-data.frame('cpr_enc'=dbds3.bmi$cpr_enc, 'bmi'=dbds3.bmi$Weight/ ((dbds3.bmi$Height/100)^2),stringsAsFactors = F) 

dbds3<-merge(dbds3.bmi,dbds3.smoke,by='cpr_enc')

#nonsmoker <- dbds3[which(!(dbds3$Smoking_1 == 'NaN') & dbds3$Smoking_1 %in% c("EC5")),]
#nonsmoker <- dbds3[which(!(dbds3$Smoking_2 == 'NaN') & dbds3$Smoking_2 %in% c("ED5")),]

#dbds3[which(!(dbds3$Smoking_1 == 'NaN' | dbds3$Smoking_2 == 'NaN') ),]

###
dbds1 <- read.table(file='/data/projects/rbc_gwas/moslemi_dbds1_questionnaire_decode_subset_sf12_smoke_drink_sf12_bmi.tsv', sep='\t', header = T, stringsAsFactors = F, fill  = T)
dbds2 <- read.table(file='/data/projects/rbc_gwas/moslemi_dbds2_questionnaire_decode_subset_sf12_smoke_drink_bmi.tsv', sep='\t', header = T, stringsAsFactors = F, fill  = T)

dbds1 <- dbds1[!is.na(dbds1$q32),]
dbds1 <- dbds1[!is.na(dbds1$q36),]

dbds1.bmi <- data.frame('cpr_enc'=dbds1$cpr_enc,'age'= as.integer(floor(difftime(as.Date(dbds1$donation_date),as.Date(dbds1$birthdate),units = 'days')/365.25)),'bmi'=dbds1$q32/((dbds1$q36/100)^2))

dbds2 <- dbds2[!is.na(dbds2$Blandet_6),]
dbds2 <- dbds2[!is.na(dbds2$Blandet_7),]

dbds2.bmi <- data.frame('cpr_enc'=dbds2$cpr_enc,'age'= as.integer(floor(difftime(as.Date(dbds2$date),as.Date(dbds2$birthdate),units = 'days')/365.25)),'bmi'=dbds2$Blandet_7/((dbds2$Blandet_6/100)^2))


### rygning
dbds1 <- read.table(file='/data/projects/rbc_gwas/moslemi_dbds1_questionnaire_decode_subset_sf12_smoke_drink_sf12_bmi.tsv', sep='\t', header = T, stringsAsFactors = F, fill  = T)
dbds2 <- read.table(file='/data/projects/rbc_gwas/moslemi_dbds2_questionnaire_decode_subset_sf12_smoke_drink_bmi.tsv', sep='\t', header = T, stringsAsFactors = F, fill  = T)

dbds1.smoker <- dbds1[which(dbds1$q17 %in% c(2,1) | dbds1$q18 %in% c(2,1)),]$cpr_enc
dbds1.nonsmoker <- dbds1[which(dbds1$q17 == 3 & dbds1$q18 == 3),]$cpr_enc

dbds2.smoker <- dbds2[which(dbds2$Blandet_1 %in% c(2,1) | dbds2$Blandet_2 %in% c(2,1)),]$cpr_enc
dbds2.nonsmoker <- dbds2[which(dbds2$Blandet_1 == 3 & dbds2$Blandet_2 == 3),]$cpr_enc
dbds2.smoker <- unique(c(dbds2.smoker,dbds2[which(dbds2$Rygning_1>0),]$cpr_enc))
dbds2.nonsmoker <- unique(c(dbds2.nonsmoker,dbds2[which(dbds2$Rygning_1==0),]$cpr_enc))

length(dbds1.smoker)/(length(dbds1.smoker)+length(dbds1.nonsmoker))
length(dbds2.smoker)/(length(dbds2.smoker)+length(dbds2.nonsmoker))
