library("survival")

reset <- FALSE

#### DBDS cohort ####
cohort <- read.table(file='cohort_overview_20200707.tsv', sep='\t', header = T, stringsAsFactors = FALSE)#, fill  = T)
cohort <- cohort[which(cohort$cohort %in% c('DBDS','DBDS_CHB')),]
cohort.male <- unique(cohort[which(cohort$sex == 'M'),]$cpr_enc)
cohort.female <- unique(cohort[which(cohort$sex == 'F'),]$cpr_enc)
cohort<-data.frame('cpr_enc'=cohort$cpr_enc, 'sex'=cohort$sex, 'birthdate'=cohort$birthdate, 'age'= as.numeric(floor(difftime(strptime(Sys.Date(), format = "%Y-%m-%d"),strptime(cohort$birthdate, format = "%Y-%m-%d"),units="days")/365.25)), stringsAsFactors = FALSE)

### LAB DATA ###
if (reset)
{
  lab <- read.table(file='lab_forsker_rbc_gwas_chb.tsv', sep='\t', header = T, stringsAsFactors = F)#, fill  = T)
  lab <- lab[which(lab$cpr_enc %in% cohort$cpr_enc),]
  save(lab,file = 'lab.dbds.Rda',compress=T)
  
  lab<-lab[,c(1,2,4,6,9,10)]
  lab<-merge(lab,cohort, by='cpr_enc')
  lab<-data.frame(lab, 'SAMPLINGAGE'= as.numeric(difftime(strptime(lab$SAMPLINGDATE, format = "%Y-%m-%d"),strptime(lab$birthdate, format = "%Y-%m-%d"),units="days")/365.25), stringsAsFactors = F)
  lab$VALUE=as.numeric(sub(",", ".", sub("<", "", lab$VALUE,fixed = TRUE), fixed = TRUE ))
  lab <- lab[which(!is.na(lab$VALUE)),]
  lab$REFERENCEINTERVAL_LOWERLIMIT=as.numeric(sub(",", ".", sub("<", "", lab$REFERENCEINTERVAL_LOWERLIMIT,fixed = TRUE), fixed = TRUE ))
  lab$REFERENCEINTERVAL_UPPERLIMIT=as.numeric(sub(",", ".", sub("<", "", lab$REFERENCEINTERVAL_UPPERLIMIT,fixed = TRUE), fixed = TRUE ))
  save(lab,file = 'lab.processed.dbds.Rda',compress=T)
} else
  #load('lab.dbds.Rda')
  load('lab.processed.dbds.Rda')

### VEL-neg cohort ###
if (reset)
{
  vel1 <- read.table(file='vel_genotypes_chb_version_20200417.tsv', sep='\t', header = TRUE, stringsAsFactors = FALSE, fill  = TRUE)
  #vel2 <- read.table(file='vel_dbds_corrected_tap_enc_20201012.tsv', sep='\t', header = TRUE, stringsAsFactors = FALSE, fill  = TRUE)
  vel2 <- read.table(file='vel_dbds_corrected_tap_enc_20201105.tsv', sep='\t', header = TRUE, stringsAsFactors = FALSE, fill  = TRUE)
  vel.neg<-unique(c(vel1[which(vel1$smim1_17bp_indel == '-:-'),]$cpr_enc, vel2[which(vel2$rs566629828 == 'Δ/Δ'),]$cpr_enc))
  vel.neg<-vel.neg[vel.neg %in% cohort$cpr_enc]
  rm(vel1)
  rm(vel2)
  save(vel.neg,file = 'vel.neg.dbds.Rda',compress=T)
} else
  load(file = 'vel.neg.dbds.Rda')

#### VEL-pos cohort ####
if (reset)
{
  #cohort.diagnosis <- read.table(file='lpr_psyk_rbc_gwas_chb.tsv', sep='\t', header = T, stringsAsFactors = FALSE)
  #cohort.diagnosis <- unique(cohort.diagnosis$cpr_enc)
  #cohort <- cohort[which(cohort$cpr_enc %in% cohort.diagnosis),]
  ## MEN ##
  vel.pos <- c()
  cohort.age<-data.frame(table("age"=cohort[which(cohort$cpr_enc %in% vel.neg & cohort$sex == 'M'),]$age), stringsAsFactors = FALSE)
  for (i in seq(nrow(cohort.age)))
  {
    age <- as.numeric(as.character(cohort.age$age[i]))
    num <- as.numeric(as.character(cohort.age$Freq[i])) * 15
    
    ids.age <- unique(cohort[which(cohort$age == age & cohort$sex == 'M' & !(cohort$cpr_enc %in% vel.neg) ),]$cpr_enc)
    picks <- c()
    while (length(picks) < num)
    {
      id <- sample(ids.age,1)
      if (! id %in% picks)
        picks <- c(picks,id)
    }
    vel.pos <- c(vel.pos,picks)
  }
  vel.pos.male <- vel.pos
  
  ## WOMEN ##
  vel.pos <- c()
  cohort.age<-data.frame(table("age"=cohort[which(cohort$cpr_enc %in% vel.neg & cohort$sex == 'F'),]$age), stringsAsFactors = FALSE)
  for (i in seq(nrow(cohort.age)))
  {
    age <- as.numeric(as.character(cohort.age$age[i]))
    num <- as.numeric(as.character(cohort.age$Freq[i])) * 15
    
    ids.age <- unique(cohort[which(cohort$age == age & cohort$sex == 'F' & !(cohort$cpr_enc %in% vel.neg)),]$cpr_enc)
    picks <- c()
    while (length(picks) < num)
    {
      id <- sample(ids.age,1)
      if (! id %in% picks)
        picks <- c(picks,id)
    }
    vel.pos <- c(vel.pos,picks)
  }
  vel.pos.female <- vel.pos
  
  vel.pos <- c(vel.pos.male,vel.pos.female)
  save(vel.pos,file = 'vel.pos.dbds.Rda',compress=FALSE)
} else
  load('vel.pos.dbds.Rda')

################### TABLE1 #####################

table1 <- data.frame(matrix(data = NA, nrow = 3, ncol = 4), stringsAsFactors = F)
names(table1)<-c('TYPE','VEL-neg cases','VEL-neg controls','remainder')
cohort.vel.neg <- cohort[which(cohort$cpr_enc %in% vel.neg),]
cohort.vel.pos <- cohort[which(cohort$cpr_enc %in% vel.pos),]
cohort.remainder <- cohort[which(!cohort$cpr_enc %in% c(vel.neg,vel.pos)),]
table1[1,] = c('size',nrow(cohort.vel.neg),nrow(cohort.vel.pos),nrow(cohort.remainder))
table1[2,] = c('mean-age',mean(cohort.vel.neg$age),mean(cohort.vel.pos$age),mean(cohort.remainder$age))
table1[3,] = c('sex-fraction (F)',
               nrow(cohort.vel.neg[which(cohort.vel.neg$sex=='F'),])/nrow(cohort.vel.neg),
               nrow(cohort.vel.pos[which(cohort.vel.pos$sex=='F'),])/nrow(cohort.vel.pos),
               nrow(cohort.remainder[which(cohort.remainder$sex=='F'),])/nrow(cohort.remainder)
)

#table1 <- data.frame(matrix(data = NA, nrow = 3, ncol = 6), stringsAsFactors = F)
#names(table1)<-c('TYPE','VEL-neg cases','VEL-neg controls','VEL-het cases','VEL-het controls','remainder')
#cohort.vel.neg <- cohort[which(cohort$cpr_enc %in% vel.neg),]
#cohort.vel.pos <- cohort[which(cohort$cpr_enc %in% vel.pos),]
#cohort.vel.het <- cohort[which(cohort$cpr_enc %in% vel.het),]
#cohort.vel.het.control <- cohort[which(cohort$cpr_enc %in% vel.het.control),]
#cohort.remainder <- cohort[which(!cohort$cpr_enc %in% c(vel.neg,vel.pos,vel.het,vel.het.control)),]
#table1[1,] = c('size',length(vel.neg),length(vel.pos))
#table1[2,] = c('mean-age',mean(cohort.vel.neg$age),mean(cohort.vel.pos$age),mean(cohort.vel.het$age),mean(cohort.vel.het.control$age),mean(cohort.remainder$age))
#table1[3,] = c('sex-fraction (F)',
               
#               nrow(cohort.vel.neg[which(cohort.vel.neg$sex=='F'),])/nrow(cohort.vel.neg),
#               nrow(cohort.vel.pos[which(cohort.vel.pos$sex=='F'),])/nrow(cohort.vel.pos),
#               nrow(cohort.vel.het[which(cohort.vel.het$sex=='F'),])/nrow(cohort.vel.het),
#               nrow(cohort.vel.het.control[which(cohort.vel.het.control$sex=='F'),])/nrow(cohort.vel.het.control),
#               nrow(cohort.remainder[which(cohort.remainder$sex=='F'),])/nrow(cohort.remainder)
#)

write.table(table1,file = 'table1.dbds.tsv', sep="\t", row.names = F, col.names = T, quote = F)

View(table1)

################### TABLE2 #####################

lab.cohort <- lab[which(lab$cpr_enc %in% c(vel.pos,vel.neg)),]
print(paste('mean sample/age diff:',mean(lab.cohort$age-lab.cohort$SAMPLINGAGE)))

table2<-c(
  c('NPU01566','cholesterol','>'),
  c('NPU01567','HDL','<'),
  c('NPU01568','LDL','>'),
  c('NPU01569','VLDL','>'), # 0.9
  #c('NPU03625','T3 (free)','<'),
  #c('NPU03624','T3','<'),
  #c('NPU03579','T4','<'),
  c('NPU03625/NPU03624/NPU03579','T3/T4-low','<'),
  c('NPU03625/NPU03624/NPU03579','T3/T4-high','>'),
  #c('NPU03577/NPU27547','TSH','>'),
  c('NPU03577/NPU27547','TSH-high','>'),
  c('NPU03577/NPU27547','TSH-low','<'),
  c('NPU27412','Glucose','>'),
  c('NPU02192','Glucose (plasma)','>'),
  c('NPU27300','HbA1c','>'),
  c('NPU03620','Triglycerides (fasting)','>'),
  c('NPU04094','Triglycerides','>'),
  c('NPU19651','ALT','>'),
  c('NPU19654','AST','>'),
  c('NPU19763','Ferritin','>'),
  c('NPU19748','hsCRP','>'),
  c('NPU03568','Platelets','<>'),
  c('NPU02636','Lymphocytes','<>'),
  c('NPU02593','Leukocytes','<>'),
  c('NPU02319','Hemoglobin','<>'),
  c('NPU26631','Myelocytes','>'),
  c('NPU01961','Erythrocytes','<>'), #ALSO NPU01960 ?!?
  c('NPU02902','ANC','<>')
)
table2<-data.frame(t(matrix(data=table2, nrow = 3, ncol = length(table2)/3 )), stringsAsFactors = F)
names(table2)<-c('NPU','TYPE','LIMIT')

test.anova <- c()
test.t <- c()
test.tM <- c()
test.tF <- c()
test.npos <- c()
test.nneg <- c()
test.abnormal.npos <- c()
test.abnormal.nneg <- c()
glm.p <- c()
glm.OR <- c()
glm.lower <- c()
glm.upper <- c()
cox.hazard <- c()
cox.p <- c()
cox.lower <- c()
cox.upper <- c()

for (i in seq(nrow(table2)))
{
  NPU <- table2$NPU[i]
  TYPE <- table2$TYPE[i]
  LIMIT <- table2$LIMIT[i]
  
  if (LIMIT == '>')
    lab.abnormal <- lab.cohort[which(lab.cohort$ANALYSISCODE %in% unlist(strsplit(NPU,'/')) & lab.cohort$VALUE > lab.cohort$REFERENCEINTERVAL_UPPERLIMIT),]
  if (LIMIT == '<')
    lab.abnormal <- lab.cohort[which(lab.cohort$ANALYSISCODE %in% unlist(strsplit(NPU,'/')) & lab.cohort$VALUE < lab.cohort$REFERENCEINTERVAL_LOWERLIMIT),]
  if (LIMIT == '<>')
    lab.abnormal <- lab.cohort[which(lab.cohort$ANALYSISCODE %in% unlist(strsplit(NPU,'/')) & ((lab.cohort$VALUE > lab.cohort$REFERENCEINTERVAL_UPPERLIMIT) | (lab.cohort$VALUE < lab.cohort$REFERENCEINTERVAL_LOWERLIMIT))),]
  
  ######### COX ########
  
  cox.normal <- lab.cohort[which(lab.cohort$ANALYSISCODE %in% unlist(strsplit(NPU,'/')) & !(lab.cohort$cpr_enc %in% lab.abnormal$cpr_enc) ),]
  cox.normal <- data.frame('cpr_enc'=unique(cox.normal$cpr_enc), 'time'=as.numeric((difftime(strptime(Sys.Date(), format = "%Y-%m-%d"),strptime('2008-01-01', format = "%Y-%m-%d"),units="days")/365.25)), stringsAsFactors = F)
  
  cox.abnormal <- data.frame('cpr_enc'=lab.abnormal$cpr_enc, 'time'=as.numeric((difftime(strptime(lab.abnormal$SAMPLINGDATE, format = "%Y-%m-%d"),strptime('2008-01-01', format = "%Y-%m-%d"),units="days")/365.25)), stringsAsFactors = F)
  cox.abnormal <- aggregate(time ~ cpr_enc, data = cox.abnormal, min)

  cox.test <- data.frame(cox.abnormal,'abnormal'='YES','VEL'=NA,stringsAsFactors = F)
  cox.test <- rbind(cox.test, data.frame(cox.normal,'abnormal'='NO','VEL'=NA,stringsAsFactors = F))

  cox.test[which(cox.test$cpr_enc %in% vel.pos), ]$VEL = '+'
  cox.test[which(cox.test$cpr_enc %in% vel.neg), ]$VEL = '-'

  cox.test <- merge(cox.test,cohort[,c(1,2,4)],by='cpr_enc')

  cox.test$VEL <- factor(cox.test$VEL, levels = c('+','-'))
  cox.test$abnormal <- factor(cox.test$abnormal)
  cox.test$sex <- factor(cox.test$sex)

  cox.res<-coxph(Surv(time,abnormal) ~ VEL + sex + age, data = cox.test, id = cox.test$cpr_enc)
  
  if (summary(cox.res)$conf.int[1,4] != Inf)
  {
    cox.hazard <- c(cox.hazard,summary(cox.res)$coefficients[1,2])
    cox.p <- c(cox.p,summary(cox.res)$coefficients[1,5])
    cox.lower <- c(cox.lower,summary(cox.res)$conf.int[1,3])
    cox.upper <- c(cox.upper,summary(cox.res)$conf.int[1,4])
  } else  {
    cox.hazard <- c(cox.hazard,NA)
    cox.p <- c(cox.p,NA)
    cox.lower <- c(cox.lower,NA)
    cox.upper <- c(cox.upper,NA)
  }
  
  ######### 2 way ANOVA ########
  lab.pos <- lab.cohort[which(lab.cohort$ANALYSISCODE %in% unlist(strsplit(NPU,'/')) & lab.cohort$cpr_enc %in% vel.pos),]
  lab.neg <- lab.cohort[which(lab.cohort$ANALYSISCODE %in% unlist(strsplit(NPU,'/')) & lab.cohort$cpr_enc %in% vel.neg),]
  lab.abnormal.pos <- lab.abnormal[which(lab.abnormal$cpr_enc %in% vel.pos),]
  lab.abnormal.neg <- lab.abnormal[which(lab.abnormal$cpr_enc %in% vel.neg),]
  
  test.npos <- c(test.npos, length(unique(lab.pos$cpr_enc)))
  test.nneg <- c(test.nneg, length(unique(lab.neg$cpr_enc)))
  test.abnormal.npos <- c(test.abnormal.npos, length(unique(lab.abnormal.pos$cpr_enc)))
  test.abnormal.nneg <- c(test.abnormal.nneg, length(unique(lab.abnormal.neg$cpr_enc)))

  if (length(unlist(strsplit(NPU,'/'))) == 1)
  {
    a<-data.frame('cpr_enc'=lab.pos$cpr_enc, 'VALUE'=lab.pos$VALUE, stringsAsFactors = F)
    b<-data.frame('cpr_enc'=lab.neg$cpr_enc, 'VALUE'=lab.neg$VALUE, stringsAsFactors = F)
    a <- a[which(!is.na(a$VALUE)),]
    b <- b[which(!is.na(b$VALUE)),]
    a <- aggregate(as.numeric(a$VALUE), list(a$cpr_enc), mean)
    b <- aggregate(as.numeric(b$VALUE), list(b$cpr_enc), mean)
    names(a)<-c('cpr_enc', 'VALUE')
    names(b)<-c('cpr_enc', 'VALUE')
    a<-merge(a,cohort[,c(1,2,4)],by='cpr_enc')
    b<-merge(b,cohort[,c(1,2,4)],by='cpr_enc')
    
    c <- rbind(
      data.frame(a, 'VEL'='+',stringsAsFactors = F),
      data.frame(b, 'VEL'='-',stringsAsFactors = F)
    )
    
    p.val<-summary(aov(VALUE ~ VEL + sex, data = c))[[1]]$'Pr(>F)'[1]
  }
  else
    p.val<-NA
  
  test.anova <- c(test.anova,p.val)

  ######### GLM ########
  lab.npu<-lab.cohort[which(lab.cohort$ANALYSISCODE %in% unlist(strsplit(NPU,'/'))),]
  a <- aggregate(as.numeric(lab.npu$VALUE), list(lab.npu$cpr_enc), mean)
  b <- aggregate(as.numeric(lab.npu$SAMPLINGAGE), list(lab.npu$cpr_enc), mean)
  
  if (nrow(c[which(c$VEL == '+'),]) > 1 & nrow(c[which(c$VEL == '-'),]) > 1 & length(unlist(strsplit(NPU,'/'))) == 1)
    test.t<-c(test.t,t.test(c[which(c$VEL == '+'),]$VALUE,c[which(c$VEL == '-'),]$VALUE)$p.value)
  else
    test.t<-c(test.t,NA)
  
  if (tail(test.abnormal.nneg,n=1) > 0)
  {
    names(a) <- c('cpr_enc','VALUE')
    names(b) <- c('cpr_enc','SAMPLINGAGE')
    
    lab.merged <- merge(a,b, by='cpr_enc')
    lab.merged<-merge(lab.merged,cohort[,c(1,2)],by='cpr_enc')
    
    lab.count <- data.frame(table('cpr_enc'=lab.npu$cpr_enc),stringsAsFactors = F)
    lab.merged <- merge(lab.merged,lab.count, by='cpr_enc')
    
    a<-data.frame('cpr_enc'=c(unique(lab.pos$cpr_enc),unique(lab.neg$cpr_enc)),'VEL'='+','abnormal'=F,stringsAsFactors = F)
    a[which(a$cpr_enc %in% vel.neg),]$VEL = '-'
    a[which(a$cpr_enc %in% unique(lab.abnormal$cpr_enc) ),]$abnormal = T
    
    lab.merged <- merge(lab.merged,a, by='cpr_enc')
    #summary(glm(abnormal ~ VEL + Freq + sex + SAMPLINGAGE,data = lab.merged, family = binomial))
    #glm.res<-glm(abnormal ~ VEL + Freq + sex + SAMPLINGAGE,data = lab.merged, family = binomial)
    #summary(glm(abnormal ~ VEL + sex + SAMPLINGAGE,data = lab.merged, family = binomial))
    lab.merged$VEL <- factor(lab.merged$VEL, levels = c('+','-'))
    lab.merged$abnormal <- factor(lab.merged$abnormal)
    lab.merged$sex <- factor(lab.merged$sex)
    
    glm.res<-glm(abnormal ~ VEL + sex + SAMPLINGAGE,data = lab.merged, family = binomial)
    glm.OR <- c(glm.OR,exp(summary(glm.res)$coefficients[2,1]))
    glm.p <- c(glm.p,summary(glm.res)$coefficients[2,4])
    glm.lower <- c(glm.lower,exp(confint((glm.res))[2,1]))
    glm.upper <- c(glm.upper,exp(confint((glm.res))[2,2]))
  
    #OR<-exp(summary(glm(abnormal ~ VEL + Freq + sex + SAMPLINGAGE,data = lab.merged, family = binomial))$coefficients[2,1])
    #p.val<-summary(glm(abnormal ~ VEL + Freq + sex + SAMPLINGAGE,data = lab.merged, family = binomial))$coefficients[2,4]
    #glm.OR<-c(glm.OR,OR)
    #glm.p<-c(glm.p,p.val)
  } else {
    glm.OR<-c(glm.OR,NA)
    glm.p<-c(glm.p,NA)
    glm.lower <- c(glm.lower,NA)
    glm.upper <- c(glm.upper,NA)
  }
}
table2<-data.frame(table2,'N.pos'=test.npos,'N.pos.abnormal'=test.abnormal.npos,'N.neg'=test.nneg,'N.neg.abnormal'=test.abnormal.nneg,'ANOVA'=test.anova, 't.test'=test.t,'GLM.OR'=glm.OR,'GLM.pval'=glm.p,'GLM.lower'=glm.lower,'GLM.upper'=glm.upper,'cox.hazard'=cox.hazard, 'cox.lower'=cox.lower, 'cox.upper'=cox.upper, 'cox.pval'=cox.p, stringsAsFactors = F)
write.table(table2,file = 'table2.dbds.tsv', sep="\t", row.names = F, col.names = T, quote = F)

table2<-table2[which(!is.na(table2$GLM.pval)),c(1,2,3,4,5,6,7,10,11,12,13,14,15)]
View(table2)

################### TABLE3 #####################

diagnosis <- read.table(file='lpr_psyk_rbc_gwas_chb.tsv', sep='\t', header = T, stringsAsFactors = FALSE)
cohort.diagnosis <- diagnosis[which(diagnosis$Source  == 'LPR' & diagnosis$cpr_enc %in% c(vel.neg,vel.pos) ),]

#ACT koder: H03AA prescription i tabel 4

#Thyroid - E05 (242) - E06 (244-245)

# BMI analyse

#ICD <- data.frame('CATEGORY'=c('Ischaemic heart disease','Thyroid','Diabetes','Obesity','lipid metabolisme'),
#                  'ICD8'=c('410,411,412,413,414','240,241,242,243,244,245,246','250,251,252,253,254,255,256,257,258','277','272'),
#                  'ICD10'=c('DI20,DI21,DI22,DI23,DI24,DI25','DE00,DE01,DE02,DE03,DE04,DE05,DE06,DE07','DE10,DE11,DE12,DE13,DE14','DE66','DE78'),
#                  stringsAsFactors = F)
ICD <- data.frame('CATEGORY'=c('Ischaemic heart disease','Thyroid1','Thyroid2','Diabetes','Obesity','lipid metabolisme'),
                  'ICD8'=c('410,411,412,413,414','242','244,245','250,251,252,253,254,255,256,257,258','277','272'),
                  'ICD10'=c('DI20,DI21,DI22,DI23,DI24,DI25','DE05','DE06','DE10,DE11,DE12,DE13,DE14','DE66','DE78'),
                  stringsAsFactors = F)

cox.npos <- c()
cox.nneg <- c()
cox.diag.npos <- c()
cox.diag.nneg <- c()
cox.hazard <- c()
cox.p <- c()
cox.lower <- c()
cox.upper <- c()
glm.coef <- c()
glm.p <-c()
glm.lower <- c()
glm.upper <- c()

for (j in seq(nrow(ICD)))
{
  ICD8 <- unlist(strsplit(ICD$ICD8[j],','))
  ICD10 <- unlist(strsplit(ICD$ICD10[j],','))
  
  cohort.data<-data.frame(cohort[which(cohort$cpr_enc %in% c(vel.neg,vel.pos)),c(1,2,3)], 'VEL'=NA,'DIAGNOSIS'=FALSE,'DIAGNOSIS.AGE'=NA,'DIAGNOSIS.YEARS'=NA, stringsAsFactors = F)
  
  cohort.data[which(cohort.data$cpr_enc %in% vel.pos),]$VEL = '+'
  cohort.data[which(cohort.data$cpr_enc %in% vel.neg),]$VEL = '-'
  
  for (i in seq(nrow(cohort.diagnosis)))
  {
    ID <- cohort.diagnosis$cpr_enc[i]

    #if (ID %in% vel.pos)
    #  cohort.data[which(cohort.data$cpr_enc == ID),]$VEL = '+'
    #if (ID %in% vel.neg)
    #  cohort.data[which(cohort.data$cpr_enc == ID),]$VEL = '-'

    version <- unlist(strsplit(cohort.diagnosis$Diagnosis[i],':'))[1]
    class <- unlist(strsplit(cohort.diagnosis$Diagnosis[i],':'))[2]
    
    if ( (version == 'ICD10' & substr(class,1,4) %in% ICD10) | (version == 'ICD8' & substr(class,1,3) %in% ICD8))
    {
      age = as.numeric(difftime(strptime(cohort.diagnosis$Date[i], format = "%Y-%m-%d"),strptime(cohort.data[which(cohort.data$cpr_enc == ID),]$birthdate, format = "%Y-%m-%d"),units="days")/365.25)
      years = as.numeric(difftime(strptime(cohort.diagnosis$Date[i], format = "%Y-%m-%d"),strptime('1977-01-01', format = "%Y-%m-%d"),units="days")/365.25)
      cohort.data[which(cohort.data$cpr_enc == ID),]$DIAGNOSIS = TRUE
      
      if (is.na(cohort.data[which(cohort.data$cpr_enc == ID),]$DIAGNOSIS.AGE) | (cohort.data[which(cohort.data$cpr_enc == ID),]$DIAGNOSIS.AGE > age))
      {
        cohort.data[which(cohort.data$cpr_enc == ID),]$DIAGNOSIS.AGE = age
        cohort.data[which(cohort.data$cpr_enc == ID),]$DIAGNOSIS.YEARS = years
      }
    }
  }
  cohort.data[which(cohort.data$DIAGNOSIS == FALSE),]$DIAGNOSIS.YEARS = as.numeric(difftime(strptime(Sys.Date(), format = "%Y-%m-%d"),strptime('1977-01-01', format = "%Y-%m-%d"),units="days")/365.25)
  cohort.data <- data.frame(cohort.data, 'age'=as.numeric(difftime(strptime(Sys.Date(), format = "%Y-%m-%d"),strptime(cohort.data$birthdate, format = "%Y-%m-%d"),units="days")/365.25), stringsAsFactors = F)

  cox.diag.npos <- c(cox.diag.npos, nrow(cohort.data[which(cohort.data$VEL == '+' & cohort.data$DIAGNOSIS == T),]))
  cox.diag.nneg <- c(cox.diag.nneg, nrow(cohort.data[which(cohort.data$VEL == '-' & cohort.data$DIAGNOSIS == T),]))

  cohort.data$sex = factor(cohort.data$sex)
  cohort.data$VEL = factor(cohort.data$VEL, levels = c('+','-'))
  cohort.data$I = factor(cohort.data$DIAGNOSIS)

  if (nrow(cohort.data[which(cohort.data$VEL == '+' & cohort.data$DIAGNOSIS == T),]) > 0 & nrow(cohort.data[which(cohort.data$VEL == '-' & cohort.data$DIAGNOSIS == T),]) > 0)
  {
    glm.res<-glm(DIAGNOSIS ~ VEL + age + sex, data = cohort.data, family=binomial)
    glm.coef <- c(glm.coef,exp(summary(glm.res)$coefficients[2,1]))
    glm.p <- c(glm.p,summary(glm.res)$coefficients[2,4])
    glm.lower <- c(glm.lower,exp(confint((glm.res))[2,1]))
    glm.upper <- c(glm.upper,exp(confint((glm.res))[2,2]))
  } else {
    glm.coef <- c(glm.coef,NA)
    glm.p <- c(glm.p,NA)
    glm.lower <- c(glm.lower,NA)
    glm.upper <- c(glm.upper,NA)
  }

  cox.res<-coxph(Surv(DIAGNOSIS.YEARS,DIAGNOSIS) ~ VEL + age + sex, data = cohort.data, id = cohort.data$cpr_enc)
  if (summary(cox.res)$conf.int[1,4] != Inf)
  {
    cox.hazard <- c(cox.hazard,summary(cox.res)$coefficients[1,2])
    cox.p <- c(cox.p,summary(cox.res)$coefficients[1,5])
    cox.lower <- c(cox.lower,summary(cox.res)$conf.int[1,3])
    cox.upper <- c(cox.upper,summary(cox.res)$conf.int[1,4])
  } else   {
    cox.hazard <- c(cox.hazard,NA)
    cox.p <- c(cox.p,NA)
    cox.lower <- c(cox.lower,NA)
    cox.upper <- c(cox.upper,NA)
  }
}

table3 <- data.frame('Category'=ICD$CATEGORY, 'N.diag.pos'=cox.diag.npos, 'N.diag.neg'=cox.diag.nneg, 'cox.hazard'=cox.hazard,'cox.CI95.lower'=cox.lower,'cox.CI95.upper'=cox.upper,'cox.p'=cox.p, 'glm.coef'=glm.coef,'glm.CI95.lower'=glm.lower,'cox.CI95.upper'=glm.upper,'cox.p'=glm.p, stringsAsFactors = F, check.names = F)
write.table(table3,file = 'table3.dbds.tsv', sep="\t", row.names = F, col.names = T, quote = F)

View(table3)

################### TABLE4 #####################

# stofskifte: H03A/H03B

statiner<-c('C10AA01','C10AA02','C10AA03','C10AA04','C10AA05','C10AA06','C10AA07','C10AA08')
antihypertensive <- c('C02A','C02B','C02C','C02D','C02K','C02L','C02N')
diabetes <- c('A10BA','A10BB','A10BC','A10BD','A10BF','A10BG','A10BH','A10BJ','A10BK','A10BX')

prescription <- read.table(file='prescriptions_chb_rbc_gwas.tsv', sep='\t', header = T, stringsAsFactors = FALSE)#, fill  = T)
cohort.prescription <- prescription[which(prescription$cpr_enc %in% c(vel.neg,vel.pos) ),]
cohort.data<-data.frame(cohort[which(cohort$cpr_enc %in% c(vel.neg,vel.pos)),c(1,2,3)], 'VEL'=NA,'statiner'=FALSE,'statiner.AGE'=NA,'statiner.YEARS'=NA,'antihypertensive'=FALSE,'antihypertensive.AGE'=NA,'antihypertensive.YEARS'=NA,'diabetes'=FALSE,'diabetes.AGE'=NA,'diabetes.YEARS'=NA, stringsAsFactors = F)

for (i in seq(nrow(cohort.prescription)))
{
  ID <- cohort.prescription$cpr_enc[i]
  
  age = as.numeric(difftime(strptime(cohort.prescription$expdato[i], format = "%d%b%Y"),strptime(cohort.data[which(cohort.data$cpr_enc == ID),]$birthdate, format = "%Y-%m-%d"),units="days")/365.25)
  years = as.numeric(difftime(strptime(cohort.prescription$expdato[i], format = "%d%b%Y"),strptime('2004-01-01', format = "%Y-%m-%d"),units="days")/365.25)
  
  if ( cohort.prescription$atc_kode[i] %in% statiner)
  {
    cohort.data[which(cohort.data$cpr_enc == ID),]$statiner = TRUE
    
    if (is.na(cohort.data[which(cohort.data$cpr_enc == ID),]$statiner.AGE) | (cohort.data[which(cohort.data$cpr_enc == ID),]$statiner.AGE > age))
    {
      cohort.data[which(cohort.data$cpr_enc == ID),]$statiner.AGE = age
      cohort.data[which(cohort.data$cpr_enc == ID),]$statiner.YEARS = years
    }
  }
  
  if (substr(cohort.prescription$atc_kode[i],1,4) %in% antihypertensive)
  {
    cohort.data[which(cohort.data$cpr_enc == ID),]$antihypertensive = TRUE
    
    if (is.na(cohort.data[which(cohort.data$cpr_enc == ID),]$antihypertensive.AGE) | (cohort.data[which(cohort.data$cpr_enc == ID),]$antihypertensive.AGE > age))
    {
      cohort.data[which(cohort.data$cpr_enc == ID),]$antihypertensive.AGE = age
      cohort.data[which(cohort.data$cpr_enc == ID),]$antihypertensive.YEARS = years
    }
  }

  if ( substr(cohort.prescription$atc_kode[i],1,5) %in% diabetes)
  {
    cohort.data[which(cohort.data$cpr_enc == ID),]$diabetes = TRUE
    
    if (is.na(cohort.data[which(cohort.data$cpr_enc == ID),]$diabetes.AGE) | (cohort.data[which(cohort.data$cpr_enc == ID),]$diabetes.AGE > age))
    {
      cohort.data[which(cohort.data$cpr_enc == ID),]$diabetes.AGE = age
      cohort.data[which(cohort.data$cpr_enc == ID),]$diabetes.YEARS = years
    }
  }
}

cohort.data[which(cohort.data$cpr_enc %in% vel.pos),]$VEL = '+'
cohort.data[which(cohort.data$cpr_enc %in% vel.neg),]$VEL = '-'

cohort.data[which(cohort.data$statiner == FALSE),]$statiner.YEARS = as.numeric(difftime(strptime(Sys.Date(), format = "%Y-%m-%d"),strptime('2004-01-01', format = "%Y-%m-%d"),units="days")/365.25)
cohort.data[which(cohort.data$antihypertensive == FALSE),]$antihypertensive.YEARS = as.numeric(difftime(strptime(Sys.Date(), format = "%Y-%m-%d"),strptime('2004-01-01', format = "%Y-%m-%d"),units="days")/365.25)
cohort.data[which(cohort.data$diabetes == FALSE),]$diabetes.YEARS = as.numeric(difftime(strptime(Sys.Date(), format = "%Y-%m-%d"),strptime('2004-01-01', format = "%Y-%m-%d"),units="days")/365.25)

cohort.data <- data.frame(cohort.data, 'age'=as.numeric(difftime(strptime(Sys.Date(), format = "%Y-%m-%d"),strptime(cohort.data$birthdate, format = "%Y-%m-%d"),units="days")/365.25), stringsAsFactors = F)

cox.prescription.npos <- c()
cox.prescription.nneg <- c()
cox.prescription.npos <- c(cox.prescription.npos,nrow(cohort.data[which(cohort.data$VEL == '+' & cohort.data$statiner == TRUE),]))
cox.prescription.nneg <- c(cox.prescription.nneg,nrow(cohort.data[which(cohort.data$VEL == '-' & cohort.data$statiner == TRUE),]))
cox.prescription.npos <- c(cox.prescription.npos,nrow(cohort.data[which(cohort.data$VEL == '+' & cohort.data$antihypertensive == TRUE),]))
cox.prescription.nneg <- c(cox.prescription.nneg,nrow(cohort.data[which(cohort.data$VEL == '-' & cohort.data$antihypertensive == TRUE),]))
cox.prescription.npos <- c(cox.prescription.npos,nrow(cohort.data[which(cohort.data$VEL == '+' & cohort.data$diabetes == TRUE),]))
cox.prescription.nneg <- c(cox.prescription.nneg,nrow(cohort.data[which(cohort.data$VEL == '-' & cohort.data$diabetes == TRUE),]))

cohort.data$sex = factor(cohort.data$sex)
cohort.data$VEL = factor(cohort.data$VEL, levels = c('+','-'))
cohort.data$statiner = factor(cohort.data$statiner)
cohort.data$antihypertensive = factor(cohort.data$antihypertensive)
cohort.data$diabetes = factor(cohort.data$diabetes)
cox.hazard <- c()
cox.p <- c()
cox.lower <- c()
cox.upper <- c()

glm.coef <- c()
glm.p <-c()
glm.lower <- c()
glm.upper <- c()

if (cox.prescription.npos[1] > 0 & cox.prescription.nneg[1] > 0)
{
  cox.res<-coxph(Surv(statiner.YEARS,statiner) ~ VEL + age + sex, data = cohort.data, id = cohort.data$cpr_enc)
  cox.hazard <- c(cox.hazard,summary(cox.res)$coefficients[1,2])
  cox.p <- c(cox.p,summary(cox.res)$coefficients[1,5])
  cox.lower <- c(cox.lower,summary(cox.res)$conf.int[1,3])
  cox.upper <- c(cox.upper,summary(cox.res)$conf.int[1,4])

  glm.res<-glm(statiner ~ VEL + sex + age ,data = cohort.data, family = binomial)
  glm.coef <- c(glm.coef,exp(summary(glm.res)$coefficients[2,1]))
  glm.p <- c(glm.p,summary(glm.res)$coefficients[2,4])
  glm.lower <- c(glm.lower,exp(confint((glm.res))[2,1]))
  glm.upper <- c(glm.upper,exp(confint((glm.res))[2,2]))
} else {
  cox.hazard <- c(cox.hazard,NA)
  cox.p <- c(cox.p,NA)
  cox.lower <- c(cox.lower,NA)
  cox.upper <- c(cox.upper,NA)
  
  glm.coef <- c(glm.coef,NA)
  glm.p <-c(glm.p,NA)
  glm.lower <- c(glm.lower,NA)
  glm.upper <- c(glm.upper,NA)
}

if (cox.prescription.npos[2] > 0 & cox.prescription.nneg[2] > 0)
{
  cox.res<-coxph(Surv(antihypertensive.YEARS,antihypertensive) ~ VEL + age + sex, data = cohort.data, id = cohort.data$cpr_enc)
  cox.hazard <- c(cox.hazard,summary(cox.res)$coefficients[1,2])
  cox.p <- c(cox.p,summary(cox.res)$coefficients[1,5])
  cox.lower <- c(cox.lower,summary(cox.res)$conf.int[1,3])
  cox.upper <- c(cox.upper,summary(cox.res)$conf.int[1,4])
  
  glm.res<-glm(antihypertensive ~ VEL + sex + age ,data = cohort.data, family = binomial)
  glm.coef <- c(glm.coef,exp(summary(glm.res)$coefficients[2,1]))
  glm.p <- c(glm.p,summary(glm.res)$coefficients[2,4])
  glm.lower <- c(glm.lower,exp(confint((glm.res))[2,1]))
  glm.upper <- c(glm.upper,exp(confint((glm.res))[2,2]))
} else {
  cox.hazard <- c(cox.hazard,NA)
  cox.p <- c(cox.p,NA)
  cox.lower <- c(cox.lower,NA)
  cox.upper <- c(cox.upper,NA)
  
  glm.coef <- c(glm.coef,NA)
  glm.p <-c(glm.p,NA)
  glm.lower <- c(glm.lower,NA)
  glm.upper <- c(glm.upper,NA)
}

if (cox.prescription.npos[3] > 0 & cox.prescription.nneg[3] > 0)
{
  cox.res<-coxph(Surv(diabetes.YEARS,diabetes) ~ VEL + age + sex, data = cohort.data, id = cohort.data$cpr_enc)
  cox.hazard <- c(cox.hazard,summary(cox.res)$coefficients[1,2])
  cox.p <- c(cox.p,summary(cox.res)$coefficients[1,5])
  cox.lower <- c(cox.lower,summary(cox.res)$conf.int[1,3])
  cox.upper <- c(cox.upper,summary(cox.res)$conf.int[1,4])
  
  glm.res<-glm(diabetes ~ VEL + sex + age ,data = cohort.data, family = binomial)
  glm.coef <- c(glm.coef,exp(summary(glm.res)$coefficients[2,1]))
  glm.p <- c(glm.p,summary(glm.res)$coefficients[2,4])
  glm.lower <- c(glm.lower,exp(confint((glm.res))[2,1]))
  glm.upper <- c(glm.upper,exp(confint((glm.res))[2,2]))
} else {
  cox.hazard <- c(cox.hazard,NA)
  cox.p <- c(cox.p,NA)
  cox.lower <- c(cox.lower,NA)
  cox.upper <- c(cox.upper,NA)
  
  glm.coef <- c(glm.coef,NA)
  glm.p <-c(glm.p,NA)
  glm.lower <- c(glm.lower,NA)
  glm.upper <- c(glm.upper,NA)
}

table4<-data.frame('prescriptions'=c('statines','antihypertensives','diabetes'), 'N.prescription.pos'=cox.prescription.npos, 'N.prescription.neg'=cox.prescription.nneg, 'cox.hazard'=cox.hazard,'cox.CI95.lower'=cox.lower,'cox.CI95.upper'=cox.upper,'cox.p'=cox.p,'glm.coef'=glm.coef,'glm.CI95.lower'=glm.lower,'glm.CI95.upper'=glm.upper,'glm.p'=glm.p, stringsAsFactors = F )
write.table(table4,file = 'table4.dbds.tsv', sep="\t", row.names = F, col.names = T, quote = F)

View(table4)

################### TABLE5 #####################

dbds1 <- read.table(file='dbds1_rbc_chb_20201021.tsv', sep='\t', header = T, stringsAsFactors = FALSE)
dbds2 <- read.table(file='dbds2_rbc_chb_20201021.tsv', sep='\t', header = T, stringsAsFactors = FALSE)
dbds3 <- read.table(file='dbds3_rbc_chb_20201021.tsv', sep='\t', header = T, stringsAsFactors = FALSE)

dbds1 <- dbds1[which(!is.na(dbds1$date) & !is.na(dbds1$height) & !is.na(dbds1$weight)),]
dbds2 <- dbds2[which(!is.na(dbds2$date) & !is.na(dbds2$height) & !is.na(dbds2$weight)),]
dbds3 <- dbds3[which(!is.na(dbds3$date) & !is.na(dbds3$height) & !is.na(dbds3$weight) & dbds3$weight != 'NaN' & dbds3$height != 'NaN'),]

dbds3[which(dbds3$weight == '>150'),]$weight = '150'
dbds3$weight <- as.integer(dbds3$weight)

dbds1 <- merge(dbds1[,c(1,3,4,5)],cohort[,c(1,2,3)],by='cpr_enc')
dbds1.bmi <- data.frame('cpr_enc'=dbds1$cpr_enc,sex=dbds1$sex,'age'= as.integer(floor(difftime(as.Date(dbds1$date),as.Date(dbds1$birthdate),units = 'days')/365.25)),'bmi'=dbds1$weight/((dbds1$height/100)^2),'weight'=dbds1$weight)
dbds2 <- merge(dbds2[,c(1,3,4,5)],cohort[,c(1,2,3)],by='cpr_enc')
dbds2.bmi <- data.frame('cpr_enc'=dbds2$cpr_enc,sex=dbds2$sex,'age'= as.integer(floor(difftime(as.Date(dbds2$date),as.Date(dbds2$birthdate),units = 'days')/365.25)),'bmi'=dbds2$weight/((dbds2$height/100)^2),'weight'=dbds2$weight)
dbds3 <- merge(dbds3[,c(1,3,4,5)],cohort[,c(1,2,3)],by='cpr_enc')
dbds3.bmi <- data.frame('cpr_enc'=dbds3$cpr_enc,sex=dbds3$sex,'age'= as.integer(floor(difftime(as.Date(dbds3$date),as.Date(dbds3$birthdate),units = 'days')/365.25)),'bmi'=dbds3$weight/((dbds3$height/100)^2),'weight'=dbds3$weight)

dbds1.bmi <- dbds1.bmi[ which(!dbds1.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds2.bmi$cpr_enc) & !dbds1.bmi$cpr_enc %in% intersect(dbds1.bmi$cpr_enc,dbds3.bmi$cpr_enc) ),]
dbds2.bmi <- dbds2.bmi[ which(!dbds2.bmi$cpr_enc %in% intersect(dbds2.bmi$cpr_enc,dbds3.bmi$cpr_enc)),]
dbds.bmi <- rbind(dbds1.bmi,dbds2.bmi,dbds3.bmi)

nrow(dbds.bmi[dbds.bmi$cpr_enc %in% vel.pos,])
nrow(dbds.bmi[dbds.bmi$cpr_enc %in% vel.neg,])

dbds.bmi <- data.frame(dbds.bmi, 'VEL'=NA, 'Normal'=FALSE, 'Overweight'=FALSE,'Obese'=FALSE, stringsAsFactors = F)

dbds.bmi<-dbds.bmi[which(dbds.bmi$cpr_enc %in% c(vel.pos,vel.neg)),]
dbds.bmi[which(dbds.bmi$cpr_enc %in% vel.pos),]$VEL = '+'
dbds.bmi[which(dbds.bmi$cpr_enc %in% vel.neg),]$VEL = '-'

dbds.bmi$sex <- factor(dbds.bmi$sex)
dbds.bmi$VEL <- factor(dbds.bmi$VEL, levels = c('+','-'))

dbds.bmi[which(dbds.bmi$bmi < 25),]$Normal = TRUE
dbds.bmi[which(dbds.bmi$bmi >= 25),]$Overweight = TRUE
dbds.bmi[which(dbds.bmi$bmi >= 30),]$Obese = TRUE

summary(aov(Normal ~ VEL + age + sex, data = dbds.bmi))
summary(aov(Overweight ~ VEL + age + sex, data = dbds.bmi))
summary(aov(Obese ~ VEL + age + sex, data = dbds.bmi))

if (TRUE)
{
  # sex-stratified weight based
  par(mfrow=c(1,1))
  
  wilcox.res <- wilcox.test(weight~VEL, data = dbds.bmi[which(dbds.bmi$sex == 'M'),])
  wilcox.p<-wilcox.res[3]$p.value
  glm.res <- glm(Normal ~ VEL + age,data = dbds.bmi[which(dbds.bmi$sex == 'M'),],family = binomial)
  glm.coef <- exp(summary(glm.res)$coefficients[2,1])
  glm.p <- summary(glm.res)$coefficients[2,4]
  glm.lower <- exp(confint((glm.res))[2,1])
  glm.upper <- exp(confint((glm.res))[2,2])
  table5.M<-data.frame('WEIGHT'=c('MEAN (wilcox) MEN','OVERWEIGHT (glm) MEN'), 'VEL.pos'=c(mean(dbds.bmi[which(dbds.bmi$VEL == '+' & dbds.bmi$sex == 'M'),]$weight),as.integer(nrow(dbds.bmi[which(dbds.bmi$VEL == '+'  & dbds.bmi$sex == 'M' & dbds.bmi$Overweight == TRUE),]))), 'VEL.neg'=c(mean(dbds.bmi[which(dbds.bmi$VEL == '-' & dbds.bmi$sex == 'M'),]$weight),as.integer(nrow(dbds.bmi[which(dbds.bmi$VEL == '-' & dbds.bmi$sex == 'M' & dbds.bmi$Overweight == TRUE),]))),'p'=c(wilcox.p,glm.p),'coef'=c(NA,glm.coef),'CI95.lower'=c(NA,glm.lower),'CI95.upper'=c(NA,glm.upper), stringsAsFactors = F )
  pdf('table5.boxplot.M.pdf')
  boxplot(dbds.bmi[which(dbds.bmi$cpr_enc %in% vel.pos  & dbds.bmi$sex == 'M'),]$weight,dbds.bmi[which(dbds.bmi$cpr_enc %in% vel.neg & dbds.bmi$sex == 'M'),]$weight, names = c("VEL.POS", "VEL.NEG"), ylim=c(50,150), main="DBDS COHORT WEIGHT: MEN")
  legend(x="topright", legend = c(paste("POS:",nrow(dbds.bmi[which(dbds.bmi$VEL == '+' & dbds.bmi$sex == 'M'),])),paste("NEG:",nrow(dbds.bmi[which(dbds.bmi$VEL == '-' & dbds.bmi$sex == 'M'),])), paste("P =",round(p.val.wilcox,2))))
  dev.off()
  wilcox.res <- wilcox.test(weight~VEL, data = dbds.bmi[which(dbds.bmi$sex == 'F'),])
  wilcox.p<-wilcox.res[3]$p.value
  glm.res <- glm(Normal ~ VEL + age,data = dbds.bmi[which(dbds.bmi$sex == 'F'),],family = binomial)
  glm.coef <- exp(summary(glm.res)$coefficients[2,1])
  glm.p <- summary(glm.res)$coefficients[2,4]
  glm.lower <- exp(confint((glm.res))[2,1])
  glm.upper <- exp(confint((glm.res))[2,2])
  table5.F<-data.frame('WEIGHT'=c('MEAN (wilcox) WOMEN','OVERWEIGHT (glm) WOMEN'), 'VEL.pos'=c(mean(dbds.bmi[which(dbds.bmi$VEL == '+' & dbds.bmi$sex == 'F'),]$weight),as.integer(nrow(dbds.bmi[which(dbds.bmi$VEL == '+'  & dbds.bmi$sex == 'F' & dbds.bmi$Overweight == TRUE),]))), 'VEL.neg'=c(mean(dbds.bmi[which(dbds.bmi$VEL == '-'  & dbds.bmi$sex == 'F'),]$weight),as.integer(nrow(dbds.bmi[which(dbds.bmi$VEL == '-' & dbds.bmi$sex == 'F' & dbds.bmi$Overweight == TRUE),]))),'p'=c(wilcox.p,glm.p),'coef'=c(NA,glm.coef),'CI95.lower'=c(NA,glm.lower),'CI95.upper'=c(NA,glm.upper), stringsAsFactors = F )
  pdf('table5.boxplot.F.pdf')
  boxplot(dbds.bmi[which(dbds.bmi$cpr_enc %in% vel.pos  & dbds.bmi$sex == 'F'),]$weight,dbds.bmi[which(dbds.bmi$cpr_enc %in% vel.neg & dbds.bmi$sex == 'F'),]$weight, names = c("VEL.POS", "VEL.NEG"), ylim=c(50,150), main="DBDS COHORT WEIGHT: WOMEN")
  legend(x="topright", legend = c(paste("POS:",nrow(dbds.bmi[which(dbds.bmi$VEL == '+' & dbds.bmi$sex == 'F'),])),paste("NEG:",nrow(dbds.bmi[which(dbds.bmi$VEL == '-' & dbds.bmi$sex == 'F'),])), paste("P =",round(p.val.wilcox,2))))
  dev.off()
  
  table5<-rbind(table5.M,table5.F)
  #write.table(table5,file = 'table5.dbds.tsv', sep="\t", row.names = F, col.names = T, quote = F)
  #View(table5)
  
  
  
} else {
  # non-stratified bmi based
  wilcox.res <- wilcox.test(bmi~VEL, data = dbds.bmi)
  wilcox.p<-wilcox.res[3]$p.value
  
  glm.res <- glm(Normal ~ VEL + age + factor(sex),data = dbds.bmi,family = binomial)
  glm.coef <- exp(summary(glm.res)$coefficients[2,1])
  glm.p <- summary(glm.res)$coefficients[2,4]
  glm.lower <- exp(confint((glm.res))[2,1])
  glm.upper <- exp(confint((glm.res))[2,2])
  
  table5<-data.frame('BMI'=c('MEAN (wilcox)','OVERWEIGHT (glm)'), 'VEL.pos'=c(mean(dbds.bmi[which(dbds.bmi$VEL == '+'),]$bmi),as.integer(nrow(dbds.bmi[which(dbds.bmi$VEL == '+' & dbds.bmi$Overweight == TRUE),]))), 'VEL.neg'=c(mean(dbds.bmi[which(dbds.bmi$VEL == '-'),]$bmi),as.integer(nrow(dbds.bmi[which(dbds.bmi$VEL == '-' & dbds.bmi$Overweight == TRUE),]))),'p'=c(wilcox.p,glm.p),'coef'=c(NA,glm.coef),'CI95.lower'=c(NA,glm.lower),'CI95.upper'=c(NA,glm.upper), stringsAsFactors = F )
  write.table(table5,file = 'table5.dbds.tsv', sep="\t", row.names = F, col.names = T, quote = F)
  View(table5)
  
  boxplot(dbds.bmi[which(dbds.bmi$cpr_enc %in% vel.pos),]$bmi,dbds.bmi[which(dbds.bmi$cpr_enc %in% vel.neg),]$bmi, names = c("VEL.POS", "VEL.NEG"), main="DBDS")
  legend(x="topright", legend = c(paste("POS:",length(vel.pos)),paste("NEG:",length(vel.neg)), paste("P =",round(p.val.wilcox,2))))
} 

################### TABLE6 #####################

#dbds1 <- read.table(file='dbds1_rbc_chb_20201021.tsv', sep='\t', header = T, stringsAsFactors = FALSE)

#dbds1$

#################### END ####################
View(read.table(file = 'table2.dbds.tsv', sep="\t", header=T)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),c(1,2,3,4,5,6,7,10,11,12,13,14,15)])
View(read.table(file = 'table3.dbds.tsv', sep="\t", header=T))
View(read.table(file = 'table4.dbds.tsv', sep="\t", header=T))
###

lab.npu<-lab.cohort[which(lab.cohort$ANALYSISCODE %in% unlist(strsplit(NPU,'/'))),]
a <- aggregate(as.numeric(lab.npu$VALUE), list(lab.npu$cpr_enc), mean)
b <- aggregate(as.numeric(lab.npu$SAMPLINGAGE), list(lab.npu$cpr_enc), mean)


summary(glm(abnormal ~ VEL + Freq + sex,data = lab.merged,family = binomial))
p.val<-summary(glm(abnormal ~ VEL + Freq + sex,data = lab.merged,family = binomial))[12]$coefficients[2,4]

boxplot(dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.pos & dbds12.bmi$age %in% ages,]$bmi,dbds12.bmi[dbds12.bmi$cpr_enc %in% vel.neg,]$bmi, names = c("VEL.POS", "VEL.NEG"), main="Imputed")
legend(x="topright", legend = c("POS: 40739","NEG: 37", "P = 0.02944"))

#