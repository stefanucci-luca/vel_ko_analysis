library("survival")

sexes <- c('M','F')

for (sex in sexes)
{
{
  reset <- FALSE
  
  cohort <- read.table(file='cohort_overview_20200707.tsv', sep='\t', header = T, stringsAsFactors = FALSE)#, fill  = T)
  cohort <- cohort[which(cohort$sex == sex),]
  cohort.sex <- unique(cohort[which(cohort$sex == sex),]$cpr_enc)
  cohort<-data.frame('cpr_enc'=cohort$cpr_enc, 'sex'=cohort$sex, 'birthdate'=cohort$birthdate, 'age'= as.numeric(floor(difftime(strptime(Sys.Date(), format = "%Y-%m-%d"),strptime(cohort$birthdate, format = "%Y-%m-%d"),units="days")/365.25)), 'cohort'=cohort$cohort, stringsAsFactors = FALSE)
  
  load('vel.neg.combo.Rda')
  load('vel.pos.combo.Rda')
  
  vel.neg <- vel.neg[vel.neg %in% cohort.sex]
  vel.pos <- vel.pos[vel.pos %in% cohort.sex]

  vel.neg.chb <- vel.neg[vel.neg %in% cohort[which(cohort$cohort == 'CHB'),]$cpr_enc]
  vel.neg.dbds <- vel.neg[vel.neg %in% cohort[which(cohort$cohort != 'CHB'),]$cpr_enc]
  vel.pos.chb <- vel.pos[vel.pos %in% cohort[which(cohort$cohort == 'CHB'),]$cpr_enc]
  vel.pos.dbds <- vel.pos[vel.pos %in% cohort[which(cohort$cohort != 'CHB'),]$cpr_enc]
}

################### TABLE1 #####################

table1 <- data.frame(matrix(data = NA, nrow = 3, ncol = 8), stringsAsFactors = F)
names(table1)<-c('TYPE','VEL-neg','VEL-pos','VEL-neg-DBDS','VEL-pos-DBDS','VEL-neg-CHB','VEL-pos-CHB','remainder')
cohort.vel.neg <- cohort[which(cohort$cpr_enc %in% vel.neg),]
cohort.vel.pos <- cohort[which(cohort$cpr_enc %in% vel.pos),]
cohort.remainder <- cohort[which(!cohort$cpr_enc %in% c(vel.neg,vel.pos)),]
table1[1,] = c('size',nrow(cohort.vel.neg),nrow(cohort.vel.pos),nrow(cohort.vel.neg[which(cohort.vel.neg$cohort == 'DBDS'),]),nrow(cohort.vel.pos[which(cohort.vel.pos$cohort == 'DBDS'),]),nrow(cohort.vel.neg[which(cohort.vel.neg$cohort == 'CHB'),]),nrow(cohort.vel.pos[which(cohort.vel.pos$cohort == 'CHB'),]),nrow(cohort.remainder))
table1[2,] = c('mean-age',mean(cohort.vel.neg$age),mean(cohort.vel.pos$age),mean(cohort.vel.neg[which(cohort.vel.neg$cohort == 'DBDS'),]$age),mean(cohort.vel.pos[which(cohort.vel.pos$cohort == 'DBDS'),]$age),mean(cohort.vel.neg[which(cohort.vel.neg$cohort == 'CHB'),]$age),mean(cohort.vel.pos[which(cohort.vel.pos$cohort == 'CHB'),]$age),mean(cohort.remainder$age))
table1[3,] = c('sex-fraction (F)',
               nrow(cohort.vel.neg[which(cohort.vel.neg$sex=='F'),])/nrow(cohort.vel.neg),
               nrow(cohort.vel.pos[which(cohort.vel.pos$sex=='F'),])/nrow(cohort.vel.pos),
               nrow(cohort.vel.neg[which(cohort.vel.neg$sex=='F' & cohort.vel.neg$cohort == 'DBDS'),])/nrow(cohort.vel.neg[which(cohort.vel.neg$cohort == 'DBDS'),]),
               nrow(cohort.vel.pos[which(cohort.vel.pos$sex=='F' & cohort.vel.pos$cohort == 'DBDS'),])/nrow(cohort.vel.pos[which(cohort.vel.pos$cohort == 'DBDS'),]),
               nrow(cohort.vel.neg[which(cohort.vel.neg$sex=='F' & cohort.vel.neg$cohort == 'CHB'),])/nrow(cohort.vel.neg[which(cohort.vel.neg$cohort == 'CHB'),]),
               nrow(cohort.vel.pos[which(cohort.vel.pos$sex=='F' & cohort.vel.pos$cohort == 'CHB'),])/nrow(cohort.vel.pos[which(cohort.vel.pos$cohort == 'CHB'),]),
               nrow(cohort.remainder[which(cohort.remainder$sex=='F'),])/nrow(cohort.remainder)
)

write.table(table1,file = paste('table1.combo.',sex,'.tsv',sep = ''), sep="\t", row.names = F, col.names = T, quote = F)

assign(paste('table1.',sex,sep=''),table1)
#View(table1)

################### TABLE2 #####################

lab.cohort <- lab[which(lab$cpr_enc %in% c(vel.pos,vel.neg)),]
#print(paste('mean sample/age diff:',mean(lab.cohort$age-lab.cohort$SAMPLINGAGE)))

{
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
  cox.hazard <- c()
  cox.p <- c()
  cox.lower <- c()
  cox.upper <- c()
}

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

  cox.test <- data.frame(cox.abnormal,'abnormal'='YES','VEL'=NA,'cohort'=NA,stringsAsFactors = F)
  cox.test <- rbind(cox.test, data.frame(cox.normal,'abnormal'='NO','VEL'=NA,'cohort'=NA,stringsAsFactors = F))

  if (sum(cox.test$cpr_enc %in% c(vel.pos.dbds,vel.neg.dbds))>0 & sum(cox.test$cpr_enc %in% c(vel.pos.chb,vel.neg.chb))>0)
  {
    cox.test[which(cox.test$cpr_enc %in% vel.pos), ]$VEL = '+'
    cox.test[which(cox.test$cpr_enc %in% vel.neg), ]$VEL = '-'
    if (sum(cox.test$cpr_enc %in% c(vel.pos.dbds,vel.neg.dbds))>0)
      cox.test[which(cox.test$cpr_enc %in% c(vel.pos.dbds,vel.neg.dbds)), ]$cohort = 'DBDS'
    if (sum(cox.test$cpr_enc %in% c(vel.pos.chb,vel.neg.chb))>0)
      cox.test[which(cox.test$cpr_enc %in% c(vel.pos.chb,vel.neg.chb)), ]$cohort = 'CHB'

    cox.test <- merge(cox.test,cohort[,c(1,2,4)],by='cpr_enc')
    
    cox.test$VEL <- factor(cox.test$VEL, levels = c('+','-'))
    cox.test$abnormal <- factor(cox.test$abnormal)
    cox.test$sex <- factor(cox.test$sex)
    cox.test$cohort <- factor(cox.test$cohort)
    
    cox.res<-coxph(Surv(time,abnormal) ~ VEL + age + cohort, data = cox.test, id = cox.test$cpr_enc)
    
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
  } else {
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
    
    p.val<-summary(aov(VALUE ~ VEL, data = c))[[1]]$'Pr(>F)'[1]
  } else {
    p.val<-NA
  }
  test.anova <- c(test.anova,p.val)

  ######### GLM ########
  lab.npu<-lab.cohort[which(lab.cohort$ANALYSISCODE %in% unlist(strsplit(NPU,'/'))),]
  a <- aggregate(as.numeric(lab.npu$VALUE), list(lab.npu$cpr_enc), mean)
  b <- aggregate(as.numeric(lab.npu$SAMPLINGAGE), list(lab.npu$cpr_enc), mean)

  if (nrow(c[which(c$VEL == '+'),]) > 1 & nrow(c[which(c$VEL == '-'),]) > 1 & length(unlist(strsplit(NPU,'/'))) == 1)
  {
    test.t<-c(test.t,t.test(c[which(c$VEL == '+'),]$VALUE,c[which(c$VEL == '-'),]$VALUE)$p.value)
  } else {
    test.t<-c(test.t,NA)
  }
  
  if (tail(test.abnormal.nneg,n=1) > 0 & length(unique(lab.npu$cohort)) > 1)
  {
    names(a) <- c('cpr_enc','VALUE')
    names(b) <- c('cpr_enc','SAMPLINGAGE')

    lab.merged <- merge(a,b, by='cpr_enc')
    lab.merged<-merge(lab.merged,cohort[,c(1,2)],by='cpr_enc')

    lab.count <- data.frame(table('cpr_enc'=lab.npu$cpr_enc),stringsAsFactors = F)
    lab.merged <- merge(lab.merged,lab.count, by='cpr_enc')

    a<-data.frame('cpr_enc'=c(unique(lab.pos$cpr_enc),unique(lab.neg$cpr_enc)),'VEL'='+','abnormal'=F,'cohort'='DBDS',stringsAsFactors = F)
    a[which(a$cpr_enc %in% vel.neg),]$VEL = '-'
    a[which(a$cpr_enc %in% unique(lab.abnormal$cpr_enc) ),]$abnormal = T
    a[which(a$cpr_enc %in% c(vel.neg.chb,vel.pos.chb)),]$cohort = 'CHB'

    lab.merged <- merge(lab.merged,a, by='cpr_enc')
    
    lab.merged$sex <- factor(lab.merged$sex)
    lab.merged$VEL <- factor(lab.merged$VEL, levels = c('+','-'))
    lab.merged$abnormal <- factor(lab.merged$abnormal)
    lab.merged$cohort <- factor(lab.merged$cohort)
    
    #summary(glm(abnormal ~ VEL + Freq + sex + SAMPLINGAGE + cohort,data = lab.merged,family = binomial))
    OR<-exp(summary(glm(abnormal ~ VEL + Freq  + SAMPLINGAGE + cohort,data = lab.merged, family = binomial))$coefficients[2,1])
    p.val<-summary(glm(abnormal ~ VEL + Freq  + SAMPLINGAGE + cohort,data = lab.merged, family = binomial))$coefficients[2,4]
    glm.OR<-c(glm.OR,OR)
    glm.p<-c(glm.p,p.val)
  } else {
    glm.OR<-c(glm.OR,NA)
    glm.p<-c(glm.p,NA)
  }
}
table2<-data.frame(table2,'N.pos'=test.npos,'N.pos.abnormal'=test.abnormal.npos,'N.neg'=test.nneg,'N.neg.abnormal'=test.abnormal.nneg,'ANOVA'=test.anova, 't.test'=test.t,'GLM.OR'=glm.OR,'GLM.pval'=glm.p, 'cox.hazard'=cox.hazard, 'cox.lower'=cox.lower, 'cox.upper'=cox.upper, 'cox.pval'=cox.p, stringsAsFactors = F)
write.table(table2,file = paste('table2.combo.',sex,'.tsv',sep = ''), sep="\t", row.names = F, col.names = T, quote = F)

table2<-table2[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),c(1,2,3,4,5,6,7,10,11,12,13,14,15)]

assign(paste('table2.',sex,sep=''),table2)
#View(table2)

################### TABLE3 #####################

diagnosis <- read.table(file='lpr_psyk_rbc_gwas_chb.tsv', sep='\t', header = T, stringsAsFactors = FALSE)
cohort.diagnosis <- diagnosis[which(diagnosis$Source  == 'LPR' & diagnosis$DiagnosisType %in% c('A','B','G') & diagnosis$cpr_enc %in% c(vel.neg,vel.pos)),]

{
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
}

for (j in seq(nrow(ICD)))
{
  ICD8 <- unlist(strsplit(ICD$ICD8[j],','))
  ICD10 <- unlist(strsplit(ICD$ICD10[j],','))
  
  cohort.data<-data.frame(cohort[which(cohort$cpr_enc %in% c(vel.neg,vel.pos)),c(1,2,3)], 'VEL'=NA,'cohort'=NA,'DIAGNOSIS'=FALSE,'DIAGNOSIS.AGE'=NA,'DIAGNOSIS.YEARS'=NA, stringsAsFactors = F)
  
  cohort.data[which(cohort.data$cpr_enc %in% vel.pos),]$VEL = '+'
  cohort.data[which(cohort.data$cpr_enc %in% vel.neg),]$VEL = '-'

  cohort.data[which(cohort.data$cpr_enc %in% c(vel.pos.dbds,vel.neg.dbds)),]$cohort = 'DBDS'
  cohort.data[which(cohort.data$cpr_enc %in% c(vel.pos.chb,vel.neg.chb)),]$cohort = 'CHB'
  
  for (i in seq(nrow(cohort.diagnosis)))
  {
    ID <- cohort.diagnosis$cpr_enc[i]

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
    glm.res<-glm(DIAGNOSIS ~ VEL + age + cohort, data = cohort.data, family=binomial)
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

  cox.res<-coxph(Surv(DIAGNOSIS.YEARS,DIAGNOSIS) ~ VEL + age + cohort, data = cohort.data, id = cohort.data$cpr_enc)
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
table3 <- data.frame('Category'=ICD$CATEGORY, 'N.diag.pos'=cox.diag.npos, 'N.diag.neg'=cox.diag.nneg, 'cox.hazard'=cox.hazard,'cox.CI95.lower'=cox.lower,'cox.CI95.upper'=cox.upper,'cox.p'=cox.p, 'glm.coef'=glm.coef,'glm.CI95.lower'=glm.lower,'glm.CI95.upper'=glm.upper,'glm.p'=glm.p, stringsAsFactors = F, check.names = F)
write.table(table3,file = paste('table3.combo.',sex,'.tsv',sep = ''), sep="\t", row.names = F, col.names = T, quote = F)

assign(paste('table3.',sex,sep=''),table3)
#View(table3)

################### TABLE4 #####################

statiner<-c('C10AA01','C10AA02','C10AA03','C10AA04','C10AA05','C10AA06','C10AA07','C10AA08')
antihypertensive <- c('C02A','C02B','C02C','C02D','C02K','C02L','C02N')
diabetes <- c('A10BA','A10BB','A10BC','A10BD','A10BF','A10BG','A10BH','A10BJ','A10BK','A10BX')

prescription <- read.table(file='prescriptions_chb_rbc_gwas.tsv', sep='\t', header = T, stringsAsFactors = FALSE)#, fill  = T)
cohort.prescription <- prescription[which(prescription$cpr_enc %in% c(vel.neg,vel.pos) ),]
cohort.data<-data.frame(cohort[which(cohort$cpr_enc %in% c(vel.neg,vel.pos)),c(1,2,3)], 'VEL'=NA,'cohort'=NA,'statiner'=FALSE,'statiner.AGE'=NA,'statiner.YEARS'=NA,'antihypertensive'=FALSE,'antihypertensive.AGE'=NA,'antihypertensive.YEARS'=NA,'diabetes'=FALSE,'diabetes.AGE'=NA,'diabetes.YEARS'=NA, stringsAsFactors = F)

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

cohort.data[which(cohort.data$cpr_enc %in% c(vel.pos.dbds,vel.neg.dbds)),]$cohort = 'DBDS'
cohort.data[which(cohort.data$cpr_enc %in% c(vel.pos.chb,vel.neg.chb)),]$cohort = 'CHB'

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
  cox.res<-coxph(Surv(statiner.YEARS,statiner) ~ VEL + age + cohort, data = cohort.data, id = cohort.data$cpr_enc)
  cox.hazard <- c(cox.hazard,summary(cox.res)$coefficients[1,2])
  cox.p <- c(cox.p,summary(cox.res)$coefficients[1,5])
  cox.lower <- c(cox.lower,summary(cox.res)$conf.int[1,3])
  cox.upper <- c(cox.upper,summary(cox.res)$conf.int[1,4])
  
  glm.res<-glm(statiner ~ VEL  + age + cohort,data = cohort.data, family = binomial)
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
  cox.res<-coxph(Surv(antihypertensive.YEARS,antihypertensive) ~ VEL + age + cohort, data = cohort.data, id = cohort.data$cpr_enc)
  cox.hazard <- c(cox.hazard,summary(cox.res)$coefficients[1,2])
  cox.p <- c(cox.p,summary(cox.res)$coefficients[1,5])
  cox.lower <- c(cox.lower,summary(cox.res)$conf.int[1,3])
  cox.upper <- c(cox.upper,summary(cox.res)$conf.int[1,4])
  
  glm.res<-glm(antihypertensive ~ VEL + age + cohort,data = cohort.data, family = binomial)
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
  cox.res<-coxph(Surv(diabetes.YEARS,diabetes) ~ VEL + age + cohort, data = cohort.data, id = cohort.data$cpr_enc)
  cox.hazard <- c(cox.hazard,summary(cox.res)$coefficients[1,2])
  cox.p <- c(cox.p,summary(cox.res)$coefficients[1,5])
  cox.lower <- c(cox.lower,summary(cox.res)$conf.int[1,3])
  cox.upper <- c(cox.upper,summary(cox.res)$conf.int[1,4])
  
  glm.res<-glm(diabetes ~ VEL + age + cohort,data = cohort.data, family = binomial)
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
write.table(table4,file = paste('table4.combo.',sex,'.tsv',sep=''), sep="\t", row.names = F, col.names = T, quote = F)

assign(paste('table4.',sex,sep=''),table4)
#View(table4)

}
#################### END ####################
View(read.table('table3.combo.F.tsv', header = T, sep = '\t')[,c(1,2,3,8,9,10,11)])
View(read.table('table2.dbds.tsv', header = T, sep = '\t'))
