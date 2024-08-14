

setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/Eric Johnson/VIDUS/candidate genes study/")
pheno = read.table("Data/master.pheno.drugs.5cellTypes.txt", header=T)




#### numbers ####
apply(pheno[,c("heroin_l6m",
               "opioid_l6m",
               "marij_noninj_l6m",
               "cocaine_l6m",
               "crack_l6m",
               "anycoc_l6m",
               "meth_l6m",
               "stim_l6m",
               "anydrug_l6m")], 2, function(x) sum(x))


#### logistic regression ####
#Age yr. (mean±SD)
mod = glm(hiv ~ ageatint, data=pheno, family="binomial")
pval = summary(mod)$coefficients[2,4]
p.table = data.frame("Covariate" = "age",
                     "P-raw" = pval)
  

#Sex no. (male, %)
mod = glm(hiv ~ female, data=pheno, family="binomial")
pval = summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("sex",pval))

#Cocaine no. (%)
mod = glm(hiv ~ cocaine_l6m, data=pheno, family="binomial")
pval=summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("cocaine",pval))


#Crack no. (%)
mod = glm(hiv ~ crack_l6m, data=pheno, family="binomial")
pval=summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("crack",pval))


#Any cocaine no. (%)
mod = glm(hiv ~ anycoc_l6m, data=pheno, family="binomial")
pval=summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("any_cocaine",pval))

#Heroin no. (%)
mod = glm(hiv ~ heroin_l6m, data=pheno, family="binomial")
pval=summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("heroin",pval))

#Meth. no. (%)
mod = glm(hiv ~ meth_l6m, data=pheno, family="binomial")
pval=summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("meth",pval))

#Opioid no. (%)
mod = glm(hiv ~ opioid_l6m, data=pheno, family="binomial")
pval=summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("opioid",pval))

#Marijuana no. (%)
mod = glm(hiv ~ marij_noninj_l6m, data=pheno, family="binomial")
pval=summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("marijuana",pval))

#Any stimulant no. (%)
mod = glm(hiv ~ stim_l6m, data=pheno, family="binomial")
pval=summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("stimulants",pval))

#Any non-prescript no. (%)
mod = glm(hiv ~ anydrug_l6m, data=pheno, family="binomial")
pval=summary(mod)$coefficients[2,4]
p.table = rbind(p.table, c("any_non_prescription",pval))



#### adjust p-values ####
p.table$P.bon = p.adjust(p.table$P.raw, method="bonferroni")
p.table

