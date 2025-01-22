library(readxl)
dt<-as.data.table(read_excel('data/High_effect_trans_qtls.xlsx'))
conversiontable<-fread("data/original/transcriptome_conversion.csv")
conversiontable[,AGI := toupper(AGI)]

setDT(dt)
setkey(dt,txname)
setDT(conversiontable)
setkey(conversiontable,AGI)
finaltable<-merge(dt,conversiontable, by.x = 'txname', by.y = 'AGI',all.x = T, all.y = F)
fwrite(finaltable,"High_effect_size_trans_eQTL_Jan22.txt")


