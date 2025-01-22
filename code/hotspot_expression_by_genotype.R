geno<-fread("data/geno_editednames.txt")
hotspots_geno<-geno[,c("ID","At2g24150","At5g46760")]
hotspots_geno[, genotype := NA_character_]
hotspots_geno[At2g24150 == "C" & At5g46760 == "C", genotype := "BB"]
hotspots_geno[At2g24150 == "L" & At5g46760 == "L", genotype := "SS"]
hotspots_geno[At2g24150 == "C" & At5g46760 == "L", genotype := "BS"]
hotspots_geno[At2g24150 == "L" & At5g46760 == "C", genotype := "SB"]
final_result <- hotspots_geno[, .(ID, At2g24150, At5g46760, genotype)]
final_result <- final_result[!is.na(genotype)]
hs_genos<-final_result[,.(ID,genotype)]

tx<-fread("data/rils_pheno.txt")
group2<-c("At1g75040", "At2g42220", "At2g43750", "At3g10940", "At3g57260", "At3g61080", "At3g62030", "At5g10760", "At5g24530", "At5g45380")
group3<-c("At1g74710", "At2g04450", "At3g52430", "At3g60420",  "At4g12720", "At4g39030", "At5g61900")

tx<-tx[ID %in% hs_genos$ID & treatment != "delta"]
tx_grp2<-tx[,.SD, .SDcols= c("ID","treatment",group2)]
tx_grp3<-tx[,.SD, .SDcols= c("ID","treatment",group3)]

grp2_final <- merge(hs_genos,tx_grp2, by = "ID", all.x = TRUE)
grp3_final <- merge(hs_genos,tx_grp3, by = "ID", all.x = TRUE)

grp2_avg<-grp2_final[,lapply(.SD,mean), by=.(genotype,treatment),.SDcols = patterns("^At")]
grp2_avg$genotype <- factor(grp2_avg$genotype, levels = c("BB", "SS", "BS", "SB"))

grp3_avg<-grp3_final[,lapply(.SD,mean), by=.(genotype,treatment),.SDcols = patterns("^At")]
grp3_avg$genotype <- factor(grp3_avg$genotype, levels = c("BB", "SS", "BS", "SB"))

at_columns <- grep("^At", names(grp2_avg), value = TRUE)

setDT(grp2_avg)

# Identify columns starting with "At"
at_columns <- grep("^At", names(grp2_avg), value = TRUE)


plot_expression_data_table <- function(data) {
  setDT(data)
  # Identify columns starting with "At" (expression columns)
  at_columns <- grep("^At", names(data), value = TRUE)
  long_data <- melt(data, id.vars = c("genotype", "treatment"), measure.vars = at_columns, 
                    variable.name = "Attribute", value.name = "Expression")
  p <- ggplot(long_data, aes(x = genotype, y = Expression, color = treatment)) +
    geom_line(aes(group = treatment)) +
    scale_color_manual(values = c("SA" = "#377EB8", "SW" = "#4DAF4A")) +
    labs(x = "Genotype", y = "Expression") +
    facet_wrap(~ Attribute, ncol = 2, nrow = 5, scales = "free_y") + 
    theme_minimal() +
    theme(strip.text = element_text(size = 6),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          text = element_text(size = 6))
  print(p)
}
plot_expression_data_table(grp2_avg)
plot_expression_data_table(grp3_avg)  
  

