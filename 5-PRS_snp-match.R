
# GWAS summary is from EU-ALS cohort consisting of 20,804 ALS cases;http://als.umassmed.edu/
# reference:Genome-wide Analyses Identify KIF5A as a Novel ALS Gene;https://doi.org/10.1016/j.neuron.2018.02.027

library(data.table)
info=fread("alsMetaSummaryStats_march21st2018.tab",data.table = F ) #EU-ALS
head(info)

dim(info)
read.table("Als_Methy_impute.bim", header = F) -> bim #PUTH-ALS
head(bim)
bim$SNP <- paste0(bim$V1, ":", bim$V4)
length(which(bim$SNP %in% info$SNP))
dim(bim)
length(which(info$Allele1 == "a" & info$Allele2 == "t"))
length(which(info$Allele1 == "c" & info$Allele2 == "g"))
length(which(info$Allele1 == "t" & info$Allele2 == "a"))
length(which(info$Allele1 == "g" & info$Allele2 == "c"))
del_pos <- c(which(info$Allele1 == "a" & info$Allele2 == "t"),
             which(info$Allele1 == "c" & info$Allele2 == "g"),
             which(info$Allele1 == "t" & info$Allele2 == "a"),
             which(info$Allele1 == "g" & info$Allele2 == "c"))
del_pos
info_del <- info[-del_pos, ]
bim_update <- bim[which(bim$SNP %in% info_del$SNP), ]
which(duplicated(bim_update$SNP)==TRUE) -> tpm_pos
snp_dup <- bim_update$SNP[tpm_pos]
bim_update <- bim_update[-which(bim_update$SNP %in% snp_dup),]
info_del_update <- info_del[which(info_del$SNP %in% bim_update$SNP), ]
info_del_update <- info_del_update[order(info_del_update$CHR, info_del_update$BP), ]
tpm_pos <- c(which(tolower(bim_update$V5) == info_del_update$Allele1),
             which(tolower(bim_update$V5) == info_del_update$Allele2))
bim_update <- bim_update[tpm_pos, ]
info_del_update <- info_del_update[tpm_pos, ]
write.table(bim_update$V2, "snp.keep.txt", col.names = F, row.names = F, sep = "\t", quote = F)
info_del_update$SNP <- bim_update$V2
info_del_update$Allele1 <- toupper(info_del_update$Allele1)
write.table(info_del_update[, c(1, 2, 4, 6)], "score.txt", col.names = F, row.names = F, sep = "\t", quote = F)

dim(info_del_update)
length(which(info_del_update$P<"0.05"))#num:253824
length(which(info_del_update$P<"0.001"))
length(which(info_del_update$P<"0.1"))
length(which(info_del_update$P<"0.0001"))

