######Codon Usage Bias Analysis coRdon package

library(coRdon)
SeqDEGs <- readSet(file="Seq.fasta")

# Codon Usage Bias Analysis by cubar
library(cubar)


s <- head(SeqDEGs, 10)
print(s)
check_cds(SeqDEGs)
cf_all <- count_codons(SeqDEGs)
dim(cf_all)

RSCU=est_rscu(cf_all, weight = 1, pseudo_cnt = 1, codon_table = get_codon_table())
write.csv(RSCU, file = "Relative Synonymous Codon Usage.csv")


CAI <- get_cai(cf_all, RSCU)
head(CAI)
hist(CAI)
write.csv(CAI, file = "Codon Adaptation Index.csv")

ENC=get_enc(cf_all, codon_table = get_codon_table())
head(ENC)
hist(ENC)
range(ENC)
write.csv(ENC, file = "Effective Number of Codons.csv")


GC <- get_gc(cf_all)
head(GC)
hist(GC)
write.csv(GC, file = "GC content.csv")

GC3S=get_gc3s(cf_all, codon_table = get_codon_table())
write.csv(GC3S, file = "GC content at synonymous third positions.csv")
