library(dplyr)
library(venn)

setwd("~/Documents/GiantClamBUSCO")
# remember deleting first 3 lines (#) in full_table.tsv
# Metazoa
Tmaxima <- read.table(file = 'Tmaxima_metazoa.tsv', sep = '\t', header = FALSE, fill=TRUE)
Tmaxima_M <- subset(filter(Tmaxima, V2 == 'Missing'), select = c(1))
write.table(Tmaxima_M, file = "Tmaxima_metazoa_missing.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

Tcrocea <- read.table(file = 'Tcrocea_metazoa.tsv', sep = '\t', header = FALSE, fill=TRUE)
Tcrocea_M <- subset(filter(Tcrocea, V2 == 'Missing'), select = c(1))
write.table(Tcrocea_M, file = "Tcrocea_metazoa_missing.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

Tgigas <- read.table(file = 'Tgigas_metazoa.tsv', sep = '\t', header = FALSE, fill=TRUE)
Tgigas_M <- subset(filter(Tgigas, V2 == 'Missing'), select = c(1))
write.table(Tgigas_M, file = "Tgigas_metazoa_missing.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")



maxima <- readLines("Tmaxima_metazoa_missing.tsv")
crocea <- readLines("Tcrocea_metazoa_missing.tsv")
gigas <- readLines("Tgigas_metazoa_missing.tsv")

# get the intersection
inter_metazoa<-intersect(maxima,intersect(crocea,gigas))

write(inter_metazoa, "inter_metazoa_busco_ids.txt")

missing_list <- list (
  Tmaxima=maxima, 
  Tcrocea=crocea, 
  Tgigas=gigas 
  )

venn(missing_list, ilabels = TRUE, 
     zcolor = "style", size = 10, cexil = 12, cexsn = 1.5, ilcs = 3, sncs = 2.5)

### Mollusca dataset
## clear the env first before running
## adding quote=""
Tmaxima <- read.table(file = 'Tmaxima_mollusca.tsv', sep = '\t', header = FALSE, fill=TRUE, quote="")
Tmaxima_M <- subset(filter(Tmaxima, V2 == 'Missing'), select = c(1))
write.table(Tmaxima_M, file = "Tmaxima_mollusca_missing.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
## complete list 
write.table(Tmaxima$V1, file = "Tmaxima_mollusca_full.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


Tcrocea <- read.table(file = 'Tcrocea_mollusca.tsv', sep = '\t', header = FALSE, fill=TRUE, quote="")
Tcrocea_M <- subset(filter(Tcrocea, V2 == 'Missing'), select = c(1))
write.table(Tcrocea_M, file = "Tcrocea_mollusca_missing.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

Tgigas <- read.table(file = 'Tgigas_mollusca.tsv', sep = '\t', header = FALSE, fill=TRUE, quote="")
Tgigas_M <- subset(filter(Tgigas, V2 == 'Missing'), select = c(1))
write.table(Tgigas_M, file = "Tgigas_mollusca_missing.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")



maxima <- readLines("Tmaxima_mollusca_missing.tsv")
crocea <- readLines("Tcrocea_mollusca_missing.tsv")
gigas <- readLines("Tgigas_mollusca_missing.tsv")

# get the intersection
inter_mollusca<-intersect(maxima,intersect(crocea,gigas))
# save to file
write(inter_mollusca, "inter_mollusca_busco_ids.txt")


missing_list <- list (
  Tmaxima=maxima, 
  Tcrocea=crocea, 
  Tgigas=gigas 
)

venn(missing_list, ilabels = TRUE, 
     zcolor = "style", size = 10, cexil = 12, cexsn = 1.5, ilcs = 3, sncs = 2.5)



