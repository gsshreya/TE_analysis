library(data.table)

gtex = read.table("eqtls_gtex_hg38.txt", header=T,comment.char="", sep="\t", quote="", nrow = 30147)

 gtex$pos = gsub("_INS_.*","", gtex$snps_hg38)
 gtex$pos = gsub("_DEL_.*","", gtex$pos)
 gtex$pos = paste0("chr", gsub("_",":", gtex$pos))



#############################
ourvars = fread("zcat /gpfs/data/user/shweta_lab/data/TE/TE_discovery/Jan_May_2026/Filtered_VCFs/ALL_TE.merged.new.vcf.gz | grep -v '#' |cut -f 1-5", data.table=F) 
ourvars = ourvars[which(grepl("INS", ourvars$V5 )), ]

ourvars$pos = paste0(ourvars[,1], ":", ourvars[,2])

indigen = fread("zcat Indigen_Alu_final_geno10_all_22K.vcf.gz | grep -v '#' | cut -f 1-5", data.table=FALSE)
indigen$pos = paste0("chr", indigen$V1, ":", indigen$V2)

popfreq = read.table("GI_TE_afs/per_population_AF_wide.tsv", header=T, stringsAsFactors=F, sep="\t")
popfreq$pos2 = paste0( popfreq$chrom, ":", popfreq$pos,popfreq$end)
popfreq$pos3 = paste0(popfreq[,1], ":", popfreq[,2])

popfreq2 <- popfreq[, c("pos2", "pos3", 'alt', grep("^AF_", names(popfreq), value = TRUE))]

#rownames(popfreq) = popfreq$pos2

spirito = read.table("spirito_2018_eqtls.txt", header=T,sep="\t", stringsAsFactors=F)
spirito$pos = paste0("chr", spirito$Chr, ":", spirito$pos)
setDT(spirito)

library(data.table)

setDT(gtex)
setDT(ourvars)

# Keep only insertions
gtex_ins <- gtex[grepl("INS", snps_hg38)]
#our_ins  <- ourvars[grepl("INS", V5)]
# Extract chromosome and position

#gtex_ins = gtex
our_ins = ourvars

# Extract chromosome and position

gtex_ins[, chr := sub(":.*", "", as.character(pos))]
gtex_ins[, gtex_pos := as.integer(sub(".*:", "", as.character(pos)))]

our_ins[, `:=`(
  chr = V1,
  our_pos = as.integer(V2)
)]

# Create search interval
our_ins[, `:=`(
  start = our_pos - 20L,
  end   = our_pos + 20L
)]

# Range overlap join
setkey(gtex_ins, chr, gtex_pos)
setkey(our_ins, chr, start, end)

gtex_query <- gtex_ins[!is.na(start) & !is.na(end) & !is.na(chr)]
gtex_ins <- gtex_query

matches <- foverlaps(
  gtex_ins[, .(chr, start = gtex_pos, end = gtex_pos, gtex_pos, snps_hg38, gene,
               tissue, pvalue, beta)],
  our_ins,
  by.x = c("chr", "start", "end"),
  type = "within",
  nomatch = 0L
)

matches


#spirito 
setDT(spirito)
spirito_ins = spirito[which(grepl("^ALU_|^L1_", spirito$TE_ID)),]

#spirito_ins = spirito

# Extract chromosome and position

spirito_ins[, chr := sub(":.*", "", as.character(pos))]
spirito_ins[, spirito_pos := as.integer(sub(".*:", "", as.character(pos)))]

setkey(spirito_ins, chr, spirito_pos)


matches_sp <- foverlaps(
  spirito_ins[, .(chr, start = spirito_pos, end = spirito_pos, spirito_pos, TE_ID, gene,
                P.value, slope)],
  our_ins,
  by.x = c("chr", "start", "end"),
  type = "within",
  nomatch = 0L
)

matches_sp


#indigen

setDT(indigen)
indigen_ins = indigen

# Extract chromosome and position

indigen_ins[, chr := sub(":.*", "", as.character(pos))]
indigen_ins[, indigen_pos := as.integer(sub(".*:", "", as.character(pos)))]

indigen_query <- indigen_ins[, .(
  chr,
  start = indigen_pos,
  end   = indigen_pos,
  indigen_pos
)]

indigen_query <- indigen_query[
  !is.na(chr) & !is.na(start) & !is.na(end)
]


gtex2 = gtex_ins[-which(gtex_ins$snps_hg38 %in% matches$snps_hg38),]
gtex3 = gtex2[-which(grepl("LINE1",gtex2$snps_hg38)),]

gtex_query <- gtex3[, .(
  chr,
  start = gtex_pos - 20L,
  end   = gtex_pos + 20L,
  gtex_pos
)]

gtex_query <- gtex_query[
  !is.na(chr) & !is.na(start) & !is.na(end)
]

setkey(gtex_query, chr, start, end)

matches_indigen <- foverlaps(
  indigen_query,
  gtex_query,
  type = "within",
  nomatch = 0L
)

matches_indigen
