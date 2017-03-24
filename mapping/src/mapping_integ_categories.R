library(GenomicRanges)

gen_path <- Sys.getenv("RNASQ_PATH", unset="")
hm_dir <- Sys.getenv("HM_DIR", unset="")
int_file <- Sys.getenv("INT_FILE", unset="")
out_file <- Sys.getenv("OUT_FILE", unset="")
chrsz_path <- Sys.getenv("CHRSZ_PATH",unset="")
unmap_path <- Sys.getenv("UNMAP_PATH",unset="")

if (gen_path == "" || hm_dir == "" || int_file == "" || out_file == "" || chrsz_path == "" || unmap_path == ""){
   print("Parameters: <integ_file> <output_file>. Additionally, must set env variables 'INT_FILE', 'OUT_FILE', 'RNASQ_PATH', 'CHRSZ_PATH', 'UNMAP_PATH' and 'HM_DIR'.")
   quit(save="no", status=1)
}

hm_sets = subset(data.frame(fname=dir(hm_dir)),grepl('.01$',fname))
hm_cnt = nrow(hm_sets)
hm_name = unlist(strsplit(as.character(hm_sets$fname),"\\."))[seq(1,2*hm_cnt,2)]

cat("\n#### Input files ####\n")
cat(paste("working directory:", getwd(), "\n"))
cat(paste("integrations:", int_file, "\n"))
cat(paste("genes:", gen_path, "\n"))

cat(paste("chr sizes:", chrsz_path, "\n\n"))
cat(paste("Histone marks found:", hm_cnt, "\n"))
cat(paste("Histone marks:", hm_name, "\n"))

# Read genome
genome = read.table(chrsz_path, as.is=TRUE)

# Load integrations
integs = subset(read.table(int_file, as.is=TRUE),!grepl("HIV|random", V2))
colnames(integs) = c("brcd", "chrom", "pos", "strand", "nread", "mapq", "rep")

# Shortcuts (suppresses the warning that function has changed).
quietDistToN = function(...) suppressWarnings(distanceToNearest(...))


# Read gene annotation and expression file.
genes = subset(read.table(gen_path, header=T), type == "protein_coding")
cutoff = quantile(genes$tpm, .4)
genact = subset(genes, tpm > cutoff)
gensil = subset(genes, tpm <= cutoff)

cat("\n#### Sanity check for genes ####\n")
cat(paste("number:", nrow(genes), "\n"))
cat(paste("active:", nrow(genact), "\n"))
cat(paste("silent:", nrow(gensil), "\n"))
cat(paste("cutoff tpm:", cutoff, "\n"))

# Read Histone Marks.
cat("\n#### Sanity check for histone marks ####\n")
hm_data = list(NA)
for (i in 1:hm_cnt) {
    hm_data[[hm_name[i]]] = read.table(paste(hm_dir,hm_sets$fname[i],sep="/"), comm="#")
    cat(paste("histone mark:", hm_name[i], "\n"))
    cat(paste("total rows:",nrow(hm_data[[hm_name[i]]]),"\n"))
    hm_data[[hm_name[i]]] = subset(hm_data[[hm_name[i]]],!grepl("random|chrUn", V1))
    cat(paste("rows without random chr:",nrow(hm_data[[hm_name[i]]]),"\n"))
    hm_dataframe = hm_data[[hm_name[i]]]
    cat(paste("number:", nrow(hm_dataframe), "\n"))
    cat(paste("coverage:", sum(hm_dataframe$V3-hm_dataframe$V2+1), "\n"))
}

######################################################################
# Create GRanges objects.

ggnom = GRanges(Rle(genome$V1), IRanges(start=1, width=genome$V2))
gsize = sum(as.numeric(sum(coverage(ggnom))))
cat("\n#### Sanity check for genome (GRanges) ####\n")
cat(paste("chromosomes:", length(ggnom), "\n"))
cat(paste("coverage:", gsize, "\n"))

# Load unmappable genome
unmap = subset(read.table(unmap_path),!grepl("HIV|random|chrUn", V1))
gumap = GRanges(Rle(unmap$V1), IRanges(start=unmap$V2, end=unmap$V3))
usize = sum(as.numeric(sum(coverage(gumap))))
cat("\n#### Sanity check for unmappable genome (GRanges) ####\n")
cat(paste("coverage:", usize, "\n"))
cat(paste("mappable genome:",(gsize-usize)/gsize*100.0,"%\n"))


# Insertions.
gins = GRanges(Rle(integs$chrom), IRanges(start=integs$pos, width=1),
	strand=Rle(integs$strand), brcd=integs$brcd,
	nread=integs$nread, mapq=integs$mapq, rep=integs$rep)

cat("\n#### Sanity check for insertions (GRanges) ####\n")
cat(paste("number:", length(gins), "\n"))
cat(paste("coverage:", sum(sum(coverage(gins))), "\n"))


# Active/silent genes
ggen = GRanges(Rle(genes$chr), IRanges(start=genes$beg, end=genes$end),
   strand=Rle(genes$strand), name=genes$gene_name, expr=genes$tpm)
actgen = GRanges(Rle(genact$chr), IRanges(start=genact$beg, end=genact$end),
   strand=Rle(genact$strand), name=genact$gene_name, expr=genact$tpm)
silgen = GRanges(Rle(gensil$chr), IRanges(start=gensil$beg, end=gensil$end),
   strand=Rle(gensil$strand), name=gensil$gene_name, expr=gensil$tpm)
acttss = resize(actgen,1)
siltss = resize(silgen,1)

cat("\n#### Sanity check for genes (GRanges) ####\n")
cat(paste("coverage genes:", sum(sum(coverage(ggen))), "\n"))
cat(paste("coverage active:", sum(sum(coverage(actgen))), "\n"))
cat(paste("coverage silent:", sum(sum(coverage(silgen))), "\n"))

# Promoters
actprom = promoters(actgen, up=2500, down=2500)
#silprom = promoters(silgen, up=2500, down=2500)
cat("\n#### Sanity check for promoters (GRanges) ####\n")
cat(paste("coverage active:", sum(sum(coverage(actprom))), "\n"))
#cat(paste("coverage silent:", sum(sum(coverage(silprom))), "\n"))

# Histone Mark Ranges
cat("\n#### Sanity check for Histone marks ####\n")
ghm = list(NA)
for (i in 1:hm_cnt) {
     hm_df = hm_data[[hm_name[i]]]
     ghm[[hm_name[i]]] = GRanges(Rle(hm_df$V1), IRanges(start=hm_df$V2, hm_df$V3))
     ghm[[hm_name[i]]] = reduce(ghm[[hm_name[i]]])
}

genh_ = ghm[['H3K27ac']]
start(genh_) <- start(genh_) - 2500
end(genh_) <- end(genh_) + 2500
cat("\n#### Sanity check for enhancers (GRanges) ####\n")
cat(paste("coverage genh_:", sum(sum(coverage(genh_))), "\n"))

######################################################################

# Record distances to closest gene, promoter and enhancer.
d_ag = quietDistToN(gins, actgen, ignore.strand = TRUE)
d_at = quietDistToN(gins, acttss, ignore.strand = TRUE)
d_sg = quietDistToN(gins, silgen, ignore.strand = TRUE)
d_st = quietDistToN(gins, siltss, ignore.strand = TRUE)
d_ap = quietDistToN(gins, actprom, ignore.strand = TRUE)
#d_sp = quietDistToN(gins, silprom)

# Assign gene names to matches
gins$gene = 'NA'
ovgen = as.data.frame(findOverlaps(gins,silgen, ignore.strand = TRUE))
gins$gene[ovgen$queryHits] = as.character(silgen$name[ovgen$subjectHits])
ovgen = as.data.frame(findOverlaps(gins,actgen, ignore.strand = TRUE))
gins$gene[ovgen$queryHits] = as.character(actgen$name[ovgen$subjectHits])

gins$d_actgen = rep(NA,length(gins))
gins$d_actgen[queryHits(d_ag)] = d_ag@elementMetadata$distance
gins$d_acttss = rep(NA,length(gins))
gins$d_acttss[queryHits(d_at)] = d_at@elementMetadata$distance
gins$d_silgen = rep(NA,length(gins))
gins$d_silgen[queryHits(d_sg)] = d_sg@elementMetadata$distance
gins$d_siltss = rep(NA,length(gins))
gins$d_siltss[queryHits(d_st)] = d_st@elementMetadata$distance
gins$d_actprom = rep(NA,length(gins))
gins$d_actprom[queryHits(d_ap)] = d_ap@elementMetadata$distance
#gins$d_silprom = d_sp@elementMetadata$distance

# Subtract overlap between active and inactive genes.
silgen_na = setdiff(silgen, actgen, ignore.strand = TRUE)

# Subtract promoter regions from genes.
actgen_np = setdiff(actgen, actprom, ignore.strand = TRUE)
silgen_nanp = setdiff(silgen_na, actprom, ignore.strand = TRUE)

# Subtract enhancers from genes and promoters
actprom_ne = setdiff(actprom, genh_, ignore.strand = TRUE)
#silprom = setdiff(silprom, genh_, ignore.strand = TRUE)
actgen_npne = setdiff(actgen_np, genh_, ignore.strand = TRUE)
silgen_nanpne = setdiff(silgen_nanp, genh_, ignore.strand = TRUE)

# Subtract unmappable regions from genes and promoters
actprom_neu = setdiff(actprom_ne, gumap, ignore.strand = TRUE)
actgen_npneu = setdiff(actgen_npne, gumap, ignore.strand = TRUE)
silgen_nanpneu = setdiff(silgen_nanpne, gumap, ignore.strand = TRUE)

# Generate intergenic GRanges
ginter = setdiff(ggnom, gumap, ignore.strand = TRUE)
ginter = setdiff(ginter, actprom_neu, ignore.strand = TRUE)
ginter = setdiff(ginter, actgen_npneu, ignore.strand = TRUE)
ginter = setdiff(ginter, silgen_nanpneu, ignore.strand = TRUE)
ginter = setdiff(ginter, genh_, ignore.strand = TRUE)

#grest = setdiff(ggnom, actgen_npne, ignore.strand = TRUE)
#grest = setdiff(grest, silgen_ne, ignore.strand = TRUE)
#grest = setdiff(grest, actprom_ne, ignore.strand = TRUE)
#grest = setdiff(grest, silprom, ignore.strand = TRUE)
#grest = setdiff(grest, genh_, ignore.strand = TRUE)

# Coverages
ag_len = sum(as.numeric(sum(coverage(reduce(actgen_npneu,ignore.strand=T)))))
ap_len = sum(as.numeric(sum(coverage(reduce(actprom_neu,ignore.strand=T)))))
sg_len = sum(as.numeric(sum(coverage(reduce(silgen_nanpneu,ignore.strand=T)))))
#sp_len = sum(sum(coverage(silprom)))
#en_len = 2*sum(sum(coverage(genh_)))
en_len = sum(as.numeric(sum(coverage(genh_))))
#in_len = sum(sum(coverage(grest)))
#in_len = 2*sum(as.numeric(genome$V2)) - ag_len - ap_len - sg_len - en_len
#in_len = sum(as.numeric(genome$V2)) - ag_len - ap_len - sg_len - en_len
in_len = sum(as.numeric(sum(coverage(ginter))))

cat("\n#### Sanity check (GRanges) ####\n")
cat(paste("coverage active genes:", ag_len, "\n"))
cat(paste("coverage silent genes:", sg_len, "\n"))
cat(paste("coverage active proms:", ap_len, "\n"))
cat(paste("coverage enhancers:", en_len, "\n"))
cat(paste("coverage intergenic:", in_len, "\n"))
cat(paste("coverage unmappable:", usize, "\n"))
cat(paste("total:", as.numeric(ag_len+sg_len+ap_len+en_len+in_len+usize)),"\n")

#ovag = countOverlaps(gins, actgen, ignore.strand=TRUE) > 0
#ovap = countOverlaps(gins, actprom, ignore.strand=TRUE) > 0
#ovsg = countOverlaps(gins, silgen, ignore.strand=TRUE) > 0
#ovsp = countOverlaps(gins, silprom, ignore.strand=TRUE) > 0
#oven = countOverlaps(gins, genh_, ignore.strand=TRUE) > 0

ovag = countOverlaps(gins, actgen, ignore.strand=TRUE) > 0
ovap = countOverlaps(gins, actprom, ignore.strand=TRUE) > 0
ovsg = countOverlaps(gins, silgen, ignore.strand=TRUE) > 0
#ovsp = countOverlaps(gins, silprom, ignore.strand=TRUE) > 0
oven = countOverlaps(gins, genh_, ignore.strand=TRUE) > 0


# Assign categories (important to write in this order to preserve EN > PROM > GENE)
gins$cat = 'IN'
gins[ovsg,]$cat = 'SG'
gins[ovag,]$cat = 'AG'
gins[ovap,]$cat = 'AP'
#gins[ovsp,]$cat = 'SP'
gins[oven,]$cat = 'EN'

# Write integ table
out_table = data.frame(
	brcd      = gins$brcd,
	chrom     = seqnames(gins),
	pos       = start(gins),
	strand    = strand(gins),
	nread     = gins$nread,
	mapq      = gins$mapq,
	cat       = gins$cat,
	gene_name = gins$gene,
	d_actgen  = gins$d_actgen,
	d_acttss  = gins$d_acttss,
	d_silgen  = gins$d_silgen,
	d_siltss  = gins$d_siltss,
	d_actprom = gins$d_actprom
#	d_silprom = gins$d_silprom
)

# Add distance to histone marks.
for (i in 1:hm_cnt) {
     cat(paste("computing distances to:",hm_name[i],"(index",i,")\n"))
     d_hm = distanceToNearest(gins, ghm[[hm_name[i]]])
     cat(paste("done, saving",hm_name[i],"\n"))
     out_table[hm_name[i]] = rep(NA,length(gins))
     out_table[queryHits(d_hm),hm_name[i]] = d_hm@elementMetadata$distance
     cat(paste("OK\n"))
}

out_table$rep = gins$rep

head(out_table)

# Write a header.
write(file=out_file, paste('#AG\t', ag_len, sep=""))
write(file=out_file, paste('#SG\t', sg_len, sep=""), append=TRUE)
write(file=out_file, paste('#AP\t', ap_len, sep=""), append=TRUE)
#write(file=out_file, paste('#SP\t', sp_len, sep=""), append=TRUE)
write(file=out_file, paste('#EN\t', en_len, sep=""), append=TRUE)
write(file=out_file, paste('#IN\t', in_len, sep=""), append=TRUE)
write(file=out_file, '# brcd     \tHIV barcode\n# chrom    \tintegration chromosome\n# pos      \tintegration locus\n# strand   \t+ forward/- reverse\n# nread    \tiPCR reads\n# mapq     \tbarcode-locus assignment score\n# cat      \tgenome region (AG: active gene, SG: silent gene, AP: active promoter, EN: enhancer, IN: intergenic)\n# gene_name\tgene name (NA: not inside a gene)\n# d_actgen \tdistance from closest active gene body\n# d_acttss \tdistance from closest active gene TSS\n# d_silgen \tdistance from closest silent gene body\n# d_siltss \tdistance from closest silent gene TSS\n# d_actprom\tdistance from promoter of closest active gene', append=TRUE)
for (i in 1:hm_cnt) {
     write(file=out_file, paste('# ',hm_name[i],'\tdistance from closest ',hm_name[i],' mark', sep=''), append=TRUE)
}
write(file=out_file, '# rep      \tbiological replicate', append=TRUE)

write.table(out_table, file = out_file, sep="\t", quote=FALSE, row.names=FALSE, append=TRUE)
