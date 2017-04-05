def print_help = {
    log.info ''
    log.info '     BHIVE BED pipeline'
    log.info '-----------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    bed.nf --integs HIV_INTEGS [optional arguments]'
    log.info ''
    log.info 'Computes HIV expression from PCR and RT-PCR reads.'
    log.info 'Additional params, such as datasets for each replicate,'
    log.info 'flanking sequences, etc. should be defined in "params.txt".'
    log.info ''
    log.info 'Required arguments:'
    log.info ''
    log.info '  --integs      Path to file containing HIV integrations.'
    log.info ''
    log.info 'Optional arguments:'
    log.info '  --options     Path to alternative mapping params file (default: params.txt)'
    log.info '  --out         Output path (default: .)'
    log.info ''
}

lib_dir = "../lib"

// Set defaults
params.help        = false
params.out         = '.'
params.options     = 'params.txt'
params.integs      = null

/**********************************************************************
**************************** PARSE  PARAMS ****************************
**********************************************************************/

if (params.help) {
   print_help()
   exit 0
}

if (params.integs == null) {
   log.info "error: HIV integration file not specified!"
   exit 1
}

ifile = file("${params.integs}")
if (!ifile.isFile()) {
   log.info "error: HIV integration file not found in '${ifile}'."
   exit 1
}
integ_file = Channel.from(ifile)

/**********************************************************************
***************************** PARSE OPTIONS ***************************
**********************************************************************/
/*
** Bed datasets format:
** 
** datasets:
** feature_name, file/URL/SRR
** 
*/

def options = ['ignore_strand': null, 'chr_filter_re': null]

mfile = file("${params.options}")
if (!mfile.isFile()) {
   log.info "error: options file not found (${mfile})!"
   exit 1
}

// Varaibles
def status = 0

// Channels
file_getref = Channel.create()
datasets    = Channel.create()

// Parse file by lines.
mfile.eachLine { line ->
   if (line =~ /^#/ || line =~ /^\s*$/) {
      return null
   } else if (line =~ /^datasets:/) {
      status = 1
      return null
   }
   switch (status) {
   case 0:
     t = line.split(/\s*=\s*/).toList()
     if (!options.containsKey(t[0])) {
        log.info "unknown option: '${t[0]}' in $params.options"
        exit 1
     }
     options[t[0]] = t[1]
     break
   case 1:
     t = line.split(/\s*,\s*/).toList()
     if (t.size() == 2) {
        if (t[1] =~ /^SRR|^ERR|^DRR/ || t[1] =~ /^http|^ftp/) {
           file_getref << [ t[1], t[0] ]
        } else {
           reads = file("${t[1]}")
           if (reads.isFile()) {
              datasets << [ reads, t[0] ]
           } else {
              log.info "WARNING: ${reads} file not found, in '$line'. Ignoring BED file."
              exit 1
           }
        }
     } else {
        log.info "error in entry '$line': rna/dna entries must specify a file path, an URL or a GEO reference starting with SRR/ERR/DRR."
        exit 1
     }
     break
   }
}

// Parse options
if (!options['ignore_strand']) {
   log.info "error: 'ignore_strand' option must be defined in $params.options before 'datasets:'."
   exit 1
}

// Group the same references/URLS to avoid redundant downloads.
file_getref.close()

process getDataset {
   // Process options
   tag "${data[0]}_${data[1]}"
   publishDir path:"${params.out}/datasets/", mode:'symlink'
   // Cluster options
   cpus 1
   memory '2GB'

   input:
   val data from file_getref
   output:
   set file('*.bed'), sample into datasets
   script:
   sample = data[1]
   ref = data[0]
   if (ref =~ /^SRR|^ERR|^DRR/) {
      """
      fastq-dump --split-files -Z -A ${ref} > ${sample}.bed
      rm -f ~/ncbi/public/sra/${ref}.sra
      """
   } else {
      """
      wget ${ref} -O | zcat -f > ${sample}.bed
      """
   }
}

process formatDatasets {
   // Process options
   tag "${sample}.bed"
   // Cluster options
   cpus 1
   memory '2GB'

   input:
   set file(bedfile), sample from datasets
   output:
   file "${sample}.bed" into datasets_
   script:
   """
   mv ${bedfile} ${sample}.bed
   """
}

/**********************************************************************
*************************  BEDFILE PIPELINE  **************************
**********************************************************************/

process ComputeGenomicDistance {
   // Process options
   tag "${bedset.size()} datasets"
   publishDir path:"${params.out}/", mode:'move'
   // Cluster options
   cpus 1
   memory '32GB'

   input:
   file bedset from datasets_.flatten().toList()
   file integs from integ_file

   output:
   file 'hiv_features.txt' into out_files

   script:
   """
   #!/usr/bin/env Rscript
   library(GenomicRanges)
   
   # Load integ file.
   integs = read.table(file='${integs}', header=TRUE, comment.char='#')

   # Integ GRange
   gins = GRanges(Rle(integs\$chr), IRanges(start=integs\$locus, width=1),
          strand=Rle(integs\$strand), brcd=integs\$brcd,
          nread=integs\$reads, mapq=integs\$mapq, rep=integs\$rep)

   # Load bed sets
   beds = Sys.glob('*.bed')
   gbed = list(NA)
   bednames <- c()
   for (file in beds) {
     sname  = strsplit(as.character(file),'.',fixed=TRUE)[[1]][1]
     write(paste(paste('detect file/name:',file),sname), stdout())
     bedset = ${options['chr_filter_re'] ?
         "subset(read.table(file, comm='#'),!grepl('${options['chr_filter_re']}', V1))" :
         "read.table(file, comm='#')"
     }
     if (ncol(bedset) == 3 || (ncol(bedset) > 3 && (bedset[1,4] != '+' && bedset[1,4] != '-'))) {
       gbed[[sname]] = reduce(GRanges(Rle(bedset\$V1), IRanges(start=bedset\$V2, end=bedset\$V3)))
     } else if (ncol(bedset) > 3 && (bedset[1,4] == '+' || bedset[1,4] == '-')) {
       gbed[[sname]] = reduce(GRanges(Rle(bedset\$V1), IRanges(start=bedset\$V2, end=bedset\$V3), strand=Rle(bedset\$V4)))
     } else {
       write(paste(paste('error: unknown format for bed dataset:',file),'. Correct format or remove from datasets.\n'), stderr())
       quit(save='no', status=1) 
     }
     bednames <- c(bednames,sname)
     write(paste(paste('read dataset file/name:',file),sname), stdout())
   }

   # Compute distanceToNearest
   for (sname in bednames) {
     dist = distanceToNearest(gins, gbed[[sname]], ignore.strand=${options['ignore_strand']})
     integs[sname] = rep(NA,length(gins))
     integs[queryHits(dist),sname] = dist@elementMetadata\$distance
   }

   # Write table
   write.table(integs, file='hiv_features.txt', sep='\\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
   """
}
