def print_help = {
   log.info ''
   log.info '     BHIVE rna-seq pipeline'
   log.info '--------------------------------'
   log.info ''
   log.info 'Usage: '
   log.info '    rnaseq.nf --index KALLISTO_INDEX [optional arguments]'
   log.info ''
   log.info 'Computes HIV integration sites from iPCR NGS reads.'
   log.info 'Additional params, such as datasets for each replicate, matching distances,'
   log.info 'mapping qualities, etc. should be defined in "params.txt".'
   log.info ''
   log.info 'Required arguments:'
   log.info ''
   log.info '  --index       Path to kallisto index file (create if it does not exist)'
   log.info ''
   log.info 'Optional arguments:'
   log.info '  --options     Path to alternative mapping params file (default: params.txt)'
   log.info '  --out         Output path (default: .)'
   log.info ''
}

src_dir = "src"
lib_dir = "../lib"

// Set defaults
params.help        = false
params.out         = '.'
params.options     = 'params.txt'
params.index       = null


if (params.help) {
   print_help()
   exit 0
}

/**********************************************************************
 ***************************** PARSE OPTIONS ***************************
 **********************************************************************/
/*
** Mapping options format:
** 
**  BFS (Barcode flanking sequence)
**  LTR (HIV LTR sequence, must be present in reverse read)
**  RFS (Restriction enzyme flanking sequence, in HIV construct)
**  DIST (Barcode clustering distance)
**  MAPQ (Minimum mapping quality, 20)
**  INTQ (Minimum assignment score, 10)
** 
** datasets:
** biological replicate,filename/URL/{SRR,ERR,DRR} reference
** 
*/

def options = ['rnalen_mean':null, 'rnalen_sd':null, 'transcripts':null, 'gene_annotation':null]

mfile = file("${params.options}")
if (!mfile.isFile()) {
   log.info "error: options file not found in '${mfile}'."
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
      log.info "unknown option: '$t[0]' in $params.options"
      exit 1
   }
   options[t[0]] = t[1]
   break
   case 1:
   t = line.split(/\s*,\s*/).toList()
   if (t.size() == 2) {
      if (t[1] =~ /^SRR|^ERR|^DRR/) {
         file_getref << tuple([t[-1]],t[0])
      } else {
         log.info "error: RNA-seq dataset entries must specify 2 files, 2 URL or a GEO reference starting with SRR/ERR/DRR. Error in params entry '$line'"
         exit 1
      }
   } else if (t.size() == 3) {
      if (t[1] =~ /^http|^ftp/ && t[2] =~ /^http|^ftp/) {
         file_getref << tuple([t[-2],t[-1]],t[0])
      } else {
         read1 = file("$t[-2]")
         read2 = file("$t[-1]")
         if (read1.isFile() && read2.isFile()) {
            datasets << tuple([read1,read2],t[0])
         } else {
            log.info "error: RNA-seq files not found, in '$line'. (URLs must start with http/ftp)"
            exit 1
         }
      }
   } else {
      log.info "error: incorrect format: '$line' in $params.options"
      exit 1
   }
   break
   }
}

// Parse options
if (!options['rnalen_mean'] || !options['rnalen_sd']) {
   log.info "error: 'rnalen_mean', 'rnalen_sd' and 'gene_annotation' options must be defined in $params.options before 'datasets:'."
   exit 1
}

/**********************************************************************
 **************************** PARSE  PARAMS ****************************
 **********************************************************************/

gtf_url = Channel.create()
annotation_file_tpm = Channel.create()
annotation_file_exp = Channel.create()

// 0. Get human gene annotations
if (options['gene_annotation']) {
   gene_path = options['gene_annotation']
   f = file("$gene_path")
   // Check whether it is url or file.
   if (gene_path =~ /^http|^ftp/) {   
      gtf_url << gene_path
   } else if (f.isFile()) {
      annotation_file_tpm << f
      annotation_file_tpm.close()
      annotation_file_exp << f
      annotation_file_exp.close()
   } else {
      log.info "error: gene annotation file (gtf) not found: ${f}"
      exit 1
   }
} else {
   log.info "error: gene annotation option 'gene_annotation' not found in metadata! Specify either a http/ftp url or a path to a local file."
   exit 1
}
gtf_url.close()

// 1. Find/Create kallisto index
kallisto_index = Channel.create()
kindex_path    = Channel.create()
cdna_url       = Channel.create()
cdna_file      = Channel.create()

// Find kallisto index
if (params.index == null) {
   log.info "error: '--index' option not specified."
   exit 1
}

// Check kallisto index files.
kindex_ref = file("${params.index}")
if (kindex_ref.isFile()) {
   kallisto_index << kindex_ref
   kallisto_index.close()
} else if (options['transcripts'] != null) {
   f = file(options['transcripts'])
   if (options['trancripts'] =~ /^http|^ftp/) {
      cdna_url << options['transcripts']
   } else if (f.isFile()) {
      cdna_file << f
      cdna_file.close()
   } else {
      log.info "error: gene annotation file (gtf) not found: ${f}"
   }
   kindex_ref.toFile().createNewFile()
   kindex_path << kindex_ref
} else {
   log.info "error: kallisto index not found in '${params.index}' and no 'transcripts' file/url was specified.."
   exit 1
}
kindex_path.close()
cdna_url.close()


/**********************************************************************
 *************************** RNASEQ PIPELINE ***************************
 **********************************************************************/

// Group the same references/URLS to avoid redundant downloads.
file_getref.close()
file_getref.groupTuple().into{gfile_ref}

process getDataset {
   // Process options
   tag "${data[0]}_${data[1]}"
   publishDir path:"${params.out}/datasets/", mode:'symlink'
   // Cluster options
   cpus 1
   memory '2GB'

   input:
   val data from gfile_ref

   output:
   set file('*.fastq.gz'), replicate into datasets

   script:
   replicate = data[1]
   ref = data[0]
   if (ref.size() == 1) {
      if (ref[0] =~ /^SRR|^ERR|^DRR/) {
         """
         fastq-dump --split-files --gzip -A ${ref[0]}
         rm -f ~/ncbi/public/sra/${ref[0]}.sra
         """
      } else {
         log.info "error: only one iPCR dataset specified (must be PE files or sigle GEO reference)!"
      }
   } else if (ref.size() == 2) {
      """
      wget ${ref[0]}
      wget ${ref[1]}
      """
   }
}

process getGTFannotation {
   // Process options
   tag "${url}"
   // Cluster options
   cpus 1
   memory '2GB'

   input:
   val url from gtf_url
   
   output:
   file '*.gtf' into annotation_file_tpm, annotation_file_exp

   script:
   """
   wget ${url} -O - | zcat -f > annotation.gtf
   """
}

process getTranscripts {
   // Process options
   tag ""
   // Cluster options
   cpus 1
   memory '4GB'
   
   input:
   val url from cdna_url

   output:
   file '*.fasta' into cdna_file
   
   script:
   """
   wget ${url} -O - | zcat -f > transcipts.fasta
   """
}


process buildKallistoIndex {
   // Process options
   tag "${kpath}"
   // Cluster options
   cpus 1
   memory '32GB'

   input:
   file kpath from kindex_path
   file cdna from cdna_file

   output:
   file kpath into kallisto_index

   script:
   """
  #!/usr/bin/env bash
  zcat -f ${cdna} | awk 'BEGIN{p=0; IGNORECASE=1;} {
    if(\$1 ~ /^>/) {
      split(\$0,chr,":");
      if (chr[4] ~ /^[MYX0-9][T0-9]*/) { p = 1; }
      else { p = 0; } 
    } 
    if (p == 1) { print \$0; }
  }' > transcripts_chr.fasta
  kallisto index -i ${kpath} transcripts_chr.fasta
  """
}


// 2. Create abundance.tsv files

process processRNASeq {
   // Process options
   tag "${reads}"
   // Cluster options
   cpus 12
   memory '32GB'

   input:
   file kindex from kallisto_index.first()
   set file(reads), sample from datasets
   
   output:
   set file('abundance.tsv'), sample into rnaseq_tsv

   script:
   """
   kallisto quant --bias -l ${options['rnalen_mean']} -s ${options['rnalen_sd']} --single -t ${task.cpus} -i ${kindex} -o . ${reads}
   """
}

// 3. Compute tpm with script
process TSVtoTPM {
   // Process options
   tag "${sample[0]}"
   publishDir path:"${params.out}/tpm/", mode:'symlink'
   // Cluster options
   cpus 1
   memory '4GB'
   
   input:
   set file(tsvfile), sample from rnaseq_tsv
   file humangtf from annotation_file_tpm.first()
   file script from Channel.fromPath("${src_dir}/rnaseq_tpm.py").first()
   output:
   file "*.tpm" into rnaseq_tpm, rnaseq_tpm_fig
   
   script:
   """
   python ${script} --gtf ${humangtf} ${tsvfile} | sort -k1,1 > rnaseq_${sample[0]}.tpm
   """
}

// 4. Merge replicates
process mergeTPMReplicates {
   // Process options
   publishDir path: "${params.out}/", mode:'move'
   // Cluster options
   cpus 1
   memory '4GB'

   input:
   file tpm from rnaseq_tpm.flatten().toList()
   output:
   file 'rnaseq.tpm' into rnaseq_data

   script:
   """
   #!/usr/bin/env Rscript
   files = Sys.glob('*.tpm')
   nf = length(files)
   data = read.table(files[1])
   i = 2
   while (i <= nf) {
      data\$V2 = data\$V2 + read.table(files[i])\$V2
      i = i + 1
   }
   data\$V2 = data\$V2/nf
   write(file='rnaseq.tpm','# gene_id  \tENSEMBL gene id\n# tpm      \ttranscripts per million\n# type     \ttranscript type\n# gene_name\tgene name\n# chr      \tchromosome\n# beg      \ttranscript start nucleotide\n# end      \ttranscript end nucleotide\n# strand   \ttranscript strand + forward, - reverse.\ngene_id\ttpm\ttype\tgene_name\tchr\tbeg\tend\tstrand')
   write.table(data, file='rnaseq.tpm', sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
   """
}

// 5. Figures
process rnaCorPlots {
   // Process options
   publishDir path:"${params.out}/figures/", mode:'move'
   // Cluster options
   cpus 1
   memory '4GB'

   input:
   file libs from Channel.fromPath("${lib_dir}/*.R").flatten().toList()
   file tpm from rnaseq_tpm_fig.flatten().toList()
   output:
   file '*.pdf' into rnaseq_figures

   script:
   """ 
   #!/usr/bin/env Rscript
   source('corplot.R')
   files = Sys.glob('*.tpm')
   nf = length(files)
   if (nf == 0) {
      quit(save='no', status=0)
   }
   rnacol = read.table(files[1])\$V2
   data = matrix(data=NA,nrow=length(rnacol),ncol=nf)
   data[,1] = rnacol
   i = 2
   while (i <= nf) {
      data[,i] = read.table(files[i])\$V2
      i = i + 1
   }
   data <- as.data.frame(data)
   pdf('rnaseq_corplot.pdf', useDingbats=F, width=18, height=18)
   corplot(data)
   dev.off()
   """
}
