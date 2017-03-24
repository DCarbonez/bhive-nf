public boolean isInteger( String input ) {
    try {
        Integer.parseInt( input );
        return true;
    }
    catch( NumberFormatException e ) {
        return false;
    }
}

def print_help = {
    log.info ''
    log.info 'GenomeArchitecture BHIVE pipeline'
    log.info '---------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info '    bhive.nf [dataset options] [genome options] [Other options]...'
    log.info ''
    log.info 'Dataset options:'
    log.info ''
    log.info '  Not implemented yet.'
}

script_dir = "src"
lib_dir = "lib"

extract_barcode_src = Channel.fromPath("${script_dir}/mapping_extract_barcodes.py")
barcode_seqnames_src = Channel.fromPath("${script_dir}/mapping_barcode_seqnames.py")
best_assignment_src = Channel.fromPath("${script_dir}/mapping_best_assignment.py")
rnaseq_tpm_src = Channel.fromPath("${script_dir}/rnaseq_tpm.py")
rlibs_src = Channel.fromPath("${lib_dir}/*.R")


// Set defaults
params.mapcpus     = 12
params.out         = 'out'
params.colorspace  = false
params.help        = false
params.metadata    = null
params.datadir     = null
params.index       = null
params.kindex      = null
params.assembly    = null
params.genome      = null


if (params.help) {
   print_help()
   exit 0
}

/**********************************************************************
**************************** PARSE  PARAMS ****************************
**********************************************************************/


bwa_index   = Channel.create()
index_files = Channel.create()
// 0. Find BWA index
if (params.index) {
   // Check bwa index files.
   fasta_ref = file("${params.index}.bwt")
   index_ref = file("${params.index}")
   if (fasta_ref.isFile()) {
      if (!index_ref.isFile()) {
         index_ref.toFile().createNewFile()
      }
      bwa_index << file("${params.index}")
      index_files << file("${index_ref}.{bwt,amb,ann,pac,sa}")
   } else if (index_ref.getExtension() in ['bwt','amb','ann','pac','sa'] && index_ref.isFile()) {
      index_str = (index_ref.getParent().equals(null) ? './' : index_ref.getParent()) + index_ref.getBaseName()
      index_ref = file(index_str)
      if (!index_ref.isFile()) {
         index_ref.toFile().createNewFile()
      }
      bwa_index << index_ref
      index_files << file("${index_ref}.{bwt,amb,ann,pac,sa}")
   } else {
      log.info "error: BWA index not found in '${params.index}'."
      exit 1
   }
} else {
   log.info "error: '--index' option not specified."
   exit 1
}
bwa_index.close()
index_files.close()


/**********************************************************************
*************************** PARSE  METADATA ***************************
**********************************************************************/
/*
** Metadata format:
** 
**  BFS (Barcode flanking sequence)
**  LTR (HIV LTR sequence, must be present in reverse read)
**  RFS (Restriction enzyme flanking sequence, in HIV construct)
**  DIST (Barcode clustering distance)
**  MAPQ (Minimum mapping quality, 20)
**  INTQ (Minimum assignment score, 10)
** 
** ipcr:
** biological replicate,filename/URL/{SRR,ERR,DRR} reference
** 
** dna:
** biological replicate,technical replicate,DNA index,file/URL/{SRR,ERR,DRR} reference
** 
** rna:
** biological replicate,technical replicate,DNA index,file/URL/{SRR,ERR,DRR} reference
*/

def options = ['bfs':null, 'ltr':null, 'rfs':null, 'dist':null, 'mapq':null, 'intq':null, 'gene_annotation':null, 'transcripts_url':null, 'rnalen_mean':null, 'rnalen_sd':null]

if (!params.metadata) {
   log.info 'error: metadata file not specified!'
   exit 1
}

// Open file.
mfile = file("$params.metadata")
// Check whether file exists.
if (!mfile.isFile()) {
   log.info "error: metadata file not found ($params.metadata)"
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
   } else if (line =~ /^ipcr:/) {
      status = 1
      return null
   } else if (line =~ /^dna:/) {
      status = 2
      return null
   } else if (line =~ /^rna:/) {
      status = 3
      return null
   } else if (line =~ /^rnaseq:/) {
      status = 4
      return null
   }
   switch (status) {
   case 0:
     t = line.split(/\s*=\s*/).toList()
     if (!options.containsKey(t[0])) {
        log.info "unknown option: '$t[0]' in $params.metadata"
        exit 1
     }
     options[t[0]] = t[1]
     break
   case 1:
     t = line.split(/\s*,\s*/).toList()
     t.add(0, 'ipcr')
     if (t.size() == 3) {
        if (t[2] =~ /^SRR|^ERR|^DRR/) {
           file_getref << tuple([t[-1]],t[0..-2])
        } else {
           log.info "error: iPCR entries must specify 2 files, 2 URL or a GEO reference starting with SRR/ERR/DRR. Error in entry '$line'"
           exit 1
        }
     } else if (t.size() == 4) {
        if (t[2] =~ /^http|^ftp/ && t[3] =~ /^http|^ftp/) {
           file_getref << tuple([t[-2],t[-1]],t[0..-3])
        } else {
           read1 = file("$t[-2]")
           read2 = file("$t[-1]")
           if (read1.isFile() && read2.isFile()) {
              datasets << tuple([read1,read2],t[0..-3])
           } else {
              log.info "error: iPCR files not found, in '$line'. (URLs must start with http/ftp)"
              exit 1
           }
        }
     } else {
        log.info "error: incorrect format of iPCR entry '$line' in $params.metadata"
        exit 1
     }
     break
   case 2..3:
     t = line.split(/\s*,\s*/).toList()
     t.add(0, status == 2 ? 'dna' : 'rna')
     if (t.size() != 4) {
        log.info "error: incorrect format of ${status == 2 ? 'DNA' : 'RNA'} entry '$line' in $params.metadata"
        exit 1
     } else {
        if (t[3] =~ /^SRR|^ERR|^DRR|^http|^ftp/) {
           file_getref << tuple([t[-1]],t[0..-2])
        } else {
           f = file("$t[-1]")
           if (f.isFile()) {
              datasets << tuple([f],t[0..-2])
           } else {
              log.info "error: ${status==2 ? 'DNA' : 'RNA'} file not found, in '$line'. (URLs must start with http/ftp)"
           }
        }
     }
     break
   case 4:
     t = line.split(/\s*,\s*/).toList()
     t.add(0,'rnaseq')
     if (t.size() != 3) {
        log.info "error: incorrect format for rnaseq entry '$line' in $params.metadata"
        exit 1
     } else {
        if (t[2] =~ /^SRR|^ERR|^DRR|^http|^ftp/) {
           file_getref << tuple([t[-1]],t[0..-2])
        } else {
           f = file("$t[-1]")
           if (f.isFile()) {
              datasets << tuple([f],t[0..-2])
           } else {
              log.info "error: rnaseq file $f not found, in '$line'. (URLs must start with http/ftp)"
           }
        }
     }
     break
   }
}

// Parse options
if (!options['bfs'] || !options['ltr'] || !options['dist'] || !options['gene_annotation']) {
   log.info "error: 'bfs', 'ltr' and 'dist' options must be defined in $params.metadata before 'ipcr:', 'dna:' and 'rna:'"
   exit 1
}
if (!isInteger(options['dist'])) {
   log.info "error: 'dist' value provided in $params.metadata must be a valid integer."
}

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
   set file('*.fastq.gz'), info into datasets
   script:
   info = data[1]
   ref = data[0]
   if (ref.size() == 1) {
      if (ref[0] =~ /^SRR|^ERR|^DRR/) {
         """
         fastq-dump --split-files --gzip -A ${ref[0]}
         rm -f ~/ncbi/public/sra/${ref[0]}.sra
         """
      }
      else if (ref[0] =~ /^http|^ftp/) {
         """
         wget ${ref[0]}
         """
      }
   } else if (ref.size() == 2) {
      """
      wget ${ref[0]}
      wget ${ref[1]}
      """
   }
}

// 1.2 Demultiplex datasets.

maps = Channel.create()
pcrs = Channel.create()
rseq = Channel.create()

datasets.subscribe onNext: { files, info ->
   for (sample in info) {
      switch (sample[0]) {
      case 'ipcr':
      maps << [files, sample]
      break

      case 'dna':
      pcrs << [files, sample]
      break

      case 'rna':
      pcrs << [files, sample]
      break
      
      case 'rnaseq':
      rseq << [files, sample]
      break
      }
   }
}, onComplete: {
   maps.close()
   pcrs.close()
}


/**********************************************************************
*************************** RNASEQ PIPELINE ***************************
**********************************************************************/

gtf_url = Channel.create()
annotation_file = Channel.create()

// 0. Get human gene annotations
if (options['gene_annotation']) {
   gene_path = options['gene_annotation']
   f = file("$gene_path")
   // Check whether it is url or file.
   if (gene_path =~ /^http|^ftp/) {   
      gtf_url << gene_path
   } else if (f.isFile()) {
      annotation_file << f
      annotation_file.close()
   } else {
      log.info "error: gene annotation file (gtf) not found: ${f}"
   }
} else {
   log.info "error: gene annotation option 'gene_annotation' not found in metadata! Specify either a http/ftp url or a path to a local file."
}
gtf_url.close()

process getGTFannotation {
   // Process options
   tag "${url}"
   // Cluster options
   cpus 1
   memory '2GB'

   input:
   val url from gtf_url
   
   output:
   file '*.gtf' into annotation_file

   script:
   """
   wget ${url} -O annotation.gtf
   """
}

// 1. Find/Create kallisto index
kallisto_index = Channel.create()
cdna_url       = Channel.create()

// Find kallisto index
if (params.kindex) {
   // Check bwa index files.
   kindex_ref = file("${params.kindex}")
   if (kindex_ref.isFile()) {
      kallisto_index << kindex_ref
      kallisto_index.close()
   } else if (options['transcripts_url'] != null) {
      kindex_ref.toFile().createNewFile()
      kindex_path << kindex_ref
      cdna_url << options['transcripts_url']
   } else {
      log.info "error: kallisto index not found in '${params.kindex}'."
      exit 1
   }
   cdna_url.close()
} else {
   log.info "error: '--kindex' option not specified."
   exit 1
}

process buildKallistoIndex {
   // Process options
   tag "${kpath}"
   // Cluster options
   cpus 1
   memory '32GB'

   input:
   file kpath from kindex_path
   val  url   from cdna_url
   output:
   file kpath into kallisto_index
   script:
     """
     wget ${url} -O - | zcat -f | awk 'BEGIN{p=0; IGNORECASE=1;} {
       if($1 ~ /^>/) {
         split($0,chr,":");
         if (chr[4] ~ /^[MYX0-9][T0-9]*/) { p = 1; }
         else { p = 0; } 
       } 
       if (p == 1) { print $0; }
     }' > transcripts.fasta
     kalisto index -i ${kindex_path} transcripts.fasta
     """
}

// 2. Create abundance.tsv files
process processRNASeq {
   // Process options
   tag "${reads}"
   // Cluster options
   cpus 12
   memory '2GB'

   input:
   file kindex from kallisto_index.first()
   set file(reads), sample from rseq
   
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
   tag "${sample[1]}"
   publishDir path:"${params.out}/expression/rnaseq/", mode:'symlink'
   // Cluster options
   cpus 1
   memory '4GB'
   
   input:
   set file(tsvfile), sample from rnaseq_tsv
   file humangtf from annotation_file
   file script from rnaseq_tpm_src.first()
   output:
   file ('*.tpm') into rnaseq_tpm

   script:
   """
   python ${script} --gtf ${humangtf} ${tsvfile} | sort -k1,1 > rnaseq_${sample[1]}.tpm
   """
}

// 4. Figures
process rnaCorPlots {
   // Process options
   tag ""
   publishDir path:"${params.out}/expression/figures/", mode:'move'
   // Cluster options
   cpus 1
   memory '4GB'

   input:
   file libs from rlibs_src.first()
   rnaseq_tpm.flatten().toList()
   output:
   file '*.pdf' into rnaseq_figures

   script:
   """ 
   #!/usr/bin/env 
   source('corplot.R')
   files = Sys.glob('*.tpm')
   nf = length(files)
   if (nf == 0) {
      quit(save='no', status=0)
   }
   rnacol = read.table(files[1])$V2
   data = matrix(data=NA,nrow=length(rnacol),ncol=nf)
   data[,1] = rnacol
   i = 2
   while (i <= nf) {
      data[,i] = read.table(files[i])$V2
   }
   data <- as.data.frame(data)
   pdf('rnaseq_corplot.pdf', useDingbats=F, width=18, height=18)
   corplot(data)
   dev.off()
   """
}

/**********************************************************************
**************************** IPCR PIPELINE ****************************
**********************************************************************/


process extractBarcodes {
   // Process options
   tag "${sample[1]}"
   publishDir path:"${params.out}/map/barcodes_raw/", mode:'symlink'
   
   // Custer options
   cpus 1
   memory '4GB'

   input:
   set file(ipcr), sample from maps
   file script from extract_barcode_src.first()

   output:
   set file("*.bcd"), file("*.fastq"), sample into raw_barcodes
   
   script:
   rfs_arg = options['rfs'] ? "-r ${options['rfs']}" : ""
   py_options = "-b ${options['bfs']} -l ${options['ltr']} ${rfs_arg}"
   """
   python ${script} ${py_options} --barcodes ipcr_${sample[1]}.bcd --reads ipcr_${sample[1]}.fastq ${ipcr[0]} ${ipcr[1]}
   """

}

process clusterBarcodes {
   // Process options
   tag "$sample[1]"
   publishDir path:"${params.out}/map/barcodes_cluster/", mode:'symlink'
   
   // Cluster options
   cpus 12
   memory '64GB'

   input:
   set file(barcodes), file(reads), sample from raw_barcodes
   file script from barcode_seqnames_src.first()
   
   output:
   set file("ipcr_integ_site_*.fastq"), sample into loci_map
   
   script:
   """
   starcode -d ${options['dist']} --seq-id -t ${task.cpus} ${barcodes} > ipcr_${sample[1]}.stc
   python ${script} ipcr_${sample[1]}.stc ${reads} > ipcr_integ_site_${sample[1]}.fastq
   """
}

process mapIntegs {
   // Process options
   tag "$sample[1]"
   publishDir path:"${params.out}/map/mapping/", mode:'symlink'

   // Cluster options
   cpus 12
   memory '32GB'

   input:
   set file(integ_sites), sample from loci_map
   file index_path from bwa_index.first()
   file index_files from index_files.first()

   output:
   set file("*.bam"), sample into raw_integs
   
   script:
   cat_head = integ_sites.getExtension() == 'gz' ? "zcat ${integ_sites}" : "cat ${integ_sites}"

   """
    OPTS_36='-k18 -B3 -O5 -T28'
    OPTS_26='-k17 -r1.3 -B2 -O4 -T22'
    OPTS_40=''
    OPTS_50=''
    SEQLEN=\$((\$(${cat_head} | head -2 | tail -1 | wc -c) - 1));
    if [ \$SEQLEN -le 30 ]; then \
      OPTS=\$OPTS_26; \
    elif [ \$SEQLEN -le 39 ]; then \
      OPTS=\$OPTS_36; \
    elif [ \$SEQLEN -le 46 ]; then \
      OPTS=\$OPTS_40; \
    else \
      OPTS=\$OPTS_50; \
    fi;
    echo \$SEQLEN
    echo \$OPTS
    echo ${index_path}
    bwa mem -t ${task.cpus} \$OPTS ${index_path} ${integ_sites} | samtools view -bS - > ipcr_${sample[1]}.bam
   """
}

process filterIntegNoise {
   // Process options
   tag "$sample[1]"
   publishDir path:"${params.out}/map/integs", mode:'symlink'

   // Cluster options
   cpus 1
   memory '4GB'
   
   input:
   set file(raw_integs), sample from raw_integs
   file script from best_assignment_src.first()

   output:
   file "*.txt" into integ_list
   
   script:
   py_args = (options['mapq'] ? "-q ${options['mapq']} " : "") + (options['intq'] ? "--min-score ${options['intq']}" : "")
   """
   python ${script} ${py_args} <(samtools view ${raw_integs}) > ipcr_integs_${sample[1]}.txt
   """
}

process mergeIntegs {
   // Process options
   tag "${finteg.size()} files"
   publishDir path:"${params.out}/map/", mode:'move'

   // Cluster options
   cpus 1
   memory '4GB'
   
   input:
   file(finteg) from integ_list.flatten().toSortedList()

   output:
   file "ipcr_integs.txt" into integ_file
   
   script:
"""
#!/usr/bin/env python
files = ${finteg.collect{'\"'+it+'\"'}}
for fname in files:
  rep = fname.split('_')[2].split('.')[0]
  with open(fname,'r') as fin, open('ipcr_integs.txt','a') as fout:
    for line in fin:
      fout.write(line.rstrip()+'\\t'+rep+'\\n')
"""
}


/**********************************************************************
************************* EXPRESSION PIPELINE *************************
**********************************************************************/


process extractPCRBarcodes {
   // Process options
   tag "${sample[0]}_${sample[1]}_${sample[2]}"
   publishDir path:"${params.out}/expression/barcodes_raw/", mode:'symlink'
   
   // Custer options
   cpus 1
   memory '4GB'

   input:
   set file(reads), sample from pcrs

   output:
   set file("*.bcd"), sample into pcr_barcodes
   
   script:
   outf = "${sample[0]}_${sample[1]}_${sample[2]}.bcd"
"""
#!/usr/bin/env python
import seeq
import gzip
MIN_BRCD_LEN = 15

# Helper functions.
def getbarcode(matcher, txt):
   return matcher.matchPrefix(txt, False) or ''

if __name__ == '__main__':
   fname = '${reads}'
   fout = '${outf}'
   outp = open(fout, 'w+')
   T7seq = '${options['bfs']}'
   T7 = seeq.compile(T7seq, int(max(1,round(0.2*len(T7seq)))))
   with gzip.open(fname) as f:
      # Read fastq file.
      for lineno,line in enumerate(f):
         if lineno % 4 != 1: continue
         barcode = getbarcode(T7, line)
         if len(barcode) < MIN_BRCD_LEN: continue
         outp.write(barcode + '\\n')

   outp.close()

"""
}

process clusterPCRBarcodes {
   // Process options
   tag "$sample[1]"
   publishDir path:"${params.out}/expression/barcodes_cluster/", mode:'symlink'
   
   // Cluster options
   cpus 12
   memory '64GB'

   input:
   set file(barcodes), sample from pcr_barcodes
   
   output:
   set file("*.stc"), sample into pcr_clusters
   
   script:
   """
   starcode -d ${options['dist']} -qt ${task.cpus} ${barcodes} > ${sample[0]}_${sample[1]}_${sample[2]}.stc
   """
}
