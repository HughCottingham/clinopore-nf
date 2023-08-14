nextflow.enable.dsl=2

log.info """\
         CLINOPORE
         ===================================
         reads        : ${params.reads}

         """
         .stripIndent()


// Assembly processes

include { FILTER } from './src/processes/assemble.nf'
include { ASSEMBLE } from './src/processes/assemble.nf'
include { MEDAKA } from './src/processes/assemble.nf'
include { POLYPOLISH } from './src/processes/assemble.nf'
include { POLCA } from './src/processes/assemble.nf'

// Assembly stats processes

include { POLISHING_STATS } from './src/processes/assembly_stats.nf'
include { POLCA_POLISHING_STATS } from './src/processes/assembly_stats.nf'
include { COMBINE_POLISHING_STATS } from './src/processes/assembly_stats.nf'
include { COMBINE_POLCA_POLISHING_STATS } from './src/processes/assembly_stats.nf'
include { COMBINE_ALL_POLISHING_STATS } from './src/processes/assembly_stats.nf'
include { ASSEMBLY_STATS } from './src/processes/assembly_stats.nf'
include { COMBINE_ASSEMBLY_STATS } from './src/processes/assembly_stats.nf'

// Read stats processes

include { RAW_ONT_STATS } from './src/processes/read_stats.nf'
include { COMBINE_RAW_ONT_STATS } from './src/processes/read_stats.nf'
include { FILTERED_ONT_STATS } from './src/processes/read_stats.nf'
include { COMBINE_FILTERED_ONT_STATS } from './src/processes/read_stats.nf'
include { ILLUMINA_STATS } from './src/processes/read_stats.nf'
include { COMBINE_ILLUMINA_STATS } from './src/processes/read_stats.nf'


def check_arguments(params) {
  // Check required input and outputs
  if (! params.reads) {
    exit 1, "ERROR: option 'reads' must be set in nextflow.config or on the command line using --reads"
  }
  if (! params.outdir) {
    exit 1, "ERROR: option 'outdir' must be set in nextflow.config or on the command line using --outdir"
  }
  if (! params.clinopore_env) {
    exit 1, "ERROR: option 'clinopore_env' must be set in nextflow.config or on the command line using --clinopore_env"
  }
  if (! params.polca_env) {
    exit 1, "ERROR: option 'polca_env' must be set in nextflow.config or on the command line using --polca_env"
  }
}

check_arguments(params)

def get_read_prefix_and_type(filepath) {
  // NOTE: nf requires escaping '$'
  regex_se = "^(.+?).fastq(?:.gz)?\$"
  regex_pe = "^(.+?)_R?[12](?:_001)?.fastq(?:.gz)?\$"
  //regex_both = "^(.+?)(_[12])?.fastq(?:.gz)?\$"

  java.util.regex.Matcher matcher;
  String read_type;

  if ((matcher = (filepath.getName() =~ /$regex_pe/))) {
    read_type = 'pe'
  } else if ((matcher = (filepath.getName() =~ /$regex_se/))) {
    read_type = 'se'
  } else {
    exit 1, "ERROR: did not find any readsets with the provided glob: ${filepath}"
  }
  return [read_type, matcher.group(1), filepath]
}



// Get reads and seperate into pe and se channels based on prefix
reads = Channel.fromPath(params.reads).ifEmpty {
    exit 1, "ERROR: did not find any read files with '${params.reads}'"
  }.map {
    get_read_prefix_and_type(it)
  }.branch {
    paired: it[0] == 'pe'
    single: it[0] == 'se'
}
reads_se = reads.single.map { it[1..-1] }.groupTuple()
reads_pe = reads.paired.map { it[1..-1] }.groupTuple()
// Check that we have the expected number of reads for each prefix in pe and se channels and flatten tuple
reads_pe = reads_pe.map {
  if (it[1].size() != 2) {
    exit 1, "ERROR: didn't get exactly two readsets prefixed with ${it[0]}:\n${it[1]}"
  }
  [it[0], *it[1]]
}
reads_se = reads_se.map {
  if (it[1].size() != 1) {
    exit 1, "ERROR: didn't get exactly one readset prefixed with ${it[0]}:\n${it[1]}"
  }
  [it[0], *it[1]]
}
//reads_se.view()
//reads_pe.view()

//combined_ch = reads_se
    //.combine(reads_pe, by: 0)



def check_host(workflow) {
  // Do not run on MASSIVE unless user specifies profile to use to avoid inadvertently using a local executor
  massive_hostnames = ['m3-login1', 'm3-login2']
  on_massive = massive_hostnames.contains(InetAddress.getLocalHost().getHostName())
  profile_explicit = workflow.commandLine.tokenize(' ').contains('-profile')
  if (on_massive && ! profile_explicit) {
    exit 1, "ERROR: to run on MASSIVE you must explicitly set -profile massive"
  }
}


// For an optional stage param variable, check that it is either a Boolean or String
// If it is a string and either 'true' or 'false', return the boolean equivalent

def check_boolean_option(option, name) {
  if (option.getClass() == java.lang.Boolean) {
    return option
  } else if (option.getClass() == java.lang.String) {
    if (option.toLowerCase() == 'true') {
      return true
    } else if (option.toLowerCase() == 'false') {
      return false
    }
  }
  exit 1, "ERROR: ${name} option must be true or false"
}

filter_reads = check_boolean_option(params.filter_reads, 'filter_reads')
long_read_polish = check_boolean_option(params.long_read_polish, 'long_read_polish')
short_read_polish = check_boolean_option(params.short_read_polish, 'short_read_polish')
read_qc = check_boolean_option(params.read_qc, 'read_qc')
polishing_stats = check_boolean_option(params.polishing_stats, 'polishing_stats')


workflow qc_and_assemble {
  take:
  reads_se
  main:
  if( filter_reads ) {
    input_ch=FILTER(reads_se)
  }
  else {
    input_ch=reads_se
  }
  if( read_qc && filter_reads ) {
    raw_stats_ch=RAW_ONT_STATS(reads_se)
    COMBINE_RAW_ONT_STATS(raw_stats_ch.collect())
    filtered_stats_ch=FILTERED_ONT_STATS(input_ch)
    COMBINE_FILTERED_ONT_STATS(filtered_stats_ch.collect())
  }
  else if( read_qc ) {
    raw_stats_ch=RAW_ONT_STATS(reads_se)
    COMBINE_RAW_ONT_STATS(raw_stats_ch.collect())
  }
  assembled_ch=ASSEMBLE(input_ch)
  emit:
  assembly_ch=ASSEMBLE.out.assembly_fasta
  reads_ch=input_ch
}

workflow long_read_polishing {
  take:
  assembly_ch
  reads_ch
  main:
  assembly_reads_ch=assembly_ch
  .combine(reads_ch, by:0)
  medaka_polished_ch=MEDAKA(assembly_reads_ch)
  if( polishing_stats ) {
    polishing_stats_ch=POLISHING_STATS(assembly_ch,MEDAKA.out.assembly_fasta,'medaka')
    combined=COMBINE_POLISHING_STATS(polishing_stats_ch.collect(),'medaka')
  }
  else {
    combined=medaka_polished_ch
  }
  emit:
  assembly_ch=medaka_polished_ch.assembly_fasta
  medaka_stats_ch=combined
}

workflow short_read_polishing {
  take:
  reads_pe
  latest_assembly_ch
  flye_ch
  main:
  if( read_qc ) {
    fastqc_ch = ILLUMINA_STATS(reads_pe)
    COMBINE_ILLUMINA_STATS(fastqc_ch.collect())
  }
  assembly_reads_ch=latest_assembly_ch
  .combine(reads_pe, by: 0)
  polypolish_input_ch=assembly_reads_ch
  .combine(flye_ch,by:0)
  polypolished_ch=POLYPOLISH(polypolish_input_ch)
  if( polishing_stats ) {
    polypolish_stats_ch=POLISHING_STATS(latest_assembly_ch,polypolished_ch,'polypolish')
    polypolish_stats_combined=COMBINE_POLISHING_STATS(polypolish_stats_ch.collect(),'polypolish')
  }
  else {
    polypolish_stats_combined=polypolished_ch
  }
  polypolish_assembly_reads_ch=polypolished_ch
  .combine(reads_pe, by: 0)
  polca_input_ch=polypolish_assembly_reads_ch
  .combine(flye_ch,by:0)  
  polca_ch=POLCA(polca_input_ch)
  if( polishing_stats ) {
    polca_stats_ch=POLCA_POLISHING_STATS(latest_assembly_ch,polca_ch,'polca')
    polca_stats_combined=COMBINE_POLCA_POLISHING_STATS(polca_stats_ch.collect(),'polca')
  }
  else {
    polca_stats_combined=polca_ch
  }
  emit:
  assembly_ch=POLCA.out.assembly_fasta
  polypolish_stats_combined_ch=polypolish_stats_combined
  polca_stats_combined_ch=polca_stats_combined
}

workflow assembly_stats {
  take:
  assembly_ch
  main:
  stats_ch = ASSEMBLY_STATS(assembly_ch)
  COMBINE_ASSEMBLY_STATS(stats_ch.collect())
}

workflow all_polishing_stats {
  take:
  medaka_stats_ch
  polypolish_stats_combined_ch
  polca_stats_combined_ch
  main:
  all1=medaka_stats_ch
  .combine(polypolish_stats_combined_ch)
  all_polishing_stats=all1
  .combine(polca_stats_combined_ch)
  COMBINE_ALL_POLISHING_STATS(all_polishing_stats)
}

workflow short_polishing_stats {
  take:
  polypolish_stats_combined_ch
  polca_stats_combined_ch
  main:
  all_short_polish_stats=polypolish_stats_combined_ch
  .combine(polca_stats_combined_ch)
  COMBINE_ALL_POLISHING_STATS(all_short_polish_stats)
}

workflow long_polishing_stats {
  take:
  medaka_stats_ch
  main:
  COMBINE_ALL_POLISHING_STATS(medaka_stats_ch)
}


//////IMPLICIT
workflow {
  if (long_read_polish && short_read_polish && polishing_stats ) {
    qc_and_assemble(reads_se)
    long_read_polishing(qc_and_assemble.out)
    short_read_polishing(reads_pe,long_read_polishing.out.assembly_ch,qc_and_assemble.out.assembly_ch)
    assembly_stats(short_read_polishing.out.assembly_ch)
    all_polishing_stats(long_read_polishing.out.medaka_stats_ch,short_read_polishing.out.polypolish_stats_combined_ch,short_read_polishing.out.polca_stats_combined_ch)

  }
  else if (long_read_polish && short_read_polish ) {
    qc_and_assemble(reads_se)
    long_read_polishing(qc_and_assemble.out)
    short_read_polishing(reads_pe,long_read_polishing.out.assembly_ch,qc_and_assemble.out.assembly_ch)
    assembly_stats(short_read_polishing.out.assembly_ch)
  }

  else if ( long_read_polish && polishing_stats ) {
    qc_and_assemble(reads_se)
    long_read_polishing(qc_and_assemble.out)
    assembly_stats(long_read_polishing.out.assembly_ch)
    long_polishing_stats(long_read_polishing.out.medaka_stats_ch)
  }

  else if ( short_read_polish && polishing_stats ) {
    qc_and_assemble(reads_se)
    short_read_polishing(reads_pe,qc_and_assemble.out.assembly_ch,qc_and_assemble.out.assembly_ch)
    assembly_stats(short_read_polishing.out.assembly_ch)
    short_polishing_stats(short_read_polishing.out.polypolish_stats_combined_ch,short_read_polishing.out.polca_stats_combined_ch)
  }

  else if ( short_read_polish ) {
    qc_and_assemble(reads_se)
    short_read_polishing(reads_pe,qc_and_assemble.out.assembly_ch,qc_and_assemble.out.assembly_ch)
    assembly_stats(short_read_polishing.out.assembly_ch)
  }

  else if ( long_read_polish ) {
    qc_and_assemble(reads_se)
    long_read_polishing(qc_and_assemble.out)
    assembly_stats(long_read_polishing.out.assembly_ch)
  }
  
  else {
    qc_and_assemble(reads_se)
    assembly_stats(qc_and_assemble.out.assembly_ch)
  }
}

