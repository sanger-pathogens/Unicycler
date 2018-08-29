# Unicycler-polish

Unicycler-polish is a program to repeatedly polish a completed assembly using all available reads. It can be given Illumina reads, long reads or (ideally) both. When both Illumina and long reads are available, Unicycler-polish can fix assembly errors repetitive parts of the genome which cannot be polished by short reads alone.

It should be considered somewhat experimental, so use with caution!


### Requirements

* If polishing with Illumina reads: [Pilon](https://github.com/broadinstitute/pilon/wiki), Java, [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/), [Samtools](http://www.htslib.org/) (version 1.0 or later)
* If polishing with PacBio reads: [pbalign](https://github.com/PacificBiosciences/pbalign), [BLASR](https://github.com/PacificBiosciences/blasr), [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus)
    * PacBio software is most easily installed using [pitchfork](https://github.com/PacificBiosciences/pitchfork).
* If polishing with Illumina and PacBio reads: all of the above dependencies plus [ALE](https://github.com/sc932/ALE).
* If polishing with both Illumina and long reads (e.g. Nanopore): all of the Illumina polishing dependencies, [Racon](https://github.com/isovic/racon), [FreeBayes](https://github.com/ekg/freebayes) and [ALE](https://github.com/sc932/ALE).


### Process

Unicycler polish uses an exhaustive iterative process that is time-consuming but can be necessary to resolve the sequence in repeat regions. For example, consider a genome with two very similar regions, A and B, and there are assembly errors in both. Polishing is initially difficult because the errors may cause reads which should map to A to instead map to B and vice versa. However, after some of these errors are fixed, more reads will map to their correct locations, allowing for more errors to be fixes, allowing more reads to map correctly, etc.

1. If Illumina reads are available:
    1. Run [Pilon](https://github.com/broadinstitute/pilon/wiki) in 'bases' mode (substitutions and small indels). If any changes were suggested, apply them and repeat this step.
    2. Run Pilon in 'local' mode (larger variants), and assess each change with ALE. If any variant improves the ALE score, apply it and go back to step 1-i.
2. If long reads are available:
    1. Run [GenomicConsensus](https://github.com/PacificBiosciences/GenomicConsensus)/[Racon](https://github.com/isovic/racon) and gather all suggested small changes.
    2. Use [FreeBayes](https://github.com/ekg/freebayes) to assess each long read-suggested change by looking for ambiguity in the Illumina read mapping. If any were found, apply them and go back to step 2-i.
3. If Illumina reads are available:
    1. Execute step 1 again.
    2. Run Pilon/GenomicConsensus/Racon again (all that apply) and assess each suggested variant with [ALE](https://github.com/sc932/ALE). If any improves the ALE score, apply it and repeat this step.


### Example commands

__Polishing with only Illumina reads:__<br>
`unicycler_polish -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -a assembly.fasta`

__Polishing with only PacBio reads:__<br>
`unicycler_polish --pb_bax path/to/*bax.h5 -a assembly.fasta`

__Hybrid read set (Illumina and PacBio) polishing:__<br>
`unicycler_polish -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz --pb_bax *bax.h5 -a assembly.fasta`

__Hybrid read set (Illumina and Nanopore) polishing:__<br>
`unicycler_polish -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz --long_reads nanopore.fastq.gz -a assembly.fasta`



### Usage

```
usage: unicycler_polish [-h] -a ASSEMBLY [-1 SHORT1] [-2 SHORT2] [--pb_bax PB_BAX [PB_BAX ...]] [--pb_bam PB_BAM] [--pb_fasta PB_FASTA]
                        [--long_reads LONG_READS] [--min_insert MIN_INSERT] [--max_insert MAX_INSERT]
                        [--min_align_length MIN_ALIGN_LENGTH] [--homopolymer HOMOPOLYMER] [--large LARGE] [--illumina_alt ILLUMINA_ALT]
                        [--freebayes_qual_cutoff FREEBAYES_QUAL_CUTOFF] [--threads THREADS] [--verbosity VERBOSITY] [--samtools SAMTOOLS]
                        [--bowtie2 BOWTIE2] [--freebayes FREEBAYES] [--pitchfork PITCHFORK] [--bax2bam BAX2BAM] [--pbalign PBALIGN]
                        [--arrow ARROW] [--pilon PILON] [--java JAVA] [--ale ALE] [--racon RACON] [--minimap MINIMAP] [--nucmer NUCMER]
                        [--showsnps SHOWSNPS]

Unicycler polish - hybrid assembly polishing

optional arguments:
  -h, --help                            show this help message and exit

Assembly:
  -a ASSEMBLY, --assembly ASSEMBLY      Input assembly to be polished

Short reads:
  To polish with short reads (using Pilon), provide two FASTQ files of paired-end reads

  -1 SHORT1, --short1 SHORT1            FASTQ file of short reads (first reads in each pair)
  -2 SHORT2, --short2 SHORT2            FASTQ file of short reads (second reads in each pair)

PacBio reads:
  To polish with PacBio reads (using Arrow), provide one of the following

  --pb_bax PB_BAX [PB_BAX ...]          PacBio raw bax.h5 read files
  --pb_bam PB_BAM                       PacBio BAM read file
  --pb_fasta PB_FASTA                   FASTA file of PacBio reads

Generic long reads:
  To polish with generic long reads, provide the following

  --long_reads LONG_READS               FASTQ/FASTA file of long reads

Polishing settings:
  Various settings for polishing behaviour (defaults should work well in most cases)

  --no_fix_local                        do not fix local misassemblies (default: False)
  --min_insert MIN_INSERT               minimum valid short read insert size (default: auto)
  --max_insert MAX_INSERT               maximum valid short read insert size (default: auto)
  --min_align_length MIN_ALIGN_LENGTH   Minimum long read alignment length (default: 1000)
  --homopolymer HOMOPOLYMER             Long read polish changes to a homopolymer of this length or greater will be ignored (default: 4)
  --large LARGE                         Variants of this size or greater will be assess as large variants (default: 10)
  --illumina_alt ILLUMINA_ALT           When assessing long read changes with short read alignments, a variant will only be applied if the
                                        alternative occurrences in the short read alignments exceed this percentage (default: 5)
  --freebayes_qual_cutoff FREEBAYES_QUAL_CUTOFF
                                        Reject Pilon substitutions from long reads if the FreeBayes quality is less than this value
                                        (default: 10.0)

Other settings:
  --threads THREADS                     CPU threads to use in alignment and consensus (default: number of CPUs)
  --verbosity VERBOSITY                 Level of stdout information (0 to 3, default: 2)
                                          0 = no stdout, 1 = basic progress indicators, 2 = extra info, 3 = debugging info

Tool locations:
  If these required tools are not available in your PATH variable, specify their location here (depending on which input reads are used,
  some of these tools may not be required)

  --samtools SAMTOOLS                   path to samtools executable (default: samtools)
  --bowtie2 BOWTIE2                     path to bowtie2 executable (default: bowtie2)
  --freebayes FREEBAYES                 path to freebayes executable (default: freebayes)
  --pitchfork PITCHFORK                 Path to Pitchfork installation of PacBio tools (should contain bin and lib directories) (default:
                                        )
  --bax2bam BAX2BAM                     path to bax2bam executable (default: bax2bam)
  --pbalign PBALIGN                     path to pbalign executable (default: pbalign)
  --arrow ARROW                         path to arrow executable (default: arrow)
  --pilon PILON                         path to pilon jar file (default: pilon*.jar)
  --java JAVA                           path to java executable (default: java)
  --ale ALE                             path to ALE executable (default: ALE)
  --racon RACON                         path to racon executable (default: racon)
  --minimap MINIMAP                     path to miniasm executable (default: minimap)
  --nucmer NUCMER                       path to nucmer executable (default: nucmer)
  --showsnps SHOWSNPS                   path to show-snps executable (default: show-snps)
```
