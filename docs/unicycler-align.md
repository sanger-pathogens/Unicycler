# Unicycler-align

Unicycler's algorithm for sensitive semi-global alignment is available as a stand-alone alignment tool with the command `unicycler_align`.


### Semi-global alignment

Semi-global alignment (a.k.a. glocal, overlap or free end-gap alignment) will not clip an alignment until one of the two sequences ends. This can be where one sequence is contained within the other or where the two sequences overlap:
```
  TAGAA        GTGCCGGAACA         GGCCACAC     AGTAAGAT
  |||||          |||||||           |||||           |||||
ACTAGAACG        GCCGGAA       GGCTGGCCA           AAGATCTTG
```

In contrast, local alignment will align only the best matching parts, clipping the alignment where the quality becomes poor:
```
      CGAACAGCATACTTG
          ||||||||
ACGTCAGACTCAGCATACGCATCTAGA
```

Semi-global alignment is appropriate when there are no structural differences between the query and reference sequences. For example, when you have a short read assembly graph and long reads from the same bacterial isolate (as is the case in the Unicycler pipeline). In this scenario, there may be small scale differences (due to read errors) but no large scale differences, and semi-global alignment is ideal.


### Versus local alignment

Semi-global alignment is probably not appropriate for mapping reads to a more distant reference genome. It does not cope with points of structural variation between the sample and the reference. For example, if the sample had a deletion relative to the reference, a read spanning that deletion would align poorly with semi-global alignment:
```
read:            AACACTAAACTTAGTCCCAA
                 |||||||||||  |   | |    
reference: GATCCCAACACTAAACTCTGGGGCGAACGGCGTAGTCCCAAGAGT
```

Local alignment (which can align only part of the read) would be more appropriate:
```
read:            AACACTAAACT               TAGTCCCAA
                 |||||||||||               |||||||||
reference: GATCCCAACACTAAACTCTGGGGCGAACGGCGTAGTCCCAAGAGT
```
Try [BWA-MEM](http://bio-bwa.sourceforge.net/), [LAST](http://last.cbrc.jp/) or [BLASR](https://github.com/PacificBiosciences/blasr) if you need a local alignment tool.


### Example commands

__Regular alignment:__<br>
`unicycler_align --reads queries.fastq --ref target.fasta --sam output.sam`

__Very sensitive (and slow) alignment:__<br>
`unicycler_align --reads queries.fastq --ref target.fasta --sam output.sam --sensitivity_level 3`

__Setting some additional thresholds:__<br>
`unicycler_align --reads queries.fastq --ref target.fasta --sam output.sam --min_len 1000 --low_score 80.0`
