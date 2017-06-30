# Unicycler-scrub

Unicycler-scrub is a tool which attempts to clean a long read set by trimming read ends and splitting chimeras. It is an experimental tool that is still under development, so use with caution!

[This blog post by Gene Myers](https://dazzlerblog.wordpress.com/2017/04/22/1344/) was a source of inspiration, and like that approach, Unicycler-scrub uses read-read alignments to intrinsically find break points in reads (not any external reference). The end result should be a cleaned set of reads with less erroneous sequence and fewer chimeras, more useful for assembly than the raw reads.



# Trimming



# Splitting



# Usage

```
usage: unicycler_scrub [-h] -i INPUT -o OUT [-r READS] [--trim TRIM] [--split SPLIT] [--min_split_size MIN_SPLIT_SIZE]
                       [--discard_chimeras] [-t THREADS] [--keep_paf] [--parameters PARAMETERS] [--verbosity VERBOSITY]

Unicycler-scrub - read trimming, chimera detection and misassembly detection

optional arguments:
  -h, --help                       show this help message and exit
  -i INPUT, --input INPUT          These are the reads or assembly to be scrubbed (can be FASTA or FASTQ format
  -o OUT, --out OUT                The scrubbed reads or assembly will be saved to this file (will have the same format as the --input
                                   file format) or use "none" to not produce an output file
  -r READS, --reads READS          These are the reads used to scrub --input (can be FASTA or FASTQ format) (default: same file as
                                   --input)
  --trim TRIM                      The aggressiveness with which the input will be trimmed (0 to 100, where 0 is no trimming and 100 is
                                   very aggressive trimming) (default: 50)
  --split SPLIT                    The aggressiveness with which the input will be split (0 to 100, where 0 is no splitting and 100 is
                                   very aggressive splitting) (default: 50)
  --min_split_size MIN_SPLIT_SIZE  Parts of split sequences will only be outputted if they are at least this big (default: 1000)
  --discard_chimeras               If used, chimeric sequences will be discarded instead of split (default: False)
  -t THREADS, --threads THREADS    Number of threads used (default: 8)
  --keep_paf                       Save the alignments to file (makes repeated runs faster because alignments can be loaded from file)
                                   (default: False)
  --parameters PARAMETERS          Low-level parameters (for debugging use only) (default: )
  --verbosity VERBOSITY            Level of stdout information (default: 1)
                                     0 = no stdout, 1 = basic progress indicators, 2 = extra info, 3 = debugging info
```
