# Unicycler-check

Unicycler-check uses long-read alignments to assess the contiguity of an assembly. It aims to highlight misassemblies or otherwise suspicious parts of the assembly.

__Important:__ This tool is experimental and no longer under development, so use with caution! It may disappear from future releases of Unicycler without warning.



### Requirements

Unicycler-check uses [Plotly](https://plot.ly/) to generate its html output, so that must be installed for Python 3.



### Usage

```
usage: unicycler_check [-h] --sam SAM --ref REF --reads READS [--min_len MIN_LEN] [--error_window_size ERROR_WINDOW_SIZE]
                       [--depth_window_size DEPTH_WINDOW_SIZE] [--error_rate_threshold ERROR_RATE_THRESHOLD] [--depth_p_val DEPTH_P_VAL]
                       [--window_tables WINDOW_TABLES] [--base_tables BASE_TABLES] [--html HTML] [--threads THREADS]
                       [--verbosity VERBOSITY]

Long read assembly checker

optional arguments:
  -h, --help                            show this help message and exit
  --sam SAM                             Input SAM file of alignments (if this file doesn't exist, the alignment will be performed with
                                        results saved to this file - you can use the aligner arguments with this script)
  --ref REF                             FASTA file containing one or more reference sequences
  --reads READS                         FASTQ file of long reads
  --min_len MIN_LEN                     Minimum alignment length (bp) - exclude alignments shorter than this length (default: 100)
  --error_window_size ERROR_WINDOW_SIZE
                                        Window size for error summaries (default: 100)
  --depth_window_size DEPTH_WINDOW_SIZE
                                        Window size for depth summaries (default: 100)
  --error_rate_threshold ERROR_RATE_THRESHOLD
                                        Threshold for high error rates, expressed as the fraction between the mean error rate and the
                                        random alignment error rate (default: 0.3)
  --depth_p_val DEPTH_P_VAL             P-value for low/high depth thresholds (default: 0.001)
  --window_tables WINDOW_TABLES         Path and/or prefix for table files summarising reference errors for reference windows (default: do
                                        not save window tables)
  --base_tables BASE_TABLES             Path and/or prefix for table files summarising reference errors at each base (default: do not save
                                        base tables)
  --html HTML                           Path for HTML report (default: do not save HTML report)
  --threads THREADS                     Number of CPU threads used to align (default: the number of available CPUs)
  --verbosity VERBOSITY                 Level of stdout information (0 to 2) (default: 1)
```
