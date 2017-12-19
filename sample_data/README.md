# Unicycler sample data

I've put together a few small read sets so users can test that Unicycler works.

The synthetic _Shigella_ plasmid reads are the smallest in size and included in the Unicycler repo – try these if you're in a hurry.

The other three are real read sets from small bacterial genomes from the [FDA-ARGOS project](https://www.ncbi.nlm.nih.gov/bioproject/231221) and are available to download via [figshare](https://figshare.com/projects/Unicycler_sample_data/23065). The _Helicobacter pylori_ and _Streptococcus pyogenes_ genomes are relatively simple and easy to assemble. The _Neisseria gonorrhoeae_ genome is complex and tougher. I subsampled each Illumina read set down to create smaller files. The PacBio read sets were subsampled based on quality (i.e. they are a high-quality subset of the original reads).

I'd recommend looking at the resulting assembly graphs in [Bandage](https://github.com/rrwick/Bandage) to get an idea of how well the assemblies completed – especially useful for comparing hybrid assemblies made with low-depth vs high-depth long reads.


### _Shigella sonnei_ plasmids (synthetic reads)

These are synthetic reads from plasmids A, B and E from the [_Shigella sonnei_ 53G genome assembly](https://www.ncbi.nlm.nih.gov/genome/417?genome_assembly_id=166795):

Download reads [from the figshare page](https://figshare.com/articles/Synthetic_Shigella_plasmid_reads/5165776) or via these direct links:
* [short_reads_1.fastq.gz](https://github.com/rrwick/Unicycler/raw/master/sample_data/short_reads_1.fastq.gz), [short_reads_2.fastq.gz](https://github.com/rrwick/Unicycler/raw/master/sample_data/short_reads_2.fastq.gz)<br>
* [long_reads_low_depth.fastq.gz](https://github.com/rrwick/Unicycler/raw/master/sample_data/long_reads_low_depth.fastq.gz), [long_reads_high_depth.fastq.gz](https://github.com/rrwick/Unicycler/raw/master/sample_data/long_reads_high_depth.fastq.gz)

These plasmids are small compared to a bacterial genome, but insertion sequences create many repeats. Only the smallest plasmid assembles completely with short reads alone. Hybrid assemblies with low-depth long reads manage to complete the medium-sized plasmid, and it takes high-depth long reads to complete all three.


### _Helicobacter pylori_

These are real Illumina and PacBio reads from [_Helicobacter pylori_ sample FDAARGOS_300](https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN06173313):

Download reads [from the figshare page](https://figshare.com/articles/Helicobacter_pylori_SAMN06173313_reads/5165782) or via these direct links:
* [short_reads_1.fastq.gz](https://ndownloader.figshare.com/files/8801860), [short_reads_2.fastq.gz](https://ndownloader.figshare.com/files/8801863)<br>
* [long_reads_low_depth.fastq.gz](https://ndownloader.figshare.com/files/8801857), [long_reads_high_depth.fastq.gz](https://ndownloader.figshare.com/files/8801854)

The _Helicobacter pylori_ genome is small and simple. It has only two copies of the RNA operon and no other large repeats, making it very easy to assemble compared to most bacterial genomes. A hybrid assembly with the high-depth long reads should produce a nice completed chromosome. A hybrid assembly with the low-depth long reads comes very close to completion, with just a couple of slightly ambiguous spots remaining.


### _Streptococcus pyogenes_

These are real Illumina and PacBio reads from [_Streptococcus pyogenes_ sample FDAARGOS_190](https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN04875527):

Download reads [from the figshare page](https://figshare.com/articles/Streptococcus_pyogenes_SAMN04875527_reads/5165788) or via these direct links:
* [short_reads_1.fastq.gz](https://ndownloader.figshare.com/files/8801875), [short_reads_2.fastq.gz](https://ndownloader.figshare.com/files/8801878)<br>
* [long_reads_low_depth.fastq.gz](https://ndownloader.figshare.com/files/8801872), [long_reads_high_depth.fastq.gz](https://ndownloader.figshare.com/files/8801869)

The _Streptococcus pyogenes_ genome is particularly small and simple and is relatively easy to assemble with Illumina reads. It does have a few repetitive elements, however, including five copies of the RNA operon and six copies of IS1548. A hybrid assembly with the high-depth long reads should produce a nice completed chromosome. A hybrid assembly with the low-depth long reads will not quite complete, leaving a bit of ambiguity around some of the RNA operons.


### _Neisseria gonorrhoeae_

These are real Illumina and PacBio reads from [_Neisseria gonorrhoeae_ sample FDAARGOS_204](https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN04875541):

Download reads [from the figshare page](https://figshare.com/articles/Neisseria_gonorrhoeae_SAMN04875541_reads/5165809) or via these direct links:
* [short_reads_1.fastq.gz](https://ndownloader.figshare.com/files/8802211), [short_reads_2.fastq.gz](https://ndownloader.figshare.com/files/8802214)<br>
* [long_reads_low_depth.fastq.gz](https://ndownloader.figshare.com/files/8802208), [long_reads_high_depth.fastq.gz](https://ndownloader.figshare.com/files/8802205)

While the _Neisseria gonorrhoeae_ genome is small, it is a difficult one to assemble, with many copies of IS<i>1016</i>, IS<i>Ngo2</i> and other repeats. A hybrid assembly with the high-depth long reads should produce a nice completed chromosome. A hybrid assembly with the low-depth long reads, while still a large improvement over the Illumina-only assembly, fails to resolve in a number of regions. This demonstrates that more complex genomes require higher long-read-depth to achieve complete assemblies.


### Assembly commands

__Illumina-only assembly:__<br>
`unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -o output_dir`

__Long-read-only assembly:__<br>
`unicycler -l long_reads_high_depth.fastq.gz -o output_dir`

__Hybrid assembly (low-depth long reads):__<br>
`unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -l long_reads_low_depth.fastq.gz -o output_dir`

__Hybrid assembly (high-depth long reads):__<br>
`unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -l long_reads_high_depth.fastq.gz -o output_dir`



