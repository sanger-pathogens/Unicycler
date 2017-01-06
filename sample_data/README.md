# Unicycler sample data

This sample dataset is made of synthetic reads generated from three plasmids in a [_Shigella sonnei reference_](https://www.ncbi.nlm.nih.gov/genome/417?genome_assembly_id=166795).


### Short read-only assembly

`unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -o output_dir`

This command will assemble the short reads alone. If you look at the resulting contigs (or view the graph in [Bandage](https://github.com/rrwick/Bandage)) you'll see that only the smallest plasmid assembled completely. The larger two plasmids contain quite a lot of repetitive sequence and short reads aren't enough to complete them.


### Low-depth long read hybrid assembly

`unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -l long_reads_low_depth.fastq.gz -o output_dir`

This command will assemble the short reads along with a small number of long reads. They have an average depth of ~1x, so some parts of the plasmids will be represented in these reads but other parts will not.

These reads are not sufficient to complete the whole assembly, but we are getting closer. The medium-sized plasmid should now be finished and only the largest plasmid remains incomplete.


### High-depth long read hybrid assembly

`unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -l long_reads_high_depth.fastq.gz -o output_dir`

This command will assemble the short reads along using long reads of ~20x depth. Now all parts of the plasmids are represented in long reads so there is sufficient information to complete the assemblies. Accordingly, you should now see that the assembly produces just three circular contigs: one for each plasmid.
