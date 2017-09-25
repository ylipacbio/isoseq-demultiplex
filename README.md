The `isoseq-demultiplex` tool tries to demultiplex an barcoded isoseq SMRTCell with two birds samples
based on consensus clustering. It is trying to assign reads to bird samples based on which cluster a
read is in and which cluster do other reads in the same cluster belong to.

Interesting discoveries: FLNC reads usually are well clustered, while many NFL reads can assign to both
birds. This makes biological sense because birds should share many transcripts in common while NFL reads
may miss critical info to be assigned unambiguously.
