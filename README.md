# Genome-assembly-and-pop-gen-stuff

This is a collection of useful perl scripts to perform various tasks for genome assembly and population genetics related projects

i) assembly_validation.pl is a perl script that validates the completeness and accurateness of a newly assembled genome by mapping all reads sequenced in the genome assembly back to the assembly, and compute the proportion of mates that are correctly oriented. 

ii) calculate_nucleotidediversity_windows.pl is a perl script that calculates the heterozygosity from a sequence alignment by chopping it into windows of flexible length, and incorporating nucleotide diversity to the alignment divergence data

iii) parse_blast.pl is a script to parse blast output and generate some useful summaries
 
