# UCLA CS CM124

Code for the Spring 2015 offering of UCLA's Computer Science CM124: Computational Genetics, taught by Professor Eskin.

[Final project blog located here.](http://www.ryanwb.com/cs-cm124)

## Final Project Description
### Haplotype Phasing (Minimum Parsimony Formulation)

A haplotype is a set of SNP "values" on a single chromosome. Sequencing technology, however, provides genotypes (not haplotypes), which do not provide information about which chromosome (maternal vs. paternal) specific SNPs are from; in other words, these "unresolved haplotypes" do not tell us anything about which SNPs were grouped with each other on a chromosome. Haplotype phasing is the process of resolving this information, grouping the SNPs according to their chromosome of origin. Studies in the early 2000s found that SNPs on each haplotype are correlated and that haplotypes have "limited diversity"; that is, local regions tend to have few unique haplotypes. Therefore, one technique for resolving haplotypes is to choose the configuration which minimizes the number of distinct haplotypes which are present ("minimum parsimony"). This project implements and analyzes algorithms for optimally (with minimum parsimony) phasing haplotypes.

Author:
* [Ryan Baker](http://github.com/ryanwb)
