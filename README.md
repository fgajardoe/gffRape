# gffRape

## Introduction

This is a collection of scripts to filter and extract relevant biological information from
GFF3 files ( https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md ).
All here still is very disorganized. Also I am learning and getting used to Git. So patience :)

Currently the next GFF features are supporteed. Any other will be skipped:

* gene
* mRNA
* transcript
* exon
* intron
* CDS
* start_codon
* stop_codon

## Dependences

This perl modules:
* Getopt::Std
* Data::Dumper


## Files

Currently there is only one file (gffRapeLib.pm) which have the backbone of all the suite.
The other perl scripts are aplications of the functions developed in gffRapeLib.pm

Here you can see a brief description of each one:

* **gffidchanger.pl**
Reset all identifiers in the gff file. New IDs are the feature name and an incremental number (example: gene1, gene2, gene3, ...).
Equivalences with old IDs are printed to STDERR; GFF with new IDs is printed to STDOUT.

* **gffextract.pl**
Extract genes from a GFF files. Outputs a new GFF with genes (and associated entries) of your interest (from a list of genes).
Also is capable of filtering isoforms, keeping the longest one, the shortest one, or all (default)

* **gffGeneStats.pl**
Calc stats from genes on gff. Outputs a table with this columns:
geneID, gene name, geneLength, number of exons, exon mean length, exon total length, number of introns, intron mean length, intron total length
At the moment, this script does not rely on gffRapeLib.pm functions (and support less features, only: gene, mRNA, exon, intron).

## Author

Felipe Gajardo, Bioinformatician at Center for Genome Regulation, Chile
