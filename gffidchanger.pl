#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: gffextract.pl
#
#        USAGE: ./gffextract.pl  
#
#  DESCRIPTION: 
#       AUTHOR: Felipe Gajardo (), fgajardoe@gmail.com
#      CREATED: 04/14/2015 05:35:05 PM
#
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Std;
use Data::Dumper;
#use Misc; #qw(gff2genemodel printGene filterGeneTranscripts);
use gffRapeLib;

our($opt_h,$opt_i,$opt_p,$opt_z,$opt_Z);


# Usage
my $usage="Usage: perl $0 -i input.gff
Options:
			-i input_gff	Provide input GFF to reset IDs (mandatory).
			-p prefix		Prefix to use on new IDs (mandatory)
			-z gene=int,mRNA=int...	Limit for zero filled counter to 'int'. for each feature. if no declared default is no 0-filling for that feature.
			-Z gene=str,mRNA=str...	Substitute feature name for 'str' on new IDs
			-h				Print this help

";

# Get options
getopts('hp:i:z:Z');

# Check if there is input
if(($opt_h)||(!$opt_i)||(!$opt_p)){
	print $usage;
	exit 0;
}

# Main
my $genes;


my $zh; #hash con limite para el largo de numeros en el identificador. (x lo del relleno con 0s)
#my $Zh;
#if($opt_Z){
#	my @c=split /,/,$opt_Z;
#	foreach my $feat(@c){
#		my @f=split /=/,$feat;
#		$Zh->{$f[0]}=$f[1];
#	}
#}

if($opt_z){

	my @c=split /,/,$opt_z;
	foreach my $feat(@c){
		my @f=split /=/,$feat;
		$zh->{$f[0]}=$f[1];
	}

	$genes=gff2genemodel_reset_id($opt_i,$opt_p,$opt_z,$zh); #,$Zh);
}
else{
	$genes=gff2genemodel_reset_id($opt_i,$opt_p);
}
my @genelst;
#if($opt_l){
#	open LST,"<$opt_l" or die "Cannot open $opt_l";
#	@genelst=<LST>;
#}
#else{
	@genelst=keys %$genes;
#}

#recorreemos la lista de genes
foreach my $idlst(@genelst){
	chomp $idlst;
	my $gene=$genes->{$idlst};
	my $filtered_gene=$gene;

#	if($opt_k ne 'all'){
#		$filtered_gene=filterGeneTranscripts($gene,$opt_k);
#	}
	if($filtered_gene!=0){
		printGene($filtered_gene);
	}
}


