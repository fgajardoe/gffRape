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

our($opt_h,$opt_i);


# Usage
my $usage="Usage: perl $0 -i input.gff
Options:
			-i input_gff	Provide input GFF to reset IDs (mandatory).
			-h				Print this help

";

# Get options
getopts('hi:');

# Check if there is input
if(($opt_h)||(!$opt_i)){
	print $usage;
	exit 0;
}

# Main

my $genes=gff2genemodel_reset_id($opt_i);

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


