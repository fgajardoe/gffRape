#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: gffstructure.pl
#
#        USAGE: ./gffstructure.pl  
#
#  DESCRIPTION: Este script es bacan. lee el gff usando reference hashes. lo que
#  				nos permite generar una estructura y extraer la informacion que
#  				nos interesa. hasta el momento soporta solo tres tipos de tipos
#  				de feature:
#
#  					- gene
#  					- mRNA
#  					- exon
#  					- intron
#
#  				Lo importante aqui es la subrutina get_gene_model. En este caso
#  				particular se utilizo para obtener los largos promedio de exones
#  				e intrones, ademas de otra informacion relacionada.
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Felipe Gajardo (), fgajardoe@gmail.com
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 11/19/2014 04:35:42 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Data::Dumper;
use gffRapeLib;

my $usage="Usage: perl $0 input.gff\n";


if(scalar @ARGV<1){
	print $usage;
	exit 0;
}

my $infile=$ARGV[0];
my $genemodel=gff2genemodel($infile);
#my $genemodel=get_gene_model($infile);

#print Dumper $genemodel;
#exit 0;

# revisamos cada gen para obtener los largos
print "geneID\tName\tgeneLength\tNexons\texonMeanLength\texonTotalLength\tNintrons\tintronMeanLength\tintronTotalLength\n";
foreach my $geneid(keys %$genemodel){

	# calculamos el largo del gen
	my $genelen=$$genemodel{$geneid}->{end} - $$genemodel{$geneid}->{start};
	my $exon_mean="NULL";
	my $exon_sum="NULL";
	my $intron_mean="NULL";
	my $intron_sum="NULL";
	my $exon_n=0;
	my $intron_n=0;

	my $acc="NULL";
	if($$genemodel{$geneid}->{attributes}->{Name}){
		$acc=$$genemodel{$geneid}->{attributes}->{Name};
	}


	# comprobamos si el gen tiene algun mRNA asociado.
	if(exists $$genemodel{$geneid}->{mRNA}){
		
		# hacemos otro hash, que correspondera a los mRNA
		my $mrna=$$genemodel{$geneid}->{mRNA};
		

		# obtenemos la isoforma mas larga
		my $genemodel_longest_iso=filterGeneTranscripts($$genemodel{$geneid}, 'longest');
		my @longest_iso_id=keys %$genemodel_longest_iso;

#		print Dumper $genemodel_longest_iso;
#		exit 1;

		# recorremos cada uno buscando el mRNA mas largo. DEPRECATED. ahora lo hacemos con la funcion (ver arriba)
#		my $longest_len=0;
#		my $longest_id;
#		foreach my $mrnaid(keys %$mrna){
#
#			# y calculamos su largo
#			my $start=$$mrna{$mrnaid}->{start};
#			my $end=$$mrna{$mrnaid}->{end};
#			my $mrnalen= $end - $start;
#			
#			if($mrnalen>=$longest_len){
#				$longest_id=$mrnaid;
#				$longest_len=$mrnalen;
#			}
#			else{
#				print STDERR "Warning: $mrnaid length ($mrnalen) < 0 ?. skipped.";
#				next;
#			}
#		}
		# ahora recorremos los exones del mRNA mas largo.
	#	if(exists $$mrna{$longest_id}->{exon}){
		if(exists $genemodel_longest_iso->{mRNA}->{$longest_iso_id[0]}->{exon}){
			$exon_sum=0;
			#my $exon=$$mrna{$longest_id}->{exon};
			my $exon=$$genemodel_longest_iso->{mRNA}->{$longest_iso_id[0]}->{exon};
			foreach my $exonid(keys %$exon){
				my $start=$$exon{$exonid}->{start};
				my $end=$$exon{$exonid}->{end};
				my $exonlen=$end-$start;
				$exon_sum=$exon_sum+$exonlen;
				$exon_n++;
			}
			$exon_mean=$exon_sum/$exon_n;
		}

		# lo mismo para los intrones
		if(exists $genemodel_longest_iso->{mRNA}->{$longest_iso_id[0]}->{intron}){
			$intron_sum=0;

			my $intron=$$genemodel_longest_iso->{mRNA}->{$longest_iso_id[0]}->{intron};
			foreach my $intronid(keys %$intron){
				my $start=$$intron{$intronid}->{start};
				my $end=$$intron{$intronid}->{end};
				my $intronlen=$end-$start;
				$intron_sum=$intron_sum+$intronlen;
				$intron_n++;
			}
			$intron_mean=$intron_sum/$intron_n;
		}
	}
	print "$geneid\t$acc\t$genelen\t$exon_n\t$exon_mean\t$exon_sum\t$intron_n\t$intron_mean\t$intron_sum\n";
}


