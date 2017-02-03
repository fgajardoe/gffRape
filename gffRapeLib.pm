#===============================================================================
#         FILE: misc.pm
#
#  DESCRIPTION: Misc fuctions for code rehuse
#       AUTHOR: Felipe Gajardo (), fgajardoe@gmail.com
#      VERSION: 1.0
#      CREATED: 04/14/2015 03:51:20 PM
#===============================================================================

#package Misc;
use strict;
use warnings;
use Data::Dumper;
#use Exporter qw(import);
#our @EXPORT_OK=qw(gff2genemodel);

# re-implementacion de gff2genemodel #########################################################
#ESTA COMENTADO PARA QUE NO WEBEE
# Clase que quiero crear para hacerlo mas dinamico
#sub element{
#	my $(start,end,strand,phase,score,attributes,type,method)=shift;
#	my $e->{start}=$start;
#	my $e->{end}=$end;
#	my $e->{strand}=$strand;
#	my $e->{phase}=$phase;
#	my $e->{score}=$score;
#	my $e->{attributes}=$attributes;
#	my $e->{type}=$type;
#	my $e->{method}=$method;
#
#	
#}


# Title		: gff2genemodel
# Usage		: my $g=gff2genemodel("path/to/gene_model.gff")
#
# Function	: Generate a hash reference from a gff3 file.
# Returns	: Hash reference
# Args		: a string containing the path to gff3 file
# Note		: only gene, mRNA, CDS, exon, intron, start_codon and stop_codon features supported.

sub gff2genemodel{

	my $ingff=shift;
	open GFF, "<$ingff" or die "Cannot open $ingff";
	
	my %gene=();
	my %mrna=();
	my %exon=();
	my %intron=();
	my %cds=();
	my %start_codon=();
	my %stop_codon=();
#	my %transposon


	my %count_type=(gene=>0,mRNA=>0,transcript=>0,exon=>0,intron=>0,CDS=>0,start_codon=>0,stop_codon=>0);

	my $dup=0; # contador necesario para renombrar las entradas con ids duplicados

	# recorremos el gff
	my @lines=<GFF>;
	foreach my $line(@lines){
		if(($line=~/^#/g)|(length($line) == 1)){
			next;
		}

		chomp $line;

		# obtenemos las columnas
		my @col=split /\t/,$line;

		# saltamos lineas con atributos no soportados (enero 2017)
		if(!exists $count_type{$col[2]}){
			next;
		}

		# obtenemos sus atributos. (ademas creamos un hash con ellos)
		my @attributes=split /;/,$col[8];
		my %att;
		foreach my $attribute(@attributes){
			my @kv=split /=/,$attribute;
			$att{$kv[0]}=$kv[1];
		}


		# hay algunas entradas sin ID. para estas se lo inventamos
		if(!$att{ID}){
			my $id=$col[2].$count_type{$col[2]};
			$att{ID}=$id;	
		}

		if($col[2] eq 'gene'){

			$gene{$att{ID}}->{phase}=$col[7];
			$gene{$att{ID}}->{type}=$col[2];
			$gene{$att{ID}}->{seqid}=$col[0];
			$gene{$att{ID}}->{source}=$col[1];
			$gene{$att{ID}}->{start}=$col[3];
			$gene{$att{ID}}->{end}=$col[4];
			$gene{$att{ID}}->{strand}=$col[6];
			$gene{$att{ID}}->{score}=$col[5];
			$gene{$att{ID}}->{attributes}=\%att;
			$count_type{$col[2]}++;
		}
		if(($col[2] eq 'mRNA')||($col[2] eq 'transcript')){
			$mrna{$att{ID}}->{type}=$col[2];
			$mrna{$att{ID}}->{seqid}=$col[0];
			$mrna{$att{ID}}->{source}=$col[1];
			$mrna{$att{ID}}->{start}=$col[3];
			$mrna{$att{ID}}->{end}=$col[4];
			$mrna{$att{ID}}->{strand}=$col[6];
			$mrna{$att{ID}}->{score}=$col[5];
			$mrna{$att{ID}}->{phase}=$col[7];
			#my %feats=$mrna{$att{ID}}

			$mrna{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
 
			#$gene{$att{Parent}}->{mrna}->{$att{ID}}=\%mrna;

			$count_type{$col[2]}++;

		}
		if($col[2] eq 'exon'){
			

			# obtenemos el gene parent.
			#my $parent=$att{'Parent'};
			#my $gene_parent=$mrna{$parent}->{attributes}->{Parent};
			$exon{$att{ID}}->{type}=$col[2];
			$exon{$att{ID}}->{seqid}=$col[0];
			$exon{$att{ID}}->{source}=$col[1];
			$exon{$att{ID}}->{start}=$col[3];
			$exon{$att{ID}}->{end}=$col[4];
			$exon{$att{ID}}->{strand}=$col[6];
			$exon{$att{ID}}->{score}=$col[5];
			$exon{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
			$exon{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%exon;

 			$count_type{$col[2]}=int($count_type{$col[2]})+1;

		}
		if($col[2] eq 'CDS'){
			# hay ocasiones donde hay ids repetidos.(me paso con el gff output de augustus). para los CDS lo admitimos ya que un cds no va a ser parent de nadie.
			if(exists $cds{$att{ID}}){
				print STDERR "Warning: ".$col[2]." ID(".$att{ID}.") duplicated. Included anyway as ".$att{ID}."_dup".$dup."\n";

				$att{ID}=$att{ID}."_dup".$dup;
				$dup++;
			}
			
			# obtenemos el gene parent.
			#my $parent=$att{'Parent'};
			#my $gene_parent=$mrna{$parent}->{attributes}->{Parent};
			$cds{$att{ID}}->{phase}=$col[7];
			$cds{$att{ID}}->{type}=$col[2];
			$cds{$att{ID}}->{seqid}=$col[0];
			$cds{$att{ID}}->{source}=$col[1];
			$cds{$att{ID}}->{start}=$col[3];
			$cds{$att{ID}}->{end}=$col[4];
			$cds{$att{ID}}->{strand}=$col[6];
			$cds{$att{ID}}->{score}=$col[5];
			$cds{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.

 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%exon;

 			$count_type{$col[2]}=int($count_type{$col[2]})+1;

		}
		if($col[2] eq 'intron'){

			# obtenemos el gene parent.
			#my $gene_parent=$mrna{$att{Parent}}->{attributes}->{Parent};
			$intron{$att{ID}}->{type}=$col[2];
			$intron{$att{ID}}->{seqid}=$col[0];
			$intron{$att{ID}}->{source}=$col[1];
			$intron{$att{ID}}->{start}=$col[3];
			$intron{$att{ID}}->{end}=$col[4];
			$intron{$att{ID}}->{strand}=$col[6];
			$intron{$att{ID}}->{score}=$col[5];
			$intron{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
			$intron{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%intron;
			$count_type{$col[2]}++;
		}

		if($col[2] eq 'start_codon'){

			# obtenemos el gene parent.
			#my $gene_parent=$mrna{$att{Parent}}->{attributes}->{Parent};
			$start_codon{$att{ID}}->{type}=$col[2];
			$start_codon{$att{ID}}->{seqid}=$col[0];
			$start_codon{$att{ID}}->{source}=$col[1];
			$start_codon{$att{ID}}->{start}=$col[3];
			$start_codon{$att{ID}}->{end}=$col[4];
			$start_codon{$att{ID}}->{strand}=$col[6];
			$start_codon{$att{ID}}->{score}=$col[5];
			$start_codon{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
			$start_codon{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%start_codon;
			$count_type{$col[2]}++;
		}

		if($col[2] eq 'stop_codon'){

			# obtenemos el gene parent.
			#my $gene_parent=$mrna{$att{Parent}}->{attributes}->{Parent};
			$stop_codon{$att{ID}}->{type}=$col[2];
			$stop_codon{$att{ID}}->{seqid}=$col[0];
			$stop_codon{$att{ID}}->{source}=$col[1];
			$stop_codon{$att{ID}}->{start}=$col[3];
			$stop_codon{$att{ID}}->{end}=$col[4];
			$stop_codon{$att{ID}}->{strand}=$col[6];
			$stop_codon{$att{ID}}->{score}=$col[5];
			$stop_codon{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
			$stop_codon{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%stop_codon;
			$count_type{$col[2]}++;
		}



	}# foreach line

	# ya recolectamos la informacion, ahora debemos unirla.
	# primero los exones
	foreach my $id(keys %exon){
		my $parent=$exon{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{exon}->{$id}=$exon{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}
	# tambien los cds
	foreach my $id(keys %cds){
		my $parent=$cds{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{CDS}->{$id}=$cds{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}


	# luego los intrones
	foreach my $id(keys %intron){
		my $parent=$intron{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			#$mrna{$parent}->{intron}->{$id}=\$intron{$id}; # aqui cambia \ por %
			$mrna{$parent}->{intron}->{$id}=$intron{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}
	
	# los codones start
	foreach my $id(keys %start_codon){
		my $parent=$start_codon{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{start_codon}->{$id}=$start_codon{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}

	# los codones stop
	foreach my $id(keys %stop_codon){
		my $parent=$stop_codon{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{stop_codon}->{$id}=$stop_codon{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}

	# finalmente los mrna a los genes
	foreach my $id(keys %mrna){
		my $parent=$mrna{$id}->{attributes}->{Parent};
		if(exists $gene{$parent}){
			#$gene{$parent}->{mrna}->{$id}=\$mrna{$id};
			$gene{$parent}->{mRNA}->{$id}=$mrna{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}

	
	
	close GFF;
	return \%gene;
}


# esta es una copia de gff2genemodel, pero modificada para renombrar todos los IDs


sub gff2genemodel_reset_id{

	my $ingff=shift;
	open GFF, "<$ingff" or die "Cannot open $ingff";
	
	my %gene=();
	my %mrna=();
	my %exon=();
	my %intron=();
	my %cds=();
	my %start_codon=();
	my %stop_codon=();
#	my %transposon

	# necesario para guardar la relacion entre los ids nuevos y los viejos
	my %eq_id;

	my %count_type=(gene=>0,mRNA=>0,transcript=>0,exon=>0,intron=>0,CDS=>0,start_codon=>0,stop_codon=>0);

	my $dup=0; # contador necesario para renombrar las entradas con ids duplicados

	# recorremos el gff
	my @lines=<GFF>;
	foreach my $line(@lines){
		if(($line=~/^#/g)|(length($line) == 1)){
			next;
		}

		chomp $line;

		# obtenemos las columnas
		my @col=split /\t/,$line;

		# saltamos lineas con atributos no soportados (enero 2017)
		if(!exists $count_type{$col[2]}){
			next;
		}

		# obtenemos sus atributos. (ademas creamos un hash con ellos)
		my @attributes=split /;/,$col[8];
		my %att;
		foreach my $attribute(@attributes){
			my @kv=split /=/,$attribute;
			$att{$kv[0]}=$kv[1];
		}


		# hay algunas entradas sin ID. para TODAS se lo inventamos
#		if(!$att{ID}){
			my $original_id;
			my $original_parent;
			my $id=$col[2].$count_type{$col[2]};
			if($att{ID}){
				$original_id=$att{ID};
			}
			else{
				$original_id=$id;
			}
			$att{ID}=$id;
			$eq_id{$original_id}=$id;
			print STDERR "$id\t$original_id\n";
			#y guardamos/reasignamos la relacion de padre
			if($col[2] ne 'gene'){
				$original_parent=$att{Parent};
				$att{Parent}=$eq_id{$original_parent};
			}
#		}

		if($col[2] eq 'gene'){

			$gene{$att{ID}}->{phase}=$col[7];
			$gene{$att{ID}}->{type}=$col[2];
			$gene{$att{ID}}->{seqid}=$col[0];
			$gene{$att{ID}}->{source}=$col[1];
			$gene{$att{ID}}->{start}=$col[3];
			$gene{$att{ID}}->{end}=$col[4];
			$gene{$att{ID}}->{strand}=$col[6];
			$gene{$att{ID}}->{score}=$col[5];
			$gene{$att{ID}}->{attributes}=\%att;
			$count_type{$col[2]}++;
		}
		if(($col[2] eq 'mRNA')||($col[2] eq 'transcript')){
			$mrna{$att{ID}}->{type}=$col[2];
			$mrna{$att{ID}}->{seqid}=$col[0];
			$mrna{$att{ID}}->{source}=$col[1];
			$mrna{$att{ID}}->{start}=$col[3];
			$mrna{$att{ID}}->{end}=$col[4];
			$mrna{$att{ID}}->{strand}=$col[6];
			$mrna{$att{ID}}->{score}=$col[5];
			$mrna{$att{ID}}->{phase}=$col[7];
			#my %feats=$mrna{$att{ID}}

			$mrna{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
 
			#$gene{$att{Parent}}->{mrna}->{$att{ID}}=\%mrna;

			$count_type{$col[2]}++;

		}
		if($col[2] eq 'exon'){
			

			# obtenemos el gene parent.
			#my $parent=$att{'Parent'};
			#my $gene_parent=$mrna{$parent}->{attributes}->{Parent};
			$exon{$att{ID}}->{type}=$col[2];
			$exon{$att{ID}}->{seqid}=$col[0];
			$exon{$att{ID}}->{source}=$col[1];
			$exon{$att{ID}}->{start}=$col[3];
			$exon{$att{ID}}->{end}=$col[4];
			$exon{$att{ID}}->{strand}=$col[6];
			$exon{$att{ID}}->{score}=$col[5];
			$exon{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
			$exon{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%exon;

 			$count_type{$col[2]}=int($count_type{$col[2]})+1;

		}
		if($col[2] eq 'CDS'){
			# hay ocasiones donde hay ids repetidos.(me paso con el gff output de augustus). para los CDS lo admitimos ya que un cds no va a ser parent de nadie.
			if(exists $cds{$att{ID}}){
				print STDERR "Warning: ".$col[2]." ID(".$att{ID}.") duplicated. Included anyway as ".$att{ID}."_dup".$dup."\n";

				$att{ID}=$att{ID}."_dup".$dup;
				$dup++;
			}
			
			# obtenemos el gene parent.
			#my $parent=$att{'Parent'};
			#my $gene_parent=$mrna{$parent}->{attributes}->{Parent};
			$cds{$att{ID}}->{phase}=$col[7];
			$cds{$att{ID}}->{type}=$col[2];
			$cds{$att{ID}}->{seqid}=$col[0];
			$cds{$att{ID}}->{source}=$col[1];
			$cds{$att{ID}}->{start}=$col[3];
			$cds{$att{ID}}->{end}=$col[4];
			$cds{$att{ID}}->{strand}=$col[6];
			$cds{$att{ID}}->{score}=$col[5];
			$cds{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.

 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%exon;

 			$count_type{$col[2]}=int($count_type{$col[2]})+1;

		}
		if($col[2] eq 'intron'){

			# obtenemos el gene parent.
			#my $gene_parent=$mrna{$att{Parent}}->{attributes}->{Parent};
			$intron{$att{ID}}->{type}=$col[2];
			$intron{$att{ID}}->{seqid}=$col[0];
			$intron{$att{ID}}->{source}=$col[1];
			$intron{$att{ID}}->{start}=$col[3];
			$intron{$att{ID}}->{end}=$col[4];
			$intron{$att{ID}}->{strand}=$col[6];
			$intron{$att{ID}}->{score}=$col[5];
			$intron{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
			$intron{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%intron;
			$count_type{$col[2]}++;
		}

		if($col[2] eq 'start_codon'){

			# obtenemos el gene parent.
			#my $gene_parent=$mrna{$att{Parent}}->{attributes}->{Parent};
			$start_codon{$att{ID}}->{type}=$col[2];
			$start_codon{$att{ID}}->{seqid}=$col[0];
			$start_codon{$att{ID}}->{source}=$col[1];
			$start_codon{$att{ID}}->{start}=$col[3];
			$start_codon{$att{ID}}->{end}=$col[4];
			$start_codon{$att{ID}}->{strand}=$col[6];
			$start_codon{$att{ID}}->{score}=$col[5];
			$start_codon{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
			$start_codon{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%start_codon;
			$count_type{$col[2]}++;
		}

		if($col[2] eq 'stop_codon'){

			# obtenemos el gene parent.
			#my $gene_parent=$mrna{$att{Parent}}->{attributes}->{Parent};
			$stop_codon{$att{ID}}->{type}=$col[2];
			$stop_codon{$att{ID}}->{seqid}=$col[0];
			$stop_codon{$att{ID}}->{source}=$col[1];
			$stop_codon{$att{ID}}->{start}=$col[3];
			$stop_codon{$att{ID}}->{end}=$col[4];
			$stop_codon{$att{ID}}->{strand}=$col[6];
			$stop_codon{$att{ID}}->{score}=$col[5];
			$stop_codon{$att{ID}}->{attributes}=\%att; # La idea de esta linea es almacenar los atributos.
			$stop_codon{$att{ID}}->{phase}=$col[7];
 			#$gene{$gene_parent}->{mrna}->{$att{Parent}}->{exon}->{$att{ID}}=\%stop_codon;
			$count_type{$col[2]}++;
		}



	}# foreach line

	# ya recolectamos la informacion, ahora debemos unirla.
	# primero los exones
	foreach my $id(keys %exon){
		my $parent=$exon{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{exon}->{$id}=$exon{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}
	# tambien los cds
	foreach my $id(keys %cds){
		my $parent=$cds{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{CDS}->{$id}=$cds{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}


	# luego los intrones
	foreach my $id(keys %intron){
		my $parent=$intron{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			#$mrna{$parent}->{intron}->{$id}=\$intron{$id}; # aqui cambia \ por %
			$mrna{$parent}->{intron}->{$id}=$intron{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}
	
	# los codones start
	foreach my $id(keys %start_codon){
		my $parent=$start_codon{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{start_codon}->{$id}=$start_codon{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}

	# los codones stop
	foreach my $id(keys %stop_codon){
		my $parent=$stop_codon{$id}->{attributes}->{Parent};
		if(exists $mrna{$parent}){
			$mrna{$parent}->{stop_codon}->{$id}=$stop_codon{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}

	# finalmente los mrna a los genes
	foreach my $id(keys %mrna){
		my $parent=$mrna{$id}->{attributes}->{Parent};
		if(exists $gene{$parent}){
			#$gene{$parent}->{mrna}->{$id}=\$mrna{$id};
			$gene{$parent}->{mRNA}->{$id}=$mrna{$id};
		}
		else{
			print STDERR "Warning: Parent($parent) for $id not found. Skipped.";
		}
	}

	
	
	close GFF;
	return \%gene;
}


# print gene and all its childs. this function is RECURSIVE
sub printGene{
	my $gene=shift;
	my @valid_childs=("mRNA","exon","intron","CDS","start_codon","stop_codon","transcript");
	
	# Imprimimos el elemento actual
		
	my $att=$gene->{attributes};
	my $attributes="";
	foreach my $k(keys %$att){
		$attributes=$attributes.";".$k."=".$$att{$k};
	}
	$attributes=~s/^;//;

	my $gene_line=join("\t",$gene->{seqid},$gene->{source},$gene->{type},$gene->{start},$gene->{end},$gene->{score},$gene->{strand},$gene->{phase},$attributes);
	print $gene_line."\n";

	# Imprimimos los elementos dentro de el (mrna,cds,exon...)
	
	foreach my $k(keys %$gene){
		foreach my $valid(@valid_childs){
			if($k eq $valid){
				my $childs=$gene->{$k};
				foreach my $child_id(keys %$childs){
					my $child=$childs->{$child_id};
					printGene($child);
				}

			}
		}
	}

}


# filter transcripts on a gene structure. if returns 0, no transcript associated to the gene
# options are, "longest", that keep only the longest transcript; and "shortest", that keep only shortest transcript. if transcripts are equally length, keep just one. the last that it finds
# usage: $filtered_gene_structure=filterGeneTranscripts($gene_structure, 'longest')
sub filterGeneTranscripts{
	my $gene=shift;
	my $opt=shift; # puede ser longest o shortest
	
	if(!exists $gene->{mRNA}){
		print STDERR "Gene not contain mRNA childs. Skipped\n";
		return 0;
	}

	my $mrnas=$gene->{mRNA};

	my $last_len=0;
	my $valid_id;
	foreach my $id(keys %$mrnas){
		my $len=$$mrnas{$id}->{end} - $$mrnas{$id}->{start};
		if($opt eq 'longest'){
			if($len>=$last_len){
				$valid_id=$id;
			}
			else{
				print STDERR "Warning: $id length ($len) < 0?\n";
			}
		}
		if($opt eq 'shortest'){
			if($last_len<=0){
				$last_len=$len;
				$valid_id=$id;
			}
			else{
				if($len<$last_len){
					$last_len=$len;
					$valid_id=$id;
				}
			}
		}
		if(($opt ne 'longest')&&($opt ne 'shortest')){
			die "$opt is not valid option for filterGeneTranscripts";
		}
	}
	# creamos la estructura del gen ya filtrada
	my $filtered_gene=$gene;
	my $valid_mrna=$gene->{mRNA}->{$valid_id};
	$gene->{mRNA}=();
	$gene->{mRNA}->{$valid_id}=$valid_mrna;

	return $gene;
}

return 1;
