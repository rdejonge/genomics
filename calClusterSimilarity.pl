#!/usr/bin/perl

####
#Perl script by Ronnie de Jonge, ronnie.dejonge@gmail.com OR ronnie.dejonge@psb.ugent.be
#Created: January 5, 2017
#Last modified: January 5, 2017
#usage: perl calClusterSimilarity.pl [gff3_1] [gff3_2] [mclOutput.iadhore] [out_prefix] [cluster_length] [allowed_gap]
#
#Description: 
#Store gene models (important is protein_id matching mcl output) and order on scaffolds and return for each cluster (of size x) the "amount" of synteny
####

use warnings;
use strict;
use experimental 'smartmatch';
use Data::Dumper;
use List::Compare;
use List::Util qw( min max );
use List::Util qw(sum);

###script usage:
my $usage = "\n\nusage:	$0 [gff3_1] [gff3_2] [mclOutput.iadhore] [out_prefix] [cluster_length (odd)] [allowed_gap] [term1] [term2]\n\n";
my $usage2 = "\n\nusage: [cluster_length] must be odd!!!\n\n";


###Input arguments
my $gff3_1 = $ARGV[0] or die $usage;
my $gff3_2 = $ARGV[1] or die $usage;
my $mcl = $ARGV[2] or die $usage;
my $prefix = $ARGV[3] or die $usage;
my $clus_len = $ARGV[4] or die $usage;
my $allow_gap = $ARGV[5] or die $usage;
my $term1 = $ARGV[6] or die $usage;
my $term2 = $ARGV[7] or die $usage;

if ($clus_len-1 & 1) {
	die $usage2;
}

#define some variables
#my (%genes1, %scaffoldToGene1, %rnaToGene1, %proteinToGene1);
#my (%genes2, %scaffoldToGene2, %rnaToGene2, %proteinToGene2);

##parse gff3 files
my ($genes1, $scaffoldToGene1, $rnaToGene1, $proteinToGene1) = &parse_gff3($gff3_1);
my ($genes2, $scaffoldToGene2, $rnaToGene2, $proteinToGene2) = &parse_gff3($gff3_2);

##dereferencing hashes and save as new hash
my %genes1b = %{$genes1};
my %genes2b = %{$genes2};

my %scaffoldToGene1b = %{$scaffoldToGene1};
my %scaffoldToGene2b = %{$scaffoldToGene2};

my @genes1;
my @genes2;

#parse $mcl and store in %genes1 %genes2
my %mclToGene;
open (MCL, $mcl) or die "can't open this file $mcl";
while (my $line = <MCL>)
{
	chomp $line;
	my ($protein_id, $mcl_group) = split ("\t", $line);
	if ($term1 ne "special") {
		$protein_id =~ s/$term1//;
	}
	if ($term2 ne "special") {
		$protein_id =~ s/$term2//;
	}	
	if (defined $proteinToGene1->{$protein_id})
	{
		my $gid = $proteinToGene1->{$protein_id};
		#print $gid,"\t",$mcl_group,"\n";
		my $scaffold = $scaffoldToGene1->{$gid};
		$genes1->{$scaffold}->{$gid}->{"mcl"} = $mcl_group;
		push (@{$mclToGene{$mcl_group}->{"member"}}, $gid);
	}
	elsif (defined $proteinToGene2->{$protein_id})
	{
		my $gid = $proteinToGene2->{$protein_id};
		#print $gid,"\t",$mcl_group,"\n";
		my $scaffold = $scaffoldToGene2->{$gid};
		$genes2->{$scaffold}->{$gid}->{"mcl"} = $mcl_group;
		push (@{$mclToGene{$mcl_group}->{"member"}}, $gid);
	}	
}

#print Dumper ($genes1);
#print Dumper (\%mclToGene);

##for each genes in %genes check whether flanks exist (minimal #flanks = $clus_len/2)
my $max_clus_len = ($clus_len-1)/2;

##get list of genes from reference and query
my @scaffolds1 = keys %genes1b;
foreach my $scaffold (@scaffolds1) { 
	my @genesOnScaffold = keys %{$genes1b{$scaffold}};
	push (@genes1, @genesOnScaffold);
}
my @scaffolds2 = keys %genes2b;
foreach my $scaffold (@scaffolds2) { 
	my @genesOnScaffold = keys %{$genes2b{$scaffold}};
	push (@genes2, @genesOnScaffold);
}

open (STATS, ">$prefix.stats") or die $usage;
open (STATSM, ">$prefix.statsm") or die $usage;
open (STATSC, ">$prefix.statsc") or die $usage;

##now proceed to reference gene by gene analyses (digging clusters and determining ortholog synteny)
my @sorted_scaffolds = sort @scaffolds1; # sort the scaffolds

foreach my $scaffold (@sorted_scaffolds) { 
	my @genesOnScaffold = keys %{$genes1b{$scaffold}};
	my @sortedGenesOnScaffold = sort {$genes1b{$scaffold}->{$a}->{"start"} <=> $genes1b{$scaffold}->{$b}->{"start"}} @genesOnScaffold;
	
	my %index;
	@index{@sortedGenesOnScaffold} = (0..$#sortedGenesOnScaffold);
	foreach my $gene (@sortedGenesOnScaffold) {
		my @flanking_genes;		
		my $index = $index{$gene};
		my $num_genes = scalar @sortedGenesOnScaffold;
		##upstream reference genome 		
		for (my $i=$index; $i <= $index+$max_clus_len; $i++) {
			if ($i < $num_genes && $i != $index) {
				my $flank_gene = $sortedGenesOnScaffold[$i];
				push (@flanking_genes, $flank_gene);
			}
		}
		##downstream reference genome
		for (my $i=$index; $i >= 0; $i--) {
			if ($i != $index && $index-$i <= $max_clus_len) {
				my $flank_gene = $sortedGenesOnScaffold[$i];	
				push (@flanking_genes, $flank_gene);
			}	
		}
		my $num_flanks = scalar @flanking_genes;		
		my @mcl_flanking_genes;
		##proceed with genes within cluster with length $clus_len
		##next if ($num_flanks != $clus_len-1);
		##push cluster_ids for each flanking gene on reference genome in large array
		foreach my $flanking_gene (@flanking_genes) {
			my $mcl_current = $genes1b{$scaffold}->{$flanking_gene}->{"mcl"};
			push (@mcl_flanking_genes, $mcl_current); 		
		}
		
		my @intersects;
		##obtain matching (orthologous) genes for gene under investigation
		my $mcl_select = $genes1b{$scaffold}->{$gene}->{"mcl"};
		my @mcl_genes = @{$mclToGene{$mcl_select}->{"member"}};
		my @mcl_ortholog_flanking_genes;
		if (scalar @mcl_genes >= 2) { 		
			my $len_mcl_genes = scalar @mcl_genes;
			##obtain ortholog(s) and get their flanking genes
			foreach my $ortholog (@mcl_genes) {
				next if ($ortholog ~~ @genes1); ##skip (for now) genes from current organism (namely the ref, first one in argument list)
				##proceed with ortholog from query organism (second in argument list)
				print STATS join ("\t", $gene, $num_flanks, $mcl_select, $len_mcl_genes, $ortholog),"\t";
				my $ortholog_scaffold = $scaffoldToGene2b{$ortholog};
				my @orthologGenesOnScaffold = keys %{$genes2b{$ortholog_scaffold}};  ###TO BE ADDED: exception when orthologous gene is 'lonely' (one gene one scaffold for crappy genomes)
				my @orthologSortedGenesOnScaffold = sort {$genes2b{$ortholog_scaffold}->{$a}->{"start"} <=> $genes2b{$ortholog_scaffold}->{$b}->{"start"}} @orthologGenesOnScaffold;
				my %ortholog_index;
				@ortholog_index{@orthologSortedGenesOnScaffold} = (0..$#orthologSortedGenesOnScaffold);
				my $ortholog_index = $ortholog_index{$ortholog};
				my @ortholog_flanking_genes;
				my $ortholog_num_genes = scalar @orthologSortedGenesOnScaffold;
				##upstream query genome
				for (my $i=$ortholog_index; $i <= $ortholog_index+$max_clus_len; $i++) {
					if ($i < $ortholog_num_genes && $i != $ortholog_index) {
						my $flank_gene = $orthologSortedGenesOnScaffold[$i];
						push (@ortholog_flanking_genes, $flank_gene);
					}
				}
				##downstream query genome
				for (my $i=$ortholog_index; $i >= 0; $i--) {
					if ($i != $ortholog_index && $ortholog_index-$i <= $max_clus_len) {
						my $flank_gene = $orthologSortedGenesOnScaffold[$i];	
						push (@ortholog_flanking_genes, $flank_gene);
					}	
				}
				my $num_ortholog_flanks = scalar @ortholog_flanking_genes;
				print STATS $num_ortholog_flanks ,"\t";
				##push cluster_ids for each flanking gene on query genome in large array
				foreach my $ortholog_flanking_gene (@ortholog_flanking_genes) {
					my $mcl_current = $genes2b{$ortholog_scaffold}->{$ortholog_flanking_gene}->{"mcl"};
					push (@mcl_ortholog_flanking_genes, $mcl_current); 		
				}
				##now check for union @mcl_flanking_genes @mcl_ortholog_flanking_genes
				my $lc = List::Compare->new(\@mcl_flanking_genes, \@mcl_ortholog_flanking_genes);
				my @intersection = $lc->get_intersection;
				my $length_intersection = scalar @intersection;
				print STATS $length_intersection,"\n";
				push (@intersects, $length_intersection);

			}
		}
		else {
			my $len_mcl_genes = scalar @mcl_genes;			
			print STATS join ("\t", $gene, $num_flanks, $mcl_select, $len_mcl_genes, "_single_"),"\n";
		}
	my $length_intersects = scalar @intersects;
	if ($length_intersects == 0) {
		my $add = 0;
		push (@intersects, $add);
	}
	my $max_intersection = max @intersects;
	$genes1->{$scaffold}->{$gene}->{"collinear"} = $max_intersection;
	print STATSM join ("\t", $gene, $genes1->{$scaffold}->{$gene}->{"collinear"}),"\n";
	}
}

##determine 'cluster score' >> summed 'collinear' score for cluster with length $clus_len
##sliding 'gene window' approach (x >> y) , (x+1 >> y+1) etc.
##only select clusters (central genes) with length = $clus_len (flanking + central/selected gene)
foreach my $scaffold (@sorted_scaffolds) { 
	my @genesOnScaffold = keys %{$genes1b{$scaffold}};
	my @sortedGenesOnScaffold = sort {$genes1b{$scaffold}->{$a}->{"start"} <=> $genes1b{$scaffold}->{$b}->{"start"}} @genesOnScaffold;
	
	my %index;
	@index{@sortedGenesOnScaffold} = (0..$#sortedGenesOnScaffold);
	foreach my $gene (@sortedGenesOnScaffold) {
		my @flanking_genes;		
		my $index = $index{$gene};
		my $num_genes = scalar @sortedGenesOnScaffold;
		##upstream reference genome 		
		for (my $i=$index; $i <= $index+$max_clus_len; $i++) {
			if ($i < $num_genes && $i != $index) {
				my $flank_gene = $sortedGenesOnScaffold[$i];
				push (@flanking_genes, $flank_gene);
			}
		}
		##downstream reference genome
		for (my $i=$index; $i >= 0; $i--) {
			if ($i != $index && $index-$i <= $max_clus_len) {
				my $flank_gene = $sortedGenesOnScaffold[$i];	
				push (@flanking_genes, $flank_gene);
			}	
		}
		my $num_flanks = scalar @flanking_genes;		
		next if ($num_flanks != $clus_len-1);
		#flanking genes stored in @flanking_genes
		#obtain for each flanking gene (+original) the 'collinear' score and sum it
		my @scores;		
		foreach my $flanked_gene (@flanking_genes) {
			my $score = $genes1->{$scaffold}->{$flanked_gene}->{"collinear"};
			push (@scores, $score);
		}
		my $current_score = $genes1->{$scaffold}->{$gene}->{"collinear"};
		push (@scores, $current_score);
		my $sum_score = sum(@scores);
		my $norm_score = $current_score * $sum_score;
		print STATSC join ("\t", $gene, $scaffold, $num_genes, $num_flanks, $sum_score, $norm_score),"\n";
	}
}


#print Dumper ($genes1);

sub parse_gff3 {
	my ($gff3) = shift;	
	my (%genes, %scaffoldToGene, %rnaToGene, %proteinToGene); 
	
	open (GFF3, $gff3) or die "can't open this file $gff3";
	my $currentGene = "";
	while (my $line = <GFF3>)
	{
		chomp $line;
		next if ($line =~ /#/); #skip commented lines
	
		my @parts = split ("\t", $line);
		if (@parts < 2) { next; }
	
		my $scaffold = $parts[0];
		my $category = $parts[1];
		my $type = $parts[2]; 
		my $start = $parts[3]; 
		my $stop = $parts[4];
		my $score = $parts[5];
		my $orientation = $parts[6]; 
		my $phase = $parts[7];
		my @attributes = split(/;/,$parts[8]);
	
		my %parsedAttributes;
		$parsedAttributes{"biotype"} = "protein_coding";

		#parse attributes
		foreach my $option (@attributes)
		{
			my @optionParts = split(/=/, $option);	
			$parsedAttributes{$optionParts[0]} = $optionParts[1];
		}
	
		if ($type eq "gene")
		{
			next if ($parsedAttributes{"biotype"} ne "protein_coding"); #skip genes that are not protein_coding		
			#new gene, create entry in %genes
			my $geneID = $parsedAttributes{"ID"};
			my $id = $parsedAttributes{"ID"};
			$scaffoldToGene{$geneID} = $scaffold;
			$genes{$scaffold}->{$geneID}->{"ID"} = $id;
			$genes{$scaffold}->{$geneID}->{"start"} = $start;
			$genes{$scaffold}->{$geneID}->{"stop"} = $stop;
			$genes{$scaffold}->{$geneID}->{"orientation"} = $orientation;
			$genes{$scaffold}->{$geneID}->{"phase"} = $phase;
			$genes{$scaffold}->{$geneID}->{"score"} = $score;
			$genes{$scaffold}->{$geneID}->{"category"} = $category;
		
			# $genes{$scaffold}->{$geneID}->{"type"} = "coding";
			if (defined $parsedAttributes{"Note"})
			{
				$genes{$scaffold}->{$geneID}->{"Note"} = $parsedAttributes{"Note"};
			}
			if (defined $parsedAttributes{"Name"})
			{
				my $old_gene_name = $parsedAttributes{"Name"};
				$old_gene_name =~ s/%20/_/g;
				my $new_gene_name = $old_gene_name;
				$genes{$scaffold}->{$geneID}->{"Name"} = $new_gene_name;
			}
			if (defined $parsedAttributes{"Alias"})
			{
				$genes{$scaffold}->{$geneID}->{"Alias"} = $parsedAttributes{"Alias"};
			}
			if (defined $parsedAttributes{"Description"})
			{
				$genes{$scaffold}->{$geneID}->{"Description"} = $parsedAttributes{"Description"};
			}
		
			$currentGene = $geneID;

		} elsif ($type eq "mRNA" || $type eq "transcript")
		{
			next if ($parsedAttributes{"biotype"} ne "protein_coding"); #skip genes that are not protein_coding		
			#new RNA gene, add data to %genes and %rnaToGene
			my $ID = $parsedAttributes{"ID"};
			my $Parent = $currentGene;
			$rnaToGene{$ID} = $Parent;
			my $protein_id;
			if (defined $parsedAttributes{"transcript_id"}) {
				$protein_id = $parsedAttributes{"transcript_id"};
				$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"transcript_id"} = $parsedAttributes{"transcript_id"};
				$proteinToGene{$protein_id} = $Parent;
			}
			else {
				$proteinToGene{$ID} = $Parent;		
			}
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"ID"} = $ID; #create placeholder, exons will be added later on
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"Name"} = $parsedAttributes{"Name"};
				
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"start"} = $start;
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"stop"} = $stop;
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"orientation"} = $orientation;
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"category"} = $category;
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"score"} = $score;
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"phase"} = $phase;
		
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"Parent"} = $currentGene;
			if (defined $parsedAttributes{"Note"})
			{
				$genes{$scaffold}->{$currentGene}->{"Note"} = $parsedAttributes{"Note"};
			}
		}else  #types of exons, CDS, 3'UTR, 5'UTR
		{	
			next if ($type ne "CDS" || $type ne "exon"); ##only work with exon and cds features		
			my $ParentRNA = $parsedAttributes{"Parent"};		
			next if (! defined $rnaToGene{$ParentRNA}); ##skip exons/cds/utr etc if parent RNA is not defined (e.g. if biotype ne protein_coding) ##assumes correctly sorted gff3 (parental features loaded PRIOR to subs)		
			my $gene = $rnaToGene{$ParentRNA};
			my $tr_feature;
			$tr_feature->{"start"} = $start;
			$tr_feature->{"stop"} = $stop;
			$tr_feature->{"phase"} = $phase;
			$tr_feature->{"score"} = $score;
			$tr_feature->{"category"} = $category;
			$tr_feature->{"ID"} = $parsedAttributes{"ID"};
		
		
			if ($type eq "CDS")
			{
				push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$ParentRNA}->{"CDS"}}, $tr_feature);
			} 
			elsif ($type eq "exon")
			{
				push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$ParentRNA}->{"exon"}}, $tr_feature);
			}
			elsif ($type eq "five_prime_UTR")
			{
				push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$ParentRNA}->{"five_prime_UTR"}}, $tr_feature);
			}
			elsif ($type eq "three_prime_UTR")
			{
				push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$ParentRNA}->{"three_prime_UTR"}}, $tr_feature);
			}
			elsif ($type eq "intron")
			{
				push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$ParentRNA}->{"intron	"}}, $tr_feature);
			}
		}
	}

	close(GFF3);
	#print Dumper (\%genes);
	return (\%genes, \%scaffoldToGene, \%rnaToGene, \%proteinToGene);
}
