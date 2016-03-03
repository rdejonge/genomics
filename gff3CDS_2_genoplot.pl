#!/usr/bin/perl
####
#By Ronnie de Jonge, 
#Created: Feb 6, 2015 
#Last modified: Feb 6, 2015 
#Usage: perl --genes genes.gff3 --target target_genes.list 
#Extracting CDS start/stop for genoPlotR.R plotting (annotations file)
####

use warnings; 
use strict; 
use Data::Dumper; 
use Getopt::Long; 
use experimental 'smartmatch'; 

my $usage = "\n\nusage: $0
				--genes <genes.gff3>
				--target <target_genes.list> \n\n"; if (@ARGV < 2) {
	print $usage;
	exit;
}

my ($gff3, $target); GetOptions("genes=s" => \$gff3,
		   "target=s" => \$target, );
##parsing genes.gff3 file
my %genes; my %rnaToGene;
#loading the gene annotations from file.gff3 (--genes option)
open (GFF3, $gff3) or die $usage; my $type_count = 0; while (my $line = <GFF3>) {
	chomp $line;
	next if ($line =~ /^#/); #skip commented lines
	
	my @parts = split ("\t", $line);
	if (@parts < 2) { next; }
	
	my $scaffold = $parts[0];
	my $source = $parts[1];
	my $type = $parts[2];
	my $start = $parts[3];
	my $stop = $parts[4];
	my $score = $parts[5];
	my $orientation = $parts[6];
	my $phase = $parts[7];
	my @attributes = split(/;/,$parts[8]);
	
	my %parsedAttributes;
	
	#parse attributes
	foreach my $option (@attributes)
	{
		my @optionParts = split(/=/, $option);
		$parsedAttributes{$optionParts[0]} = $optionParts[1];
	}
	
	if ($type eq "gene")
	{
		#new gene, create entry in %genes
		
		my $geneID = $parsedAttributes{"ID"};
		my $id = $parsedAttributes{"ID"};
		$genes{$scaffold}->{$geneID}->{"ID"} = $id;
		$genes{$scaffold}->{$geneID}->{"start"} = $start;
		$genes{$scaffold}->{$geneID}->{"stop"} = $stop;
		$genes{$scaffold}->{$geneID}->{"orientation"} = $orientation;
		if (defined $parsedAttributes{"Note"})
		{
			$genes{$scaffold}->{$geneID}->{"Note"} = $parsedAttributes{"Note"};
		}
		if (defined $parsedAttributes{"Name"})
		{
			$genes{$scaffold}->{$geneID}->{"Name"} = $parsedAttributes{"Name"};
		}
		if (defined $parsedAttributes{"Alias"})
		{
			$genes{$scaffold}->{$geneID}->{"Alias"} = $parsedAttributes{"Alias"};
		}
		if (defined $parsedAttributes{"Description"})
		{
			$genes{$scaffold}->{$geneID}->{"Description"} = $parsedAttributes{"Description"};
		}
		
	} elsif ($type eq "mRNA" || $type eq "transcript")
	{
		#new RNA gene, add data to %genes and %rnaToGene
		my $ID = $parsedAttributes{"ID"};
		my $Parent = $parsedAttributes{"Parent"};
		$rnaToGene{$ID} = $Parent;
		$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"ID"} = $ID;
		$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"Name"} = $ID;
		$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"start"} = $start;
		$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"stop"} = $stop;
		$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"orientation"} = $orientation;
		$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"length"} = $stop - $start + 1;
		
		$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"Parent"} = $Parent;
		if (defined $parsedAttributes{"Note"})
		{
			$genes{$scaffold}->{$Parent}->{"Note"} = $parsedAttributes{"Note"};
		}
	}elsif ($type eq "exon" || $type eq "CDS" || $type eq "five_prime_UTR" || $type eq "three_prime_UTR")
	{
		my $ParentRNA = $parsedAttributes{"Parent"};
		my $gene = $rnaToGene{$ParentRNA};
		my $tr_feature;
		$tr_feature->{"start"} = $start;
		$tr_feature->{"stop"} = $stop;
		
		if (defined $parsedAttributes{"ID"})
		{
			$tr_feature->{"ID"} = $parsedAttributes{"ID"};
		}
		else
		{
			$tr_feature->{"ID"} = $ParentRNA . ":" . $type . ":" . $type_count;
			++$type_count;
		}
		
		if ($type =~ /codon/ || $type =~ "splice") { next; }
		if ($type eq "CDS")
		{
			#CDS IDs are not always unique; I need them to be unique, so add specifying number
			$tr_feature->{"ID"} = "CDS:" . $parsedAttributes{"ID"} . ":" . $type_count;
			++$type_count;
			push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$ParentRNA}->{"CDS"}}, $tr_feature);
		} 
		elsif ($type eq "exon")
		{
			push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$ParentRNA}->{"exon"}}, $tr_feature);
		}
		elsif (($type eq "five_prime_UTR") || ($type eq "5'UTR"))
		{
			push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$ParentRNA}->{"five_prime_UTR"}}, $tr_feature);
		}		
		elsif (($type eq "three_prime_UTR") || ($type eq "3'UTR"))
		{
			push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$ParentRNA}->{"three_prime_UTR"}}, $tr_feature);
		}		
	}
	else{
		#other features >> ignore ...
	}
}
close(GFF3);
##extracting start/stop per gene
my @scaffolds = keys %genes; my @sorted_scaffolds = sort @scaffolds; my %transcripts; foreach my $scaffold (@sorted_scaffolds) {
	my @genesOnScaffold = keys %{$genes{$scaffold}};
	my $n_genes = @genesOnScaffold;
		
	foreach my $gene (@genesOnScaffold)
	{
		my @rnas = keys %{$genes{$scaffold}->{$gene}->{"RNA"}};
		foreach my $rna (@rnas)
		{
			my $strand = $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"orientation"};
			my @sorted_CDS; my @CDSs; my %cdss;
			@CDSs = @{$genes{$scaffold}->{$gene}->{"RNA"}->{"$rna"}->{"CDS"}};
			foreach my $cds (@CDSs)
			{
				my $cds_start = $cds->{"start"};
				my $cds_stop = $cds->{"stop"};
				my $cds_id = $cds->{"ID"};
				
				$cdss{$cds_id}->{"ID"} = $cds->{"ID"};
				$cdss{$cds_id}->{"start"} = $cds->{"start"};
				$cdss{$cds_id}->{"stop"} = $cds->{"stop"};
			}
			#sort the cds hash using start values
			foreach my $cds_id (sort {$cdss{$a}->{"start"} <=> $cdss{$b}->{"start"}} keys %cdss)
			{
				push (@sorted_CDS, $cds_id);
			}
			my $first_cds = $sorted_CDS[0];
			my $last_cds = $sorted_CDS[-1];
			my $start_coding = $cdss{$first_cds}->{"start"};
			my $stop_coding = $cdss{$last_cds}->{"stop"};
			$transcripts{$rna}->{"startCDS"} = $start_coding;
			$transcripts{$rna}->{"stopCDS"} = $stop_coding;
			$transcripts{$rna}->{"strandCDS"} = $strand;
		}
	}
}
print join ("\t", "name", "start", "stop", "strand", "col"),"\n";
##extracting target start/stops
open (TARGET, $target) or die $usage; while (my $targets = <TARGET>) {
	chomp $targets;
	print $targets,"\t";
	print join ("\t", $transcripts{$targets}->{"startCDS"},
						$transcripts{$targets}->{"stopCDS"},
						$transcripts{$targets}->{"strandCDS"},
						"black"
				),"\n";
}
