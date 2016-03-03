## projectGOterms.pl

#!/usr/bin/perl
#project GO terms from Aspergillus nidulans to C. beticola based on orthoMCL clusters 
#Ronnie de Jonge, VIB, BEG, 
#Created on August 14, 2014
#Last modified on August 14, 2014

use warnings;
use strict;
use Data::Dumper;

my $GO = $ARGV[0];
my $fam = $ARGV[1];
my $cbeList = $ARGV[2];
my $GOtargetGff3 = $ARGV[3];
my $prefixOut = $ARGV[4];

###script usage:
my $usage = "\n\nusage:	$0 <GO_set_FromAMIGO> <FamilyFile(mclConvertedOutput)> <ProteinListCBE> <GOtarget.gff3 (for name conversion)> <prefixForOutput> \n\n";

my @cbe; my %famRev; my %fam;
my %ids; my %ids_rev; 
my %go_terms;

&run_idConverter;

#parse cbeList
open (CL, $cbeList) or die $usage;
while (my $line = <CL>)
{
	chomp $line;
	push (@cbe, $line);
}
close (CL);
my @sort_cbe = sort @cbe;

#parseFam
open (FAM, $fam) or die $usage;
while (my $fam_line = <FAM>)
{
	chomp $fam_line;
	my ($protein_id, $famID) = split ("\t", $fam_line);
	push(@{$fam{$famID}->{"proteins"}}, $protein_id);
	$famRev{$protein_id}->{"family"} = $famID;
}
close (FAM);

#AspGD	ASPL0000117821	5S-rRNA		GO:0002181	AspGD_REF:ASPL0000082571|PMID:6453331	TAS		P		5S-rRNA|5S rRNA|5S RNA	gene_product	taxon:162425	20130211	AspGD		
#parse asp_go term file
open (GO, $GO) or die $usage;
while (my $go_line = <GO>)
{
	next if ($go_line =~ /^!/); #skip lines starting with a exclamation mark
	next if ($go_line =~ /^#/); #skip commented lines
	
	my @parts = split ("\t", $go_line);
	if (@parts < 2) { next; } #skip lines with less than 2 parts
	
	my $prot_id = $parts[2];
	my $go_term = $parts[4];
	my $tax_id = $parts[12];
	
	if ($tax_id eq "taxon:162425")
	{
		#filter for target species; Aspergillus nidulans
		#print join ("\t", $prot_id, $go_term, $tax_id),"\n";
		if (defined $ids_rev{$prot_id})
		{
			my $mRNAid = $ids_rev{$prot_id}->{"rna_id"};
			push(@{$go_terms{$mRNAid}->{"go_terms"}}, $go_term);		
		}
	}
	else
	{
	}	
}
close (GO);

#print Dumper (\%go_terms);

#parse target species gene/protein list and print for each protein
#1) ID
#2) fam (orthoMCL clusters)
#3) Aspergillus protein(s) in fam (if they are present)
#4) GO-terms of Aspergillus proteins from fam (if they are present)

open (BYGO, ">$prefixOut.byGO.list") or die $usage;
open (BYPROT, ">$prefixOut.byPROT.list") or die $usage;

open (LIST, $cbeList) or die $usage;
while (my $list_line = <LIST>)
{
	my @asp_prots; my @goTerms; my $famProt=""; my @famProts;
	chomp $list_line;
	my $prot_id = $list_line;
	if (defined $famRev{$prot_id})
	{
		$famProt = $famRev{$prot_id}->{"family"};
		@famProts = @{$fam{$famProt}->{"proteins"}};
		
		foreach my $famProt (@famProts)
		{
			if ($famProt =~ /ASP/)
			{
				push (@asp_prots, $famProt);
			}
		}
		foreach my $asp_prot (@asp_prots)
		{
			if (defined $go_terms{$asp_prot}->{"go_terms"})
			{
				@goTerms = @{$go_terms{$asp_prot}->{"go_terms"}};
				foreach my $goTerm (@goTerms)
				{
					print BYGO join ("\t", $goTerm, "ISA", $prot_id),"\n";
				}
			}
		}
	}
	print BYPROT join ("\t", $prot_id, $famProt, "@asp_prots", "@goTerms"),"\n";
}
close (LIST);
close (BYPROT);
close (BYGO);


sub run_idConverter
{
	#the go-term database use the original AN identifiers
	#these identifiers are stored in my parse gff3 files (that were used to build protein sequences, compliantFasta db and consequently orthoMCL output)
	#now get a converter list and store in hash
	#need to parse gff3 and get two things:
	#1) geneID (new) with name (oldId, targetID!)
	#2) mRNAid > parent geneID 
	#The final list should contain a hash of oldGeneName linked to new_mRNAids
	my %genes;
	my %rnaToGene;
	
	open (GFF3, $GOtargetGff3) or die $usage;
	while (my $line = <GFF3>)
	{
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
			my $geneID = $parsedAttributes{"ID"};
			my $id = $parsedAttributes{"ID"};
			$genes{$scaffold}->{$geneID}->{"ID"} = $id;
			$genes{$scaffold}->{$geneID}->{"start"} = $start;
			$genes{$scaffold}->{$geneID}->{"stop"} = $stop;
			$genes{$scaffold}->{$geneID}->{"orientation"} = $orientation;
			if (defined $parsedAttributes{"Name"})
			{
				$genes{$scaffold}->{$geneID}->{"Name"} = $parsedAttributes{"Name"};
			}

			
		} 
		elsif ($type eq "mRNA")
		{
			#new RNA gene, add data to %genes and %rnaToGene
			my $ID = $parsedAttributes{"ID"};
			my $Parent = $parsedAttributes{"Parent"};
			$rnaToGene{$ID} = $Parent;
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"ID"} = $ID;					
			$genes{$scaffold}->{$Parent}->{"RNA"}->{$ID}->{"Parent"} = $Parent;
		}
		else
		{#do nothing
		}
	}	
	#the dat is now stored in %genes and %rnaToGene	
	#now get the useful parts and store in %ids
	
	my @scaffolds = keys %genes;
	foreach my $scaffold (@scaffolds)
	{
		my @genesOnScaffold = keys %{$genes{$scaffold}};
		foreach my $geneOnScaffold (@genesOnScaffold)
		{
			my $gene_id = $genes{$scaffold}->{$geneOnScaffold}->{"ID"};
			my $gene_name = $genes{$scaffold}->{$geneOnScaffold}->{"Name"};
			my @rnas = keys %{$genes{$scaffold}->{$geneOnScaffold}->{"RNA"}};
			foreach my $rna (@rnas)
			{
				my $rna_id = $genes{$scaffold}->{$geneOnScaffold}->{"RNA"}->{$rna}->{"ID"};
				#print $rna_id,"\t",$gene_name,"\n";
				$ids{$rna_id}->{"oldName"} = $gene_name;
				$ids_rev{$gene_name}->{"rna_id"} = $rna_id;
			}
		}
	}
}



