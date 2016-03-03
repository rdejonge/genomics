## GFF3stats.pl

#!/usr/bin/perl

####
#Perl script by Ronnie de Jonge, ronnie.dejonge@gmail.com OR ronnie.dejonge@psb.ugent.be
#Created: March 13, 2014
#Last modified: March 18, 2014
#usage: perl file.gff3
#Depends on: File::Basename , Data::Dumper , Math::Round , Statistics::Basic
#Description: 
#Calculate stats on a GFF3 file; such as:
# number of genes
# number of transcripts
# number of exons
# number of CDS-exons
# length of concatenated CDS exons
# distributions of the above
# number of introns
# number of introns/gene
# length of introns
# distribution of intron length 
####

$| = 1;

use warnings;
use strict;
use Data::Dumper; #we can use print Dumper (\%hash_ref) to print the content of hash tables using the Data::Dumper module
use File::Basename;
use Math::Round;
use Statistics::Basic qw(:all);

###script usage:
my $usage = "\n\nusage:	$0 file.gff3 \n\n";


###Input arguments, file.gff3
my $gff3 = $ARGV[0] or die $usage;

my $core_annot_filename = basename($gff3);
my ($gffstatsfh, $genesfh, $intronsfh, $exonsfh, $CDSsfh, $five_prime_UTRsfh, $three_prime_UTRsfh);

open ($gffstatsfh, ">$core_annot_filename.GFF3stats") or die $!;
open ($genesfh, ">$core_annot_filename.genes") or die $!;
open ($intronsfh, ">$core_annot_filename.introns") or die $!;
open ($exonsfh, ">$core_annot_filename.exons") or die $!;
open ($CDSsfh, ">$core_annot_filename.CDSs") or die $!;
open ($five_prime_UTRsfh, ">$core_annot_filename.five_prime_UTRs") or die $!;
open ($three_prime_UTRsfh, ">$core_annot_filename.three_prime_UTRs") or die $!;

# stats interested in:
my $gene_count = 0;
my $mRNA_count = 0;
my $alt_spliced_gene_count = 0;
my $intron_containing_gene_count = 0;
my $single_exon_genes = 0;
my $unique_exon_count = 0;
my $unique_cds_count = 0;
my $unique_intron_count = 0;

my $alt_splice_diff_CDSs_count = 0;
my $diff_splice_CDS_count = 0;

my $geneCounter = 0;

my $num_intergenic_regions = 0;
my $sum_intergenic_lengths = 0;
my $sum_gene_lengths = 0;

my $sum_exon_lengths = 0;
my $sum_CDS_lengths = 0;

#Initialize some hashes
my %genes;
my %rnaToGene;

#loading the gene annotations from file.gff3
open (GFF3, $gff3) or die $usage;

my $type_count = 0;
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

my @genes; my @transcripts; 
my @genes_all_introns; my @genes_all_exons; my @genes_all_CDSs; my @genes_all_fpUTRs; my @genes_all_tpUTRs; my @genes_all; my @genes_all_rnas;
my $n_introns=0;
my $sum_intron_lengths=0;

my (@fpUTR_full_lengths, @tpUTR_full_lengths, @CDS_full_lengths);


#print Dumper (\%genes);

#analyse each gene, write output to tab-delimited genefile like so (take the longest transcript if multiple available)
#geneX	scaffoldX	#transcripts	longest-transcriptID	length_gene	length_transcript	#exons	#CDSs	#introns	#length_exons	#length_CDSs	total_length_CDS
#	length_introns(I,II,III,n)	#five_prime_utr_exons	length_five_prime_utr	#three_prime_utr_exons	length_three_prime_utr	distance_to_nearest_upstream_gene	distance_to_nearest_downstream_gene

#and write output of each transcript to tab-delimited transcript file like so:
#transcriptX	scaffoldX	geneX	#transcripts	#exons	#CDSs	#introns	#length_exons	#length_CDSs	total_length_CDS
#	length_introns(I,II,III,n)	#five_prime_utr_exons	length_five_prime_utr	#three_prime_utr_exons	length_three_prime_utr

print $genesfh join ("\t", "gene_id", "scaffold_id", "gene_length", "no. RNAs", "longestRNA_id", "length_longestRNA_id", "no. exons", "no. CDSs", "no. 5'UTRs", "no. 3'UTRs", "total exon length", "total CDS length", "total 5'UTR length", "total 3'UTR length", "5'UTR present", "3'UTR present", "a UTR present", "No. of Introns", "mean Intron length", "3pUTR-Intron"), "\n";

my @scaffolds = keys %genes;
my @sorted_scaffolds = sort @scaffolds;

my $n_scaffolds = @sorted_scaffolds; #no. of scaffolds (with genes)

foreach my $scaffold (@sorted_scaffolds) #do for each scaffold
{ 
	my @genesOnScaffold = keys %{$genes{$scaffold}};
	my $n_genes = @genesOnScaffold; #no. of genes on scaffold
		
	foreach my $gene (@genesOnScaffold) #do for each gene
	{
		my $gene_id = $genes{$scaffold}->{$gene}->{"ID"};

		push (@genes, $gene_id);		
		
		my $gene_start = $genes{$scaffold}->{$gene}->{"start"};
		my $gene_stop = $genes{$scaffold}->{$gene}->{"stop"};
		my $gene_length = $gene_stop - $gene_start + 1;
		my $gene_orient = $genes{$scaffold}->{$gene}->{"orientation"};
		push (@genes_all, $gene_length);

		my @rnas = keys %{$genes{$scaffold}->{$gene}->{"RNA"}};
		my $gene_n_rnas = @rnas;
		my @gene_rna_lengths;
		my $longest_rna = 0;
		my @current_transcripts; #initialize current transcript array
		foreach my $rna (@rnas)
		{
			push (@transcripts, $rna); #add to big transcripts array
			push (@current_transcripts, $rna); #push to current transcript array
			my @sorted_current_transcripts = sort @current_transcripts; #sort current transcript array
			my $rna_length = $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"length"};
						
			if ($rna_length > $longest_rna)
			{ #if this rna is longer, we set it as the longest and change the length of the $longest_rna variable
				$longest_rna = $rna_length;
				$genes{$scaffold}->{$gene}->{"longest_rna"} = $rna;
			}
			elsif ($rna_length == $longest_rna)
			{ #if the length is identical: take first one alphabetically (no need to update $longest_rna variable)
				$genes{$scaffold}->{$gene}->{"longest_rna"} = $sorted_current_transcripts[0];				
			}
			else
			{ #if shorter; do nothing
			}		
			push (@gene_rna_lengths, $rna_length);
			push (@genes_all_rnas, $rna_length);
		}
		my $longest_rnaID = $genes{$scaffold}->{$gene}->{"longest_rna"};
		my @sorted_gene_rna_lengths = reverse sort { $a <=> $b } @gene_rna_lengths;

		#initialize some variables
		my @exons; my @CDSs; my @fpUTRs; my @tpUTRs;
		my @exon_lengths; my @CDS_lengths; my @fpUTR_lengths; my @tpUTR_lengths;
		my ($fpUTRpresence, $tpUTRpresence, $UTRpresence);
		my @introns = ();
		
		
		##setting some start variables
		my $gene_n_rna_exons=0; my $gene_n_rna_CDSs=0; my $gene_n_rna_fpUTR=0; my $gene_n_rna_tpUTR=0;
		my $exon_length=0; my $CDS_length=0; my $fpUTR_length=0; my $tpUTR_length=0;

		@exons = @{$genes{$scaffold}->{$gene}->{"RNA"}->{"$longest_rnaID"}->{"exon"}};
		$gene_n_rna_exons = @exons;
		
		@CDSs = @{$genes{$scaffold}->{$gene}->{"RNA"}->{"$longest_rnaID"}->{"CDS"}}; 
		$gene_n_rna_CDSs = @CDSs;
		
		#we don't always have UTRs, so only make updates/changes to @fpUTRs/@tpUTRs when we actually have them defined in the %genes
		if (defined $genes{$scaffold}->{$gene}->{"RNA"}->{"$longest_rnaID"}->{"five_prime_UTR"})
		{
			@fpUTRs = @{$genes{$scaffold}->{$gene}->{"RNA"}->{"$longest_rnaID"}->{"five_prime_UTR"}}; 
			$gene_n_rna_fpUTR = @fpUTRs;
		}
		
		if (defined $genes{$scaffold}->{$gene}->{"RNA"}->{"$longest_rnaID"}->{"three_prime_UTR"})
		{
			@tpUTRs = @{$genes{$scaffold}->{$gene}->{"RNA"}->{"$longest_rnaID"}->{"three_prime_UTR"}}; 
			$gene_n_rna_tpUTR = @tpUTRs;
		}
		
		
		
		foreach my $exon (@exons)
		{
			my $exon_length = $exon->{"stop"} - $exon->{"start"} + 1;
			push (@exon_lengths, $exon_length);
			push (@genes_all_exons, $exon_length);
			print $exonsfh join ("\t", "exon", $gene_id, $longest_rnaID, $scaffold, $gene_orient, $exon->{"start"}, $exon->{"stop"}, $exon_length),"\n";
		}
		foreach my $CDS (@CDSs)
		{
			my $CDS_length = $CDS->{"stop"} - $CDS->{"start"} + 1;
			push (@CDS_lengths, $CDS_length);
			push (@genes_all_CDSs, $CDS_length);
			print $CDSsfh join ("\t", "CDS", $gene_id, $longest_rnaID, $scaffold, $gene_orient, $CDS->{"start"}, $CDS->{"stop"}, $CDS_length),"\n";
		}
		foreach my $fpUTR (@fpUTRs)
		{
			my $fpUTR_length = $fpUTR->{"stop"} - $fpUTR->{"start"} + 1;
			push (@fpUTR_lengths, $fpUTR_length);
			push (@genes_all_fpUTRs, $fpUTR_length);
			print $five_prime_UTRsfh join ("\t", "five_prime_UTR", $gene_id, $longest_rnaID, $scaffold, $gene_orient, $fpUTR->{"start"}, $fpUTR->{"stop"}, $fpUTR_length),"\n";
		}
		foreach my $tpUTR (@tpUTRs)
		{
			my $tpUTR_length = $tpUTR->{"stop"} - $tpUTR->{"start"} + 1;
			push (@tpUTR_lengths, $tpUTR_length);
			push (@genes_all_tpUTRs, $tpUTR_length);
			print $three_prime_UTRsfh join ("\t", "three_prime_UTR", $gene_id, $longest_rnaID, $scaffold, $gene_orient, $tpUTR->{"start"}, $tpUTR->{"stop"}, $tpUTR_length),"\n";
		}
		
		foreach my $exon_lengthh (@exon_lengths)
		{
			$exon_length = $exon_length + $exon_lengthh;
		}
		foreach my $CDS_lengthh (@CDS_lengths)
		{
			$CDS_length = $CDS_length + $CDS_lengthh;
		}	
		foreach my $fpUTR_lengthh (@fpUTR_lengths)
		{
			$fpUTR_length = $fpUTR_length + $fpUTR_lengthh;
		}		
		foreach my $tpUTR_lengthh (@tpUTR_lengths)
		{
			$tpUTR_length = $tpUTR_length + $tpUTR_lengthh;
		}
		
		push (@CDS_full_lengths, $CDS_length);	#for full length statistics	
		
		#do we have five prime (fp) or three prime (tp) UTRs?
		#if we do, set the presence switches AND add the length to the @___full_lengths arrays for final statistics)
		if ($fpUTR_length > 0) {$fpUTRpresence = 1; push (@fpUTR_full_lengths, $fpUTR_length);} else {$fpUTRpresence = 0;}
		if ($tpUTR_length > 0) {$tpUTRpresence = 1; push (@tpUTR_full_lengths, $tpUTR_length);} else {$tpUTRpresence = 0;}
		if (($fpUTR_length > 0) || ($tpUTR_length > 0)) {$UTRpresence = 1;} else {$UTRpresence = 0;}
		
		#number of introns is the number of exons - 1
		my $gene_n_rna_introns = $gene_n_rna_exons - 1;
		my $gene_avgIntron_length;
		
		#total intronic length of a gene/transcript = length_mRNA(gene) - length_summed_exons(transcript)
		my $gene_intron_length = $longest_rna - $exon_length;
		if ($gene_n_rna_introns > 0)
		{
			#this gene/transcript has an intron (or more)
			$gene_avgIntron_length = $gene_intron_length/$gene_n_rna_introns;
			$gene_avgIntron_length = nearest(0.01, $gene_avgIntron_length ); #rounding the intron density for tab format clarity
			
			#finding the coordinates of the introns, by examining start/stop in exons and start/stop transcript
			#I use a temporary hash to store the exon data, sort it, and get the gaps == introns
			my %exons;
			my @sorted_exons;
			
			#get start and stop coordinates from the extremeties of the CDS coordinates
			my %cdss;
			my @sorted_CDS;
			
			foreach my $cds (@CDSs)
			{
				my $cds_length = $cds->{"stop"} - $cds->{"start"} + 1;
				my $cds_start = $cds->{"start"};
				my $cds_stop  = $cds->{"stop"};
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
	
			#print join ("\t", $first_cds, $last_cds, $start_coding, $stop_coding),"\n";
								
			foreach my $exon (@exons)
			{
				my $exon_length = $exon->{"stop"} - $exon->{"start"} + 1;
				my $exon_start = $exon->{"start"};
				my $exon_stop  = $exon->{"stop"};
				my $exon_id = $exon->{"ID"};
				
				
				$exons{$exon_id}->{"ID"} = $exon->{"ID"};
				$exons{$exon_id}->{"start"} = $exon->{"start"};
				$exons{$exon_id}->{"stop"} = $exon->{"stop"};
			}
			
			#now sort the hash using the start values
			foreach my $exon_id (sort {$exons{$a}->{"start"} <=> $exons{$b}->{"start"}} keys %exons)
			{
				#print "$exon_id\n";
				push (@sorted_exons, $exon_id);
			}

			my $num_exons = $#sorted_exons + 1;
			if ($gene_orient eq "+")
			{
				#gene on positive strand
				my $first_exon = shift @sorted_exons;

				
				while (@sorted_exons)
				{
					my $next_exon = shift @sorted_exons;
					my $first_end5 = $exons{$first_exon}->{"start"};
					my $first_end3 = $exons{$first_exon}->{"stop"};
					my $next_end5 = $exons{$next_exon}->{"start"};
					my $next_end3 = $exons{$next_exon}->{"stop"};
					my $intron_end5 = $first_end3 + 1;
					my $intron_end3 = $next_end5 - 1;
					if ($intron_end5 < $intron_end3)
					{			
						if ($intron_end5 < $start_coding)
						{
							#upstream coding, thus 5'UTR
							push (@introns, [$intron_end5, $intron_end3, "fpUTR"]);
						}
						elsif ( ($intron_end5 > $start_coding) && ($intron_end5 < $stop_coding) )
						{
							#within CDS, thus coding intron
							push (@introns, [$intron_end5, $intron_end3, "CDS"]);
						}
						elsif ($intron_end5 > $stop_coding)
						{
							#downstream coding, thus 3'UTR
							push (@introns, [$intron_end5, $intron_end3, "tpUTR"]);
							$genes{$scaffold}->{$gene}->{"RNA"}->{"$longest_rnaID"}->{"tpUTR_intron"} = "yes";							
						}
					}
					$first_exon = $next_exon;
				}
			}
			elsif ($gene_orient eq "-")
			{
				#gene on negative strand
				my $first_exon = shift @sorted_exons;
				while (@sorted_exons)
				{
					my $next_exon = shift @sorted_exons;
					my $first_end5 = $exons{$first_exon}->{"start"}; 
					my $first_end3 = $exons{$first_exon}->{"stop"};
					my $next_end5 = $exons{$next_exon}->{"start"}; 
					my $next_end3 = $exons{$next_exon}->{"stop"};
					my $intron_end5 = $first_end3 + 1; 
					my $intron_end3 = $next_end5 - 1;
					if ($intron_end3 > $intron_end5) 
					{
						if ($intron_end3 > $stop_coding)
						{
							#upstream coding, thus 5'UTR
							push (@introns, [$intron_end5, $intron_end3, "fpUTR"]);
						}
						elsif ( ($intron_end3 < $stop_coding) && ($intron_end3 > $start_coding) )
						{
							#within CDS, thus coding intron
							push (@introns, [$intron_end5, $intron_end3, "CDS"]);
						}
						elsif ($intron_end3 < $start_coding)
						{
							#downstream coding, thus 3'UTR
							push (@introns, [$intron_end5, $intron_end3, "tpUTR"]);
							
							#adding a flag here to tell the system that the current gene has a 3' UTR intron 
							$genes{$scaffold}->{$gene}->{"RNA"}->{"$longest_rnaID"}->{"tpUTR_intron"} = "yes";
						}
					}
					$first_exon = $next_exon;
				}
			}
			
			my $intron_n = 1;
			
			foreach my $intron_coordset (@introns) 
			{
				my $intron_token = join ("_", @$intron_coordset);
				my ($intron_lend, $intron_rend, $type) = split (/_/, $intron_token);
				print $intronsfh join ("\t", "intron", $gene_id, $longest_rnaID, $scaffold, $gene_orient, $intron_lend, $intron_rend, $intron_n, $type),"\n"; 
				my $intron_len = $intron_rend - $intron_lend + 1;
				push (@genes_all_introns, $intron_len);
				$sum_intron_lengths += $intron_len;
				++$n_introns;
				++$intron_n;
			}

			#print Dumper (\%exons);
		}
		else
		{
			$gene_avgIntron_length = "-";
		}
		
		#I'll calculate the average length of the introns in a gene, as well as for all genes together (using @genes_all_introns); same for exons (@genes_all_exons), CDSs (@genes_all_CDSs), fpUTRs (@genes_all_fpUTRs) and tpUTRs (@genes_all_tpUTRs)	
		
		my $tpUTR_intronFlag;
		if (defined $genes{$scaffold}->{$gene}->{"RNA"}->{"$longest_rnaID"}->{"tpUTR_intron"})
		{ $tpUTR_intronFlag = 1; } else{ $tpUTR_intronFlag = 0; }
		
		print $genesfh join ("\t", $gene_id, $scaffold, $gene_length, $gene_n_rnas, $longest_rnaID, $longest_rna, $gene_n_rna_exons, $gene_n_rna_CDSs, $gene_n_rna_fpUTR, $gene_n_rna_tpUTR, $exon_length, $CDS_length, $fpUTR_length, $tpUTR_length, $fpUTRpresence, $tpUTRpresence, $UTRpresence, $gene_n_rna_introns, $gene_avgIntron_length, $tpUTR_intronFlag), "\n";
	}
}
		
###now that the output is written to $input.genes we will calculate some overall statistics and report them

my $n_genes = @genes;
my $n_transcripts = @transcripts;
my $n_exons = @genes_all_exons;
my $n_CDSs = @genes_all_CDSs;
my $n_fpUTRs = @genes_all_fpUTRs;
my $n_tpUTRs = @genes_all_tpUTRs;
my $n_full_CDSs = @CDS_full_lengths;
my $n_full_fpUTRs = @fpUTR_full_lengths;
my $n_full_tpUTRs = @tpUTR_full_lengths;

my ($mean_genes, $mean_transcripts, $mean_exons, $mean_CDSs, $mean_fpUTRs, $mean_tpUTRs, $mean_introns) = (0,0,0,0,0,0,0);
my ($median_genes, $median_transcripts, $median_exons, $median_CDSs, $median_fpUTRs, $median_tpUTRs, $median_introns) = (0,0,0,0,0,0,0);


$mean_genes = mean(@genes_all); $median_genes = median(@genes_all);
$mean_transcripts = mean(@genes_all_rnas); $median_transcripts = median(@genes_all_rnas);
$mean_exons = mean(@genes_all_exons); $median_exons = median(@genes_all_exons);
$mean_CDSs = mean(@genes_all_CDSs); $median_CDSs = median(@genes_all_CDSs);
$mean_fpUTRs = mean(@genes_all_fpUTRs); $median_fpUTRs = median(@genes_all_fpUTRs);
$mean_tpUTRs = mean(@genes_all_tpUTRs); $median_tpUTRs = median(@genes_all_tpUTRs);
$mean_introns = mean(@genes_all_introns); $median_introns = median(@genes_all_introns);

my ($mean_CDS_full, $mean_fpUTRs_full, $mean_tpUTRs_full) = (0,0,0);
my ($median_CDS_full, $median_fpUTRs_full, $median_tpUTRs_full) = (0,0,0);
$mean_CDS_full = mean(@CDS_full_lengths); $median_CDS_full = median(@CDS_full_lengths);
$mean_fpUTRs_full = mean(@fpUTR_full_lengths); $median_fpUTRs_full = median(@fpUTR_full_lengths);
$mean_tpUTRs_full = mean(@tpUTR_full_lengths); $median_tpUTRs_full = median(@tpUTR_full_lengths);



print $gffstatsfh "GFF3 overall stats","\t",$gff3,"\n\n";
print $gffstatsfh "Total number of ..","\n";
print $gffstatsfh "..genes","\t",$n_genes,"\n";
print $gffstatsfh "..mRNAs","\t",$n_transcripts,"\n";
print $gffstatsfh "..exons","\t",$n_exons,"\n";
print $gffstatsfh "..CDSs frag","\t",$n_CDSs,"\n";
print $gffstatsfh "..5'UTRs frag","\t",$n_fpUTRs,"\n";
print $gffstatsfh "..3'UTRs frag","\t",$n_tpUTRs,"\n";
print $gffstatsfh "..CDSs","\t",$n_full_CDSs,"\n";
print $gffstatsfh "..5'UTRs","\t",$n_full_fpUTRs,"\n";
print $gffstatsfh "..3'UTRs","\t",$n_full_tpUTRs,"\n";
print $gffstatsfh "..introns","\t",$n_introns,"\n\n";
print $gffstatsfh "Mean|Median length of ..","\n";
print $gffstatsfh "..genes","\t",$mean_genes,"\t",$median_genes,"\n";
print $gffstatsfh "..mRNAs","\t",$mean_transcripts,"\t",$median_transcripts,"\n";
print $gffstatsfh "..exons","\t",$mean_exons,"\t",$median_exons,"\n";
print $gffstatsfh "..CDSs","\t",$mean_CDSs,"\t",$median_CDSs,"\n";
print $gffstatsfh "..5'UTRs","\t",$mean_fpUTRs,"\t",$median_fpUTRs,"\n";
print $gffstatsfh "..3'UTRs","\t",$mean_tpUTRs,"\t",$median_tpUTRs,"\n";
print $gffstatsfh "..introns","\t",$mean_introns,"\t",$median_introns,"\n\n";
print $gffstatsfh "Mean|Median length of full-length ..","\n";
print $gffstatsfh "..CDSs","\t",$mean_CDS_full,"\t",$median_CDS_full,"\n";
print $gffstatsfh "..5'UTRs","\t",$mean_fpUTRs_full,"\t",$median_fpUTRs_full,"\n";
print $gffstatsfh "..3'UTRs","\t",$mean_tpUTRs_full,"\t",$median_tpUTRs_full,"\n";

