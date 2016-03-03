## JGI-GTF_2_GFF3.pl

#!/usr/bin/perl

####
#Perl script by Ronnie de Jonge, ronnie.dejonge@gmail.com OR ronnie.dejonge@psb.ugent.be
#Created: March 26, 2014
#Last modified: April 4, 2014
#usage: perl JGI_TO_GFF3_1.pl [jgi gff_input] [prefix]
#
#
#Description: 
#The script will first store each gene model in a hash table and then print each gene model
#Adds UTR information using _addUTR subroutine (works on the hash structure with CDS and EXON arrays)
#Adds alternative splice variants (if they exist; and based on the combination 'gene name' <> 'transcriptId' (for exons) or <> 'proteinId' (for CDSs) )
####

use warnings;
use strict;
use Data::Dumper;
use File::Basename;

=item input JGI gff

stop_codon is contained in last CDS (simply skipping these lines)

scaffold_1	JGI	exon	4845	6785	.	-	.	name "estExt_fgenesh1_kg.C_1_t10001"; transcriptId 67451
scaffold_1	JGI	CDS	5450	6571	.	-	0	name "estExt_fgenesh1_kg.C_1_t10001"; proteinId 67409; exonNumber 1
scaffold_1	JGI	start_codon	6569	6571	.	-	0	name "estExt_fgenesh1_kg.C_1_t10001"
scaffold_1	JGI	stop_codon	5450	5452	.	-	0	name "estExt_fgenesh1_kg.C_1_t10001"

to: (where Dotse1_ is a set prefix by $ARGV[1] and t1 is the first mRNA for a particular gene (subsequent transcripts will be .t2 .t3 etc..))

scaffold_1	JGI_MOD	gene	4845	6785	.	-	.	ID=Dotse1_67451;Name=estExt_fgenesh1_kg.C_1_t10001;
scaffold_1	JGI_MOD	mRNA	5450	6571	.	-	.	ID=Dotse1_67451.t1;Parent=Dotse1_67451;Alias=Dotse1_pID_67409;
scaffold_1	JGI_MOD	exon	4845	6785	.	-	.	ID=exon:1:Dotse1_67409.t1;Parent=Dotse1_67409.t1;
scaffold_1	JGI_MOD	CDS 	5450	6571	.	-	0	ID=CDS:1:Dotse1_67409.t1;Parent=Dotse1_67409.t1;

Important notices:
The name should remain in the file (as addition). 
proteinId will be added as alias to mRNA track
transcriptId will serve as second part of the new ID

=cut

###script usage:
my $usage = "\n\nusage:	$0 jgi.gff genemodel_ID \n\n";


###Input arguments, GFF3_IN and GFF3_OUT 
my $gff = $ARGV[0] or die $usage;
my $core_annot_filename = basename($gff);

my $gff3_out = $core_annot_filename . ".gff3";
my $idprefix = $ARGV[1] or die $usage;
my %genes;
my $countTotalExon = 1;
my $countTotalCDS = 1;
my $countGenes = 1;
my $countTranscripts = 1;
my $countTotalUTR = 1;
my $countTotal5UTR = 1;
my $countTotal3UTR = 1;
my $tn = -1;

open (GFF, $gff) or die "can't open this file $gff";

#first check all exon features; in the second round we look at CDS (to prevent problems with initiation and more)
while (my $line = <GFF>)
{
	chomp $line;
	next if ($line =~ /^#/); #skip commented lines
	
	my @parts = split ("\t", $line);
	if (@parts < 2) { next; }	
	
	my $scaffold = $parts[0];
	my $source = $parts[1]; 
	my $new_source = $source . "_MOD";
	my $type = $parts[2]; 
	my $start = $parts[3]; 
	my $stop = $parts[4];
	my $score = $parts[5];
	my $orientation = $parts[6]; 
	my $phase = $parts[7];
	my @attributes = split(/; /,$parts[8]);	
	
	my %parsedAttributes;
	
	#parse attributes
	foreach my $option (@attributes)
	{
		my @optionParts = split(/ /, $option);
		$parsedAttributes{$optionParts[0]} = $optionParts[1];
	}	
	
	if (($type eq "startcodon") || ($type eq "stop_codon") || ($type eq "CDS"))
	{
		#do nothing, skip these features (for now)
	}
	
	elsif ($type eq "exon")
	{	
		#a exon from a gene >> find gene by "name attribute" lookup; if it doesn't exist yet create it.
		my $geneName = $parsedAttributes{"name"};
		$geneName =~ s/\"//g;
		my $geneID = $geneName;
		my $tid = $parsedAttributes{"transcriptId"};
		my $new_tid = $idprefix . "_" . $tid . "T";
		my $geneIDnew = $idprefix . "_" . $countGenes;
		my $tr_feature;
		if (defined $genes{$scaffold}->{$geneID}->{"ID"})
		{#this gene is in the database, let's add the current line of information to it
			if ($start < $genes{$scaffold}->{$geneID}->{"start"})
			{
				$genes{$scaffold}->{$geneID}->{"start"} = $start;
			}
			if ($stop > $genes{$scaffold}->{$geneID}->{"stop"})
			{
				$genes{$scaffold}->{$geneID}->{"stop"} = $stop;
			}			
			
			my $n_rnas = scalar keys %{$genes{$scaffold}->{$geneID}->{"RNA"}};
			
			$tn = $n_rnas + 1;

			$new_tid = $new_tid . $tn;
			print $geneID,"\t",$n_rnas,"\t",$tn,"\t",$new_tid,"\n";

			if (defined $genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"ID"})
			{#this rna (transcriptId) is in the database; does add the information also here
				
				#update the gene/mRNA borders
				if ($start < $genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"start"})
				{
					$genes{$scaffold}->{$geneID}->{"RNA"}->{$tid}->{"start"} = $start;
				}
				if ($stop > $genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"stop"})
				{
					$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"stop"} = $stop;
				}
			}
			else
			{#this rna (transcriptId) is new to the database, thus make it
				
				$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"ID"} = $tid;
				$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"start"} = $start;
				$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"stop"} = $stop;
				$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"new_name"} = $new_tid;
				$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"old_name"} = $tid;
				++$countTranscripts;
			}
				
				
			
			#add exon information (to the transcript)
			$tr_feature->{"start"} = $start;
			$tr_feature->{"stop"} = $stop;
			$tr_feature->{"phase"} = $phase;
			$tr_feature->{"score"} = $score;
			$tr_feature->{"ID"} = "exon:" . $new_tid . ":" . $countTotalExon;
			++$countTotalExon;
			
			push(@{$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"exon"}}, $tr_feature);		
			
		}
		else
		{
			#this gene is new to the database
			#first: initiate it and add the current line of information to it
			$tn = -1;
			
			$genes{$scaffold}->{$geneID}->{"ID"} = $geneID;
			$genes{$scaffold}->{$geneID}->{"start"} = $start;
			$genes{$scaffold}->{$geneID}->{"stop"} = $stop;
			$genes{$scaffold}->{$geneID}->{"orientation"} = $orientation;
			$genes{$scaffold}->{$geneID}->{"phase"} = $phase;
			$genes{$scaffold}->{$geneID}->{"score"} = $score;
			$genes{$scaffold}->{$geneID}->{"source"} = $source;
			$genes{$scaffold}->{$geneID}->{"new_source"} = $new_source;
			$genes{$scaffold}->{$geneID}->{"name"} = $geneName;
			$genes{$scaffold}->{$geneID}->{"new_ID"} = $geneIDnew;
			
			++$countGenes;
			
			#add transcript information
			$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"ID"} = $tid;
			$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"start"} = $start;
			$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"stop"} = $stop;
			$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"new_name"} = $new_tid;
			$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"old_name"} = $tid;
			++$countTranscripts;
			
			#add exon information (to the transcript)
			$tr_feature->{"start"} = $start;
			$tr_feature->{"stop"} = $stop;
			$tr_feature->{"phase"} = $phase;
			$tr_feature->{"score"} = $score;
			$tr_feature->{"ID"} = "exon:" . $new_tid . ":" . $countTotalExon;
			++$countTotalExon;
			
			push(@{$genes{$scaffold}->{$geneID}->{"RNA"}->{$new_tid}->{"exon"}}, $tr_feature);		
		}
	}
}
	
close(GFF);

open (GFF, $gff) or die "can't open this file $gff";

#EXON features are already loaded; in this second round I append CDS information
while (my $line = <GFF>)
{
	chomp $line;
	next if ($line =~ /^#/); #skip commented lines
	
	my @parts = split ("\t", $line);
	if (@parts < 2) { next; }	
	
	my $scaffold = $parts[0];
	my $source = $parts[1]; 
	my $new_source = $source . "_MOD";
	my $type = $parts[2]; 
	my $start = $parts[3]; 
	my $stop = $parts[4];
	my $score = $parts[5];
	my $orientation = $parts[6]; 
	my $phase = $parts[7];
	my @attributes = split(/; /,$parts[8]);	
	
	my %parsedAttributes;
	
	#parse attributes
	foreach my $option (@attributes)
	{
		my @optionParts = split(/ /, $option);
		$parsedAttributes{$optionParts[0]} = $optionParts[1];
	}	
	
	if (($type eq "startcodon") || ($type eq "stop_codon") || ($type eq "exon"))
	{
		#do nothing
	}

	elsif ($type eq "CDS")
	{
		#a CDS from a gene >> find gene by "name attribute" lookup (if it doesn't exist, print warning message and drop)
		my $geneName = $parsedAttributes{"name"};
		$geneName =~ s/\"//g;
		my $geneID = $geneName;
		my $id = $geneID;
		my $pid = $parsedAttributes{"proteinId"};
		my $tr_feature;
		if (defined $genes{$scaffold}->{$geneID}->{"ID"})
		{
			#print "warning: gene exists; add CDS information to it","\t",$geneName,"\n";
			
			if (! defined $genes{$scaffold}->{$geneID}->{"proteinId"})
			{
				#print "warning: proteinId is new, add this information to the genelevel, add CDS information to the mRNA level/exons","\n";
				$genes{$scaffold}->{$geneID}->{"proteinId"} = $parsedAttributes{"proteinId"};
			}
			
			my @rnas = keys %{$genes{$scaffold}->{$geneID}->{"RNA"}};
			my $rnaCount = 1;
			my $n_rnas = @rnas;
			if ($n_rnas > 1)
			{
				print "warning: multiple RNAs/transcripts are listed in the original gff for this gene, to which should we add the CDS!?","\t",$geneName, " ", "@rnas","\n";
			}
			else
			{
				#print "warning: one RNA/transcript for this gene","\t",$geneName," ", "@rnas","\n";
				#continue
				my $rnaID = $rnas[0];
				
					
				#the gene is in the database (so just add the CDS)
				$tr_feature->{"start"} = $start;
				$tr_feature->{"stop"} = $stop;
				$tr_feature->{"phase"} = $phase;
				$tr_feature->{"score"} = $score;
				$tr_feature->{"ID"} = "CDS:" . $rnaID . ":" . $countTotalCDS;
				++$countTotalCDS;
				
				push(@{$genes{$scaffold}->{$geneID}->{"RNA"}->{$rnaID}->{"CDS"}}, $tr_feature);
				
				
				#let's add a translation_start and translation_stop to RNA
				if (! defined $genes{$scaffold}->{$geneID}->{"RNA"}->{$rnaID}->{"translation_start"})
				{
					$genes{$scaffold}->{$geneID}->{"RNA"}->{$rnaID}->{"translation_start"} = $start;
					$genes{$scaffold}->{$geneID}->{"RNA"}->{$rnaID}->{"translation_stop"} = $stop;
				}
				else
				{
					if ($start < $genes{$scaffold}->{$geneID}->{"RNA"}->{$rnaID}->{"translation_start"})
					{
						$genes{$scaffold}->{$geneID}->{"RNA"}->{$rnaID}->{"translation_start"} = $start;
					}
					elsif ($stop > $genes{$scaffold}->{$geneID}->{"RNA"}->{$rnaID}->{"translation_stop"})
					{
						$genes{$scaffold}->{$geneID}->{"RNA"}->{$rnaID}->{"translation_stop"} = $stop;
					}
				}
			}
		}
		else
		{
			#this gene doesn't exist : ERROR
			print "warning: no exon features for this gene; drop it","\t", $geneName,"\n";
		}
	}
}

close (GFF);

&add_utr;

#print Dumper (\%genes);	

sub add_utr
{
	#add UTR information per transcript based on current exon/CDS descriptions
	
	my @scaffolds = keys %genes; # getting an array of all scaffolds in %genes
	my @sorted_scaffolds = sort @scaffolds; # sort the scaffolds

	foreach my $scaffold (@sorted_scaffolds)
	{ 
		my @genesOnScaffold = keys %{$genes{$scaffold}};
		my @sortedGenesOnScaffold = sort {$genes{$scaffold}->{$a}->{"start"} <=> $genes{$scaffold}->{$b}->{"start"}} @genesOnScaffold;
		foreach my $gene (@sortedGenesOnScaffold)
		{	
			my @rnas = keys %{$genes{$scaffold}->{$gene}->{"RNA"}};
			my $strand = 
			my $rnaCount = 1;
			foreach my $rna (@rnas)
			{
				my @CDS = @{$genes{$scaffold}->{$gene}->{"RNA"}->{"$rna"}->{"CDS"}};
				my @exons = @{$genes{$scaffold}->{$gene}->{"RNA"}->{"$rna"}->{"exon"}};
				my @sorted_CDS = sort @CDS;
				my @sorted_exons = sort @exons;
				
				my (@CDSs, @exonss);				
				foreach my $CDS (@sorted_CDS)
				{
					my $CDS_start = $CDS->{"start"};
					my $CDS_stop = $CDS->{"stop"};
					push (@CDSs, [$CDS_start, $CDS_stop]);			
				}
				foreach my $exon (@sorted_exons)
				{
					my $exon_start = $exon->{"start"};
					my $exon_stop = $exon->{"stop"};
					push (@exonss, [$exon_start, $exon_stop]);			
				}
				
				#now check for each exon_coordset, whether it overlaps a CDS coordset from @CDSs
				foreach my $exon_coordset (@exonss)
				{
					
					my $exon_CDS_overlap_switch = 0; #set a switch so to check in the end whether a certain exon has CDS overlap (at all); could be full UTRs!
					my $exon_token = join ("_", @$exon_coordset);
					my ($exon_lend, $exon_rend) = split (/_/, $exon_token);
					
					my $utr_type; my $tr_feature;
					
					foreach my $cds_coordset (@CDSs)
					{
						my $cds_token = join ("_", @$cds_coordset);
						my ($cds_lend, $cds_rend) = split (/_/, $cds_token);
						if ( ($cds_lend > $exon_lend) && ($cds_lend < $exon_rend) ) #left partial overlap
						{
							$exon_CDS_overlap_switch = 1; #partial overlap (the rest = UTR)
							
							#utr from exon_lend to $cds_lend-1
							
							#now check whether this is a 5'UTR or a 3'UTR
							#if $strand eq "+" then this is upstream (does 5'UTR)
							#if $strand eq "-" then this is downstream (does 3'UTR)
							
							if ($genes{$scaffold}->{$gene}->{"orientation"} eq "+")
							{
								$utr_type = "five_prime_UTR";
								++$countTotal5UTR;
							}
							else
							{
								$utr_type = "three_prime_UTR";
								++$countTotal3UTR;
							}
							
							#print join ("\t",$rna,$exon_lend,$exon_rend,$cds_lend,$cds_rend,$exon_lend,$cds_lend-1),"\n";
							$tr_feature->{"start"} = $exon_lend;
							$tr_feature->{"stop"} = $cds_lend-1;
							$tr_feature->{"ID"} = "UTR:" . $rna . ":" . $countTotalUTR;
							$tr_feature->{"type"} = $utr_type;
							++$countTotalUTR;
							push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"UTR"}}, $tr_feature);				
						}
						elsif ( ($cds_rend < $exon_rend) && ($cds_rend > $exon_lend) ) #right partial overlap
						{
							$exon_CDS_overlap_switch = 1; #partial overlap (the rest = UTR)
							
							#utr from cds_rend+1 to exon_rend
							
							#now check whether this is a 5'UTR or a 3'UTR
							#if $strand eq "+" then this is downstream (does 3'UTR)
							#if $strand eq "-" then this is upstream (does 5'UTR)
							
							if ($genes{$scaffold}->{$gene}->{"orientation"} eq "+")
							{
								$utr_type = "three_prime_UTR";
								++$countTotal3UTR;
							}
							else
							{
								$utr_type = "five_prime_UTR";
								++$countTotal5UTR;
							}
							
							#print join ("\t",$rna,$exon_lend,$exon_rend,$cds_lend,$cds_rend,$cds_rend+1,$exon_rend),"\n";
							$tr_feature->{"start"} = $cds_rend + 1;
							$tr_feature->{"stop"} = $exon_rend;
							$tr_feature->{"ID"} = "UTR:" . $rna . ":" . $countTotalUTR;
							$tr_feature->{"type"} = $utr_type;
							++$countTotalUTR;
							push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"UTR"}}, $tr_feature);
						}
						elsif ( ($cds_lend == $exon_lend) && ($cds_rend == $exon_rend) ) #full overlap exon - cds
						{
							#no utr, fully overlapping CDS and exon features
							$exon_CDS_overlap_switch = 1;
						}
					}
					
					if ($exon_CDS_overlap_switch == 0) #no CDS overlap detected for this exon coordset >> full UTR
					{
						
						my $test_orient = $genes{$scaffold}->{$gene}->{"orientation"};
						my $test_tstart = $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"translation_start"};
						my $test_tstop = $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"translation_stop"};
						
						
						if ( ( ($test_orient eq "+") && ($exon_lend < $test_tstart) ) || ( ($test_orient eq "-") && ($exon_lend > $test_tstart) ) )
						{
							$utr_type = "five_prime_UTR";
						}
						elsif ( ( ($test_orient eq "+") && ($exon_lend > $test_tstart) ) || ( ($test_orient eq "-") && ($exon_lend < $test_tstart) ) )
						{
							$utr_type = "three_prime_UTR";
						}
						
						
						$tr_feature->{"start"} = $exon_lend;
						$tr_feature->{"stop"} = $exon_rend;
						$tr_feature->{"ID"} = "UTR:" . $rna . ":" . $countTotalUTR;
						$tr_feature->{"type"} = $utr_type;
						++$countTotalUTR;
						push(@{$genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"UTR"}}, $tr_feature);
					}
					else
					{
						#this exon had CDS overlap, do nothing here
					}
				}
			}
		}
	}				
}



open (GFFOUT, ">$gff3_out");
print GFFOUT "##gff-version 3","\n";
print GFFOUT "##created using JGI_TO_GFF3_1.pl on ",scalar localtime(),"\n";
print GFFOUT "##updated $gff","\n","###","\n";

my @scaffolds = keys %genes; # getting an array of all scaffolds in %genes
my @sorted_scaffolds = sort @scaffolds; # sort the scaffolds

foreach my $scaffold (@sorted_scaffolds)
{ 
	my @genesOnScaffold = keys %{$genes{$scaffold}};
	my @sortedGenesOnScaffold = sort {$genes{$scaffold}->{$a}->{"start"} <=> $genes{$scaffold}->{$b}->{"start"}} @genesOnScaffold;

	foreach my $gene (@sortedGenesOnScaffold)
	{
		&print_gene($gene, $scaffold);
	}
}

sub print_gene  	# subroutine to print gene structure from gene tree in hash structure (gene > rna > exon > cds || utr) 
{
	my ($gene, $scaffold) = @_;
	my $orient = $genes{$scaffold}->{$gene}->{"orientation"};
	my $source = $genes{$scaffold}->{$gene}->{"new_source"};
	
	#print gene line (including name, description and notes if they exist)
	print GFFOUT join ("\t", $scaffold, $source,"gene",$genes{$scaffold}->{$gene}->{"start"},$genes{$scaffold}->{$gene}->{"stop"},".",$orient,".","");
	print GFFOUT "ID=",$genes{$scaffold}->{$gene}->{"new_ID"},";Name=",$genes{$scaffold}->{$gene}->{"name"},"\n";
	
	#now for each transcript/mRNA, print mRNA+exons+CDS+UTR
	my @rnas = keys %{$genes{$scaffold}->{$gene}->{"RNA"}};
	foreach my $rna (@rnas)
	{
		my $rnaID = $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"new_name"};
		print GFFOUT $scaffold,"\t",$source,"\t","mRNA","\t";
		print GFFOUT $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"start"}, "\t", $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"stop"},"\t",".","\t";
		print GFFOUT $orient,"\t",".","\t","ID=",$rnaID,";Parent=",$genes{$scaffold}->{$gene}->{"new_ID"},";jgiTranscriptId=",$genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"ID"}, ";jgiProteinId=",$genes{$scaffold}->{$gene}->{"proteinId"},"\n";
		
		my @exon = @{$genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"exon"}};
		my @CDS = @{$genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"CDS"}};
		
		my @sorted_UTR; my @UTR;
		
		if (defined $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"UTR"})
		{
			@UTR = @{$genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"UTR"}};
		}
		
		my @sorted_exon = sort @exon;
		my @sorted_CDS = sort @CDS;
		@sorted_UTR = sort @UTR;
		
		foreach my $exon (@sorted_exon)
		{
			my $exon_start = $exon->{"start"};
			my $exon_stop = $exon->{"stop"};
			my $exon_id = $exon->{"ID"};
			print GFFOUT join ("\t", $scaffold, $source, "exon", $exon_start, $exon_stop, ".", $orient, ".", "ID=");
			print GFFOUT $exon_id,";Parent=",$rnaID,"\n";
		}
		foreach my $CDS (@sorted_CDS)
		{
			my $CDS_start = $CDS->{"start"};
			my $CDS_stop = $CDS->{"stop"};
			my $CDS_phase = $CDS->{"phase"};
			my $CDS_id = $CDS->{"ID"};
			print GFFOUT join ("\t", $scaffold, $source, "CDS", $CDS_start, $CDS_stop, ".", $orient, $CDS_phase, "ID=");
			print GFFOUT $CDS_id,";Parent=",$rnaID,"\n";
		}
		foreach my $UTR (@sorted_UTR)
		{
			my $UTR_start = $UTR->{"start"};		
			my $UTR_stop = $UTR->{"stop"};
			my $UTR_id = $UTR->{"ID"};
			my $UTR_type = $UTR->{"type"};
			print GFFOUT join ("\t", $scaffold, $source, $UTR_type, $UTR_start, $UTR_stop, ".", $orient, ".", "ID=");
			print GFFOUT $UTR_id, ";Parent=",$rnaID,"\n";
		}
	}
	print GFFOUT "###","\n";
}
