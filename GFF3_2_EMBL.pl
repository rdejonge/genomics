## GFF3_2_EMBL.pl

#!/usr/bin/perl

# perl script to convert gff3+fasta into EMBL format (for antismash analyses)
# Author: Ronnie de Jonge, rojon@psb.ugent.be
# Last update on August 29, 2014

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq;
use File::Basename;
use Data::Dumper;

#Input arguments:
my $gff3 = $ARGV[0];
my $fasta = $ARGV[1];
my $organism = $ARGV[2];

my $base = basename($gff3, ".gff3");

#insert the gene annotations from GFF3
my %genes; my %rnaToGene; my %new_genes; my %new_rnaToGene;
open (GFF3, $gff3) or die "can't open this file $gff3";
my $currentGene = "";
while (my $line = <GFF3>)
{
	chomp $line;
	next if ($line =~ /^#/); #skip commented lines (STARTING WITH HASHTAG !!!, some @#@$&Q#$ JGI ID's have ### too)
	
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
		
		$currentGene = $geneID;

	} elsif ($type eq "mRNA" || $type eq "transcript")
	{
		#new RNA gene, add data to %genes and %rnaToGene
		my $ID = $parsedAttributes{"ID"};
		my $Parent = $currentGene;
		$rnaToGene{$ID} = $Parent;
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
		my $ParentRNA = $parsedAttributes{"Parent"};
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

my %scaffolds;
#load fasta db
my $seqin = Bio::SeqIO->new(-file=>$fasta);

while (my $seqio = $seqin->next_seq() ) 
{
	my $scaffold_id = $seqio->display_id;
	my $scaffold_length = $seqio->length();
	my $scaffold_seq = $seqio->seq;
	print $seqio->display_id, "\t",$seqio->length(),"\n";
	$scaffolds{$scaffold_id}->{"length"} = $scaffold_length;
	$scaffolds{$scaffold_id}->{"seq"} = $scaffold_seq;
}


#write EMBL files, one for each scaffold/chromosome 
##FIX: For now only works for scaffolds that have predicted genes



my @scaffolds = keys %genes;
my @sorted_scaffolds = sort @scaffolds;

foreach my $scaffold (@sorted_scaffolds)
{
	open (OUT, ">$base.$scaffold.embl") or die "$!";
	my $sLength = $scaffolds{$scaffold}->{"length"};
	print OUT "ID   $scaffold    standard; DNA;  HTG; $sLength BP.\n";
	print OUT "XX\n";
	print OUT "FH   Key             Location\/Qualifiers\n";
	print OUT "FH\n";
	print OUT "FT   source          1..";
	print OUT $sLength,"\n";
	print OUT "FT                   \/mol_type\=\"genomic DNA\"","\n";
	print OUT "FT                   \/organism\=\"$organism\"","\n";

	my @genesOnScaffold = keys %{$genes{$scaffold}};
	my @sortedGenesOnScaffold = sort {$genes{$scaffold}->{$a}->{"start"} <=> $genes{$scaffold}->{$b}->{"start"}} @genesOnScaffold;
	
	foreach my $gene (@sortedGenesOnScaffold)
	{	
		my $strand = $genes{$scaffold}->{$gene}->{"orientation"};
		if ($strand eq "+"){ $strand = "1"; }
		elsif ($strand eq "-"){ $strand = "-1"; }
		print OUT "FT   gene            ";
		if ($strand eq "1")
		{
			print OUT "\<",$genes{$scaffold}->{$gene}->{"start"},"..\>",$genes{$scaffold}->{$gene}->{"stop"},"\n";
		}
		elsif ($strand eq  "-1")
		{	
			print OUT "complement\(\<",$genes{$scaffold}->{$gene}->{"start"},"..\>",$genes{$scaffold}->{$gene}->{"stop"},"\)","\n";
		}
		print OUT "FT                   /note=\"";
		print OUT $genes{$scaffold}->{$gene}->{"ID"},"\"","\n";

		print OUT "FT                   /gene=\"";
		print OUT $genes{$scaffold}->{$gene}->{"ID"},"\"","\n";

		
		my @rnas = keys %{$genes{$scaffold}->{$gene}->{"RNA"}};
		foreach my $rna (@rnas)
		{
			my $rna_start = $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"start"};
			my $rna_stop = $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"stop"};
			my @CDS = @{$genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"CDS"}};
			my $protein_id = $genes{$scaffold}->{$gene}->{"RNA"}->{$rna}->{"ID"};

			my %cdss;
			my @sorted_CDS;
			#get start and stop coordinates from the extremeties of the CDS coordinates
			foreach my $cds (@CDS)
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

			my @CDSs;	
			foreach my $CDS (@sorted_CDS)
			{
				my $CDS_start = $cdss{$CDS}->{"start"};
				my $CDS_stop = $cdss{$CDS}->{"stop"};

				push (@CDSs, [$CDS_start, $CDS_stop]);	
			}
			
			my @tokens;
			foreach my $CDS_coordset (@CDSs)
			{
				my $CDS_token = join ("..", @$CDS_coordset);
				push (@tokens, $CDS_token);		
			}
			
			my $mRNA_string = join ',', @tokens;

			if(scalar @CDS > 1)
			{
				$mRNA_string = 'join('.$mRNA_string.')';
			}

			if($strand eq '-1'){
				$mRNA_string = 'complement('.$mRNA_string.')';
			}
			
			#replace the last occurrence of '..' by '..>'
			$mRNA_string =~ s/\.\.(?!.*\.\..*)/\.\.>/;
			
			#add an opening < before the first number, but after the 'join(' statement (if it exists)
			if ( ($mRNA_string =~ /complement\(join/) || ($mRNA_string =~ /join/) )
			{
				$mRNA_string =~ s/join\(/join\(</;
			}
			elsif ($mRNA_string =~ /complement\(\d+/)
			{
				$mRNA_string =~ s/complement\(/complement\(</;
			}
			else
			{
				$mRNA_string = "<" . $mRNA_string;
			}
			#FIXME: now hardcoded for maximum number of 4*59 characters...
			
			#if the mRNA_string is larger than 59 characters (including complement( , complement(join( or join( 
			#split to multiple lines
			
			#complement(join(<52083..52836,52899..53022,53073..53199,53263..>53418))
			
			my ($mRNA_string1, $mRNA_string2, $mRNA_string3, $mRNA_string4) = ("","","","");
			
			my @mRNA_strings = split(',',$mRNA_string);
			my $n_strings = $#mRNA_strings + 1;
			my $flag = 1;
			my $count;
			my $characters = 59;
			$mRNA_string1 = $mRNA_strings[0];
			
			if ($n_strings > 1)
			{			
				for ($count = 1 ; $count < $n_strings ; $count++)
				{ 
					if ( (( length($mRNA_string1) + length($mRNA_strings[$count]) + 1 ) < $characters ) && ($flag == 1))
					{
						$mRNA_string1 = $mRNA_string1 . "," . $mRNA_strings[$count];
						$flag = 1;
					}
					else
					{
						if ( ( ( length($mRNA_string2) + length($mRNA_strings[$count]) + 1 ) < ($characters*2) ) && ($flag == 1 || $flag == 2) )
						{
							if ($flag == 1)
							{
								$mRNA_string1 .= ",";
								$mRNA_string2 = $mRNA_strings[$count];
								$flag = 2;
							}
							elsif ($flag == 2)
							{
								$mRNA_string2 = $mRNA_string2 . "," . $mRNA_strings[$count];
								$flag = 2;
							}
						}
						else
						{
							if ( (( length($mRNA_string3) + length($mRNA_strings[$count]) + 1 ) < ($characters*3)) && ($flag == 2 || $flag == 3) )
							{
								if ($flag == 2)
								{
									$mRNA_string2 .= ",";
									$mRNA_string3 = $mRNA_strings[$count];
									$flag = 3;
								}
								else
								{
									$mRNA_string3 = $mRNA_string3 . "," . $mRNA_strings[$count];
									$flag = 3;
								}	
							}
							else
							{
								if ( ( length($mRNA_string4) + length($mRNA_strings[$count]) + 1 ) < ($characters*4) )
								{
									if ($flag == 3)
									{
										$mRNA_string3 .= ",";
										$mRNA_string4 = $mRNA_strings[$count];
										$flag = 4;
									}
									else
									{									
										$mRNA_string4 = $mRNA_string4 . "," . $mRNA_strings[$count];
										$flag = 4;
									}
								}
							}
						}
					}
				}
			}
			else
			{
				$flag = 1;
			}
			
			my ($CDS_string1, $CDS_string2, $CDS_string3, $CDS_string4) = ($mRNA_string1, $mRNA_string2, $mRNA_string3, $mRNA_string4);
			
			$CDS_string1 =~ s/<//; $CDS_string1 =~ s/>//;
			$CDS_string2 =~ s/<//; $CDS_string2 =~ s/>//;
			$CDS_string3 =~ s/<//; $CDS_string3 =~ s/>//;
			$CDS_string4 =~ s/<//; $CDS_string4 =~ s/>//;
			
			if ($flag == 1)
			{
				print OUT "FT   mRNA            $mRNA_string1\n";
				print OUT "FT                   \/gene\=\"$gene\"\n";
				print OUT "FT   CDS             $CDS_string1\n";
			}
			
			if ($flag == 2)
			{
				print OUT "FT   mRNA            $mRNA_string1\n";
				print OUT "FT                   $mRNA_string2\n";
				print OUT "FT                   \/gene\=\"$gene\"\n";
				print OUT "FT   CDS             $CDS_string1\n";
				print OUT "FT                   $CDS_string2\n";
			}
			elsif ($flag == 3)
			{
				print OUT "FT   mRNA            $mRNA_string1\n";
				print OUT "FT                   $mRNA_string2\n";
				print OUT "FT                   $mRNA_string3\n";
				print OUT "FT                   \/gene\=\"$gene\"\n";
				print OUT "FT   CDS             $CDS_string1\n";
				print OUT "FT                   $CDS_string2\n";
				print OUT "FT                   $CDS_string3\n";
			}
			elsif ($flag == 4)
			{
				print OUT "FT   mRNA            $mRNA_string1\n";
				print OUT "FT                   $mRNA_string2\n";
				print OUT "FT                   $mRNA_string3\n";
				print OUT "FT                   $mRNA_string4\n";
				print OUT "FT                   \/gene\=\"$gene\"\n";
				print OUT "FT   CDS             $CDS_string1\n";
				print OUT "FT                   $CDS_string2\n";
				print OUT "FT                   $CDS_string3\n";
				print OUT "FT                   $CDS_string4\n";
			}
			print OUT "FT                   \/gene\=\"$gene\"\n";
			print OUT "FT                   \/locus_tag\=\"$protein_id\"\n";
			print OUT "FT                   \/codon_start\=1\n";
			print OUT "FT                   \/protein_id\=\"$protein_id\"\n";
		}
	}
	
	print OUT "XX\n";
	print OUT "SQ   Sequence ", $sLength, " BP\;\n";
	
	my $seq_lines = roundup($sLength / 60);
	my $embl_format = '';
	my $end_flag = 0;
	foreach my $line (1..$seq_lines)
	{
		$embl_format .= '     ';

		foreach my $col(1..6)
		{
			next if $end_flag;

			my $start = (($line - 1) * 60 ) + (($col - 1) * 10);

			my $seq_length = 10;

			if(($start + $seq_length) > $sLength)
			{
				$seq_length = $sLength - $start;
				$end_flag = 1;
			}

			$embl_format .= substr($scaffolds{$scaffold}->{"seq"}, $start, $seq_length);
			$embl_format .= ' ';

		}
	
		my $bp_label;
		my $spacing;

		if($end_flag)
		{
			$bp_label = $sLength;
			my $label_length = length $bp_label;
			my $left_over = ($sLength % 60);
			my $extra_gaps = int($left_over / 10);
			$spacing = 79 - 5 - $extra_gaps - $left_over - $label_length;

		}
		else
		{
			$bp_label = $line * 60;
			my $label_length = length $bp_label;
			$spacing = 79 - 10 - 60 - $label_length;
		}

		$embl_format .= ' ' x $spacing ."$bp_label\n";

		my $embl_length = length $embl_format;
	}
	$embl_format .= "\/\/\n";

	print OUT "$embl_format";
	print OUT "\n";
	close OUT;
}


sub roundup 
{
	my $n = shift;
	return(($n == int($n)) ? $n : int($n + 1))
}



