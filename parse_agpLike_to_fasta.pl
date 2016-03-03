## parse_agpLike_to_fasta.pl

#!/usr/bin/perl -w

#script to parse super-basic AGP-like files to chromosome FASTA files
#extensive usage of BioPerl

use warnings;
use strict;

use Bio::SeqIO;
use Bio::Seq;
use Bio::Index::Fasta;

#Chr1	scaffold118	-
#Chr1	gap	13000

my $agp = $ARGV[0];
my $contigs = $ARGV[1];
my $out = $ARGV[2];
my $prev_chr = 0;


unless (-e "$contigs.inx")
{
        my $fasta_index = Bio::Index::Fasta->new(
                                                -filename => "$contigs.inx",
                                                -write_flag => 1 );
        $fasta_index->make_index($contigs);
}
my $index = Bio::Index::Fasta->new("$contigs.inx");

open (OUT, ">$out");

open (FILE, $agp)
   or die "impossible to open file $agp for reading !";

while (my $line = <FILE>) 
{
	chomp $line;
	next if ($line =~ /^#/);
	
	my ($chr, $type, $info) = split("\t", $line);
	
	if ($chr ne $prev_chr)
	{
		print "\n";
		print join ("\t", $chr, $type, $info),"\n";		
		print OUT "\n",">",$chr,"\n";
		if ($type eq "gap") 
		{
			print "writing a gap of length $info","\n";
			my $string;
			my $N = "N";
			for (my $i = 0; $i < $info; $i++)
			{
				$string = $string . $N;
			}
			$string =~ s/(.{1,60})/$1\n/gs;
			print OUT $string;
		}
		
		elsif ($type =~ /scaffold/)
		{			
			my $seq = $index->fetch($type);
			print $seq->display_id,"\n";
			print "new or first chromosome, printing first scaffold $type","\n";
			if ($info eq "+")
			{
				print "positive strand sequence $type","\n";
				my $print_seq = $seq->seq;
				$print_seq =~ s/(.{1,60})/$1\n/gs;
				print OUT $print_seq;
			}
			elsif ($info eq "-")
			{
				print "reverse complementing sequence $type","\n";
				my $seq_rev = $seq->revcom;
				my $print_seq = $seq_rev->seq;
				$print_seq =~ s/(.{1,60})/$1\n/gs;
				print OUT $print_seq;
			}
			else
			{
				print "problem?","\n";
			}

		}
	$prev_chr = $chr;
	}
	else
	{
		print join ("\t", $chr, $type, $info),"\n";		
		if ($type eq "gap")
		{
			print "writing a gap of length $info","\n";			
			my $string;
			my $N = "N";
			for (my $i = 0; $i < $info; $i++)
			{
				$string = $string . $N;
			}
			$string =~ s/(.{1,60})/$1\n/gs;
			print OUT $string;
		}
	
		elsif ($type =~ /scaffold/)
		{
			my $seq = $index->fetch($type);
			print $seq->display_id,"\n";
			print "same chromosome, printing next scaffold $type","\n";
			if ($info eq "+")
			{
				print "positive strand sequence $type","\n";
				my $print_seq = $seq->seq;
                        	$print_seq =~ s/(.{1,60})/$1\n/gs;
                        	print OUT $print_seq;
			}
			elsif ($info eq "-")
			{
				print "reverse complementing sequence $type","\n";
				my $seq_rev = $seq->revcom;
				my $print_seq = $seq_rev->seq;
                                $print_seq =~ s/(.{1,60})/$1\n/gs;
                        	print OUT $print_seq;
			}
			else
			{
				print "problem?","\n";
			}
		}
	}
}

close FILE;
close OUT;
exit;
