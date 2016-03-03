## analyseOrthoMCL.pl

#!/usr/bin/perl

####
#Perl script by Ronnie de Jonge, ronnie.dejonge@gmail.com OR ronnie.dejonge@psb.ugent.be
#Created: October 20, 2014
#Last modified: October 20, 2014
####

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use experimental 'smartmatch';
use List::MoreUtils qw(uniq);
use Bio::SeqIO;
use Bio::Seq;
use Bio::Index::Fasta;

my $usage = "\n\nusage:	$0 
				--clusters <orthoMCL cluster file; mclGroups.txt>
				--goodproteins <goodProteins.fasta>
				--poorproteins <poorProteins.fasta>
				--mclswitch [0/1] (optional, default=0)
\n\n";


if (@ARGV < 3)
{
	print $usage;
	exit;
}

my ($clusters, $good, $poor);
my $mcl = 0;
my %group;
my %proteinToGroup;


GetOptions("clusters=s" => \$clusters, 
	   "goodproteins=s" => \$good,
	   "poorproteins=s" => \$poor,
	   "mclswitch=i" 	=> \$mcl,
);


#what do I need (for each input species)?
#1. list of unique proteins (including proteins from poorProteins)
#2. list of core orthologs (1:1:1:1 orthologs, present in single copy in each species)
#3. list of proteins in paralog clusters (no orthologs, one species, multi-protein)
#4. list of proteins in paralog/ortholog clusters (no single copy core, multi-species, multi-protein)

#get protein lists
my @goodProteins;
my @poorProteins;
my @allProteins;

my $goodIn = Bio::SeqIO->new(-file => "$good", -format => "fasta"); 
while (my $goodSeq = $goodIn->next_seq())
{
	push (@goodProteins, $goodSeq->display_id);
	push (@allProteins, $goodSeq->display_id);
}
my $poorIn = Bio::SeqIO->new(-file => "$poor", -format => "fasta"); 
while (my $poorSeq = $poorIn->next_seq())
{
	push (@poorProteins, $poorSeq->display_id);
	push (@allProteins, $poorSeq->display_id);
}

my $goodProteinNo = scalar @goodProteins;
my $poorProteinNo = scalar @poorProteins;
my $allProteinNo = scalar @allProteins;

print "\nNumber of proteins that went into orthoMCL (good and poor proteins)\n";
print "good: ", $goodProteinNo, " poor: ", $poorProteinNo, " all: ", $allProteinNo, "\n";

#parse orthoMCL groups file

my @genomeids; my $count = 1;
open (CLUSTER, $clusters) or die $usage;
while (my $line = <CLUSTER>)
{
	my @current_genomeids;
	chomp $line;
	my ($groupid, $proteinlist); my @group_proteins;
	if ($mcl == 0)
	{ #mclGroups.txt (final output from orthoMCL pipeline) with groupid
		($groupid, $proteinlist) = split (/: /, $line);
		@group_proteins = split (/ /, $proteinlist);
	}
	elsif ($mcl == 1)
	{ #if no groupid has been defined, we do it now here (mcl output case)
		$proteinlist = $line;
		@group_proteins = split (/\t/, $proteinlist);
		$groupid = "MCL" . $count;
		++$count;
	}	
	my $num_prots = scalar @group_proteins;
	$group{$groupid}->{"ID"} = $groupid;
	$group{$groupid}->{"n_proteins"} = $num_prots;
	
	foreach my $protein (@group_proteins)
	{
		my ($genomeid, $proteinid) = split (/\|/, $protein);
		push (@genomeids, $genomeid);
		push (@current_genomeids, $genomeid);
		my $protein_member;
		$protein_member->{"ID"} = $proteinid;
		$protein_member->{"genome"} = $genomeid;
		push (@{$group{$groupid}->{"member"}}, $protein_member);
		$proteinToGroup{$proteinid}->{"group"} = $groupid;
	}
	my @unique_current_genomeids = uniq @current_genomeids;
	$group{$groupid}->{"uniqId"} = scalar @unique_current_genomeids;
}

my @unique_genomeids = uniq @genomeids;
my @sorted_unique_genomeids = sort @unique_genomeids;
my $no_unique_genomeids = scalar @unique_genomeids;

print "Uniq genomeIDs: ", "@sorted_unique_genomeids", "\nNo. of unique genomeIDs: ", $no_unique_genomeids,"\n";

my %prots;
foreach my $unique_genomeid (@sorted_unique_genomeids)
{
	$prots{$unique_genomeid}->{"core"} = 0;
	$prots{$unique_genomeid}->{"unique"} = 0;
	$prots{$unique_genomeid}->{"shared"} = 0;
	$prots{$unique_genomeid}->{"only_paralogs"} = 0;
}


my @clusters = keys %group;
my $no_clusters = scalar @clusters;

print "No. of clusters: ", $no_clusters,"\n";

my (@paraClusters, @coreClusters);
my $paraClusterProteins = 0;

foreach my $cluster (@clusters)
{
	if ($group{$cluster}->{"uniqId"} == $no_unique_genomeids)
	{ #cluster proteins present in all species
		$group{$cluster}->{"allId"} = 1;
		$group{$cluster}->{"onlyParalogs"} = 0;
	}
	elsif ($group{$cluster}->{"uniqId"} == 1)
	{ #cluster proteins present in only one species (only paralogs)
		$group{$cluster}->{"allId"} = 0;
		$group{$cluster}->{"onlyParalogs"} = 1;
		push (@paraClusters, $cluster);
		$paraClusterProteins = $paraClusterProteins + $group{$cluster}->{"n_proteins"};
	}
	else
	{ #cluster proteins present in more than one species, but not all (>1 , <total)
		$group{$cluster}->{"allId"} = 0;
		$group{$cluster}->{"onlyParalogs"} = 0;
	}

	if ( ($group{$cluster}->{"n_proteins"} == $no_unique_genomeids) && ($group{$cluster}->{"allId"} == 1) ) 
	{ #clusters with all Ids as well as only one copy/species (single copy core clusters)
		$group{$cluster}->{"core"} = 1;
		push (@coreClusters, $cluster);
	}
	else
	{ #all other clusters
		$group{$cluster}->{"core"} = 0;
	}
}

print "Core clusters: ", scalar @coreClusters, " Para clusters: ", scalar @paraClusters, " Para cluster proteins: ", $paraClusterProteins, "\n";

my @sortgoodProteins = sort @goodProteins;

open (TAGS, ">$clusters.TAGs") or die $usage;

foreach my $protein (@sortgoodProteins)
{
	my ($genomeid, $proteinid) = split (/\|/, $protein);
	print TAGS $genomeid, "\t", $protein, "\t";
	my $set;
	if (defined $proteinToGroup{$proteinid}->{"group"})
	{
		
		my $current_cluster = $proteinToGroup{$proteinid}->{"group"};
		if ($group{$current_cluster}->{"core"} == 1)
		{ #core cluster
			my $prevCount = $prots{$genomeid}->{"core"};
			my $newCount = $prevCount + 1;
			$prots{$genomeid}->{"core"} = $newCount;
			$set = "core";
		}
		elsif ($group{$current_cluster}->{"onlyParalogs"} == 1)
		{ #non-core , only-paralog clusters
			my $prevCount = $prots{$genomeid}->{"only_paralogs"};
			my $newCount = $prevCount + 1;
			$prots{$genomeid}->{"only_paralogs"} = $newCount;
			$set = "only_paralogs";	
		}
		else
		{ #non-core , shared cluster
			my $prevCount = $prots{$genomeid}->{"shared"};
			my $newCount = $prevCount + 1;
			$prots{$genomeid}->{"shared"} = $newCount;
			$set = "shared";				
		}
	}
	else
	{ #protein not in an cluster; unique (add id to array and count to %prots
		my $prevCount = $prots{$genomeid}->{"unique"};
		my $newCount = $prevCount + 1;
		$prots{$genomeid}->{"unique"} = $newCount;
		$set = "unique";
	}
	print TAGS $set,"\n";
}

close (TAGS);
print Dumper (\%prots);





