## get_seq_by_id.pl

#!/usr/bin/perl

#get sequence by ID 
#simple and slow ... (no indexing)

use strict;
use warnings;

use Bio::SeqIO;

my $db = $ARGV[0];
my $seq_idq = $ARGV[1];
my $region_switch = $ARGV[2];
my $region_start = $ARGV[3];
my $region_stop = $ARGV[4];

my $seq_db = Bio::SeqIO->new(-file=>$db);

while (my $seqio = $seq_db->next_seq()) {
  my $seq_id = $seqio->display_id();
  if ($seq_idq eq $seq_id) {
  	if ($region_switch == 1) {
  		my $subseq = $seqio->subseq($region_start, $region_stop);
  		print ">",$seq_id,"\t","$region_start-$region_stop","\n",$subseq,"\n";
  	}
	elsif ($region_switch == 0) {
  		print ">",$seq_id,"\n",$seqio->seq,"\n";
  	}
  }
  else {
  }


}

exit;
