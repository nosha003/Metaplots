#! /usr/bin/perl -w
# use bedtools intersect from bins peaks and gene gff output... pull out data to use for metaplot
# TEfamily_macs.pl
# 2_Jan_2017
# Jaclyn_Noshay

use warnings;
use strict;
use Getopt::Std;

#set soft coded files
my $usage = "\n0 -i in -o out\n";
our ($opt_i, $opt_o, $opt_h);
getopts("i:o:h") || die "$usage";

#check that all files are defined
if ( (!(defined $opt_i)) || (!(defined $opt_o)) || (defined $opt_h) ) {
  print "$usage";
}

#read in bedtools intersect file 
open (my $in_fh, '<', $opt_i) || die;
open (my $out_fh, '>', $opt_o) || die;

print $out_fh "chr\tbinstart\tbinstop\tbinid\tchr\tgene\tgenestart\tgenestop\tgenesize\tstrand\tdistance\tstrand_distance\trelative_distance\treal_distance\tcount\n";
##print $out_fh "chr\tbinstart\tbinstop\tbinid\tchr\tTE\tTEstart\tTEstop\tTEsize\tstrand\tdistance\tstrand_distance\trelative_distance\treal_distance\tcount\n";

my $strand_dist;
my $relative_dist;
my $real_dist;
my $genesize;
my $header = <$in_fh>;
while (my $line = <$in_fh>) {
  chomp $line;
  # 1	0	100	1	1	LTRharvest	LTR_retrotransposon	30293	38691	.	-	.	ID=RLC00004Zm00004b00001	30193
  #my ($chr, $start, $end, $count, $TEchr, undef, $TE, $TEstart, $TEend, undef, $strand, undef, $TEID, $dist_bin) = split ("\t", $line);

  #my ($chr, $start, $end, $count, $genechr, undef, $gene, $genestart, $geneend, undef, $strand, undef, $geneID, $dist_bin) = split ("\t", $line);
  # 1	0	100	1	1	ENSEMBLPEP	gene	12601	18309	.	+	.	ID=Zm00004b000001;Name=Zm00004b000001;biotype=protein_coding	12501
  
  #1	0	100	1	1	ENSEMBLPEP	gene	12601	18309	.	+	.	ID=Zm00004b000001;Name=Zm00004b000001;biotype=protein_coding	-12501
  #1	0	100	1	1	ENSEMBLPEP	gene	12601	18309	.	+	.	ID=Zm00004b000001;Name=Zm00004b000001;biotype=protein_coding	-12501

  my ($chr, $start, $end, $count, $genechr, undef, $gene, $genestart, $geneend, undef, $strand, undef, $geneID, $dist_bin) = split ("\t", $line);

  #my ($chr, $start, $end, undef, undef, undef, $count, $genechr, $undef, $gene, $genestart, $geneend, undef, $strand, undef, $geneID, $dist_bin) = split ("\t", $line);
  #chr09	0	100	6	25	29	0.862068965517241	chr09	phytozomev12	gene	200205	207396	.	-	.	ID=Zm00008a033253.v1.1;Name=Zm00008a033253	20
0105

  my ($gene_name, undef) = split (";", $geneID);
  $gene_name =~ s/ID=gene://;
  #my ($TE_name, undef) = split("Zm", $TEID);
  #$TE_name =~ s/ID=//;

  my $dist = $dist_bin * -1;

  if (!($strand eq "+")) {
    $strand_dist = $dist * -1;
  }
  else {
    $strand_dist = $dist;
  }
  
  if ($strand_dist > 0) {
    $relative_dist = $strand_dist + 1000;
  }
  else {
    $relative_dist = $strand_dist;
  }

  $genesize = abs($geneend - $genestart);
  #$genesize = abs($TEend - $TEstart);

  if ($strand_dist == 0) {
    my $diff = $end - $genestart;
    #my $diff = $end - $TEstart;
    my $decimal = $diff / $genesize;
    $relative_dist = $decimal * 1000;
    $real_dist = $diff;
  }
  else {
    $real_dist = $strand_dist;
  }
  my $binid = $chr . "_" . $start . "_" . $end;

  print $out_fh "$chr\t$start\t$end\t$binid\t$genechr\t$gene_name\t$genestart\t$geneend\t$genesize\t$strand\t$dist\t$strand_dist\t$relative_dist\t$real_dist\t$count\n";
  #print $out_fh "$chr\t$start\t$end\t$binid\t$TEchr\t$TE_name\t$TEstart\t$TEend\t$genesize\t$strand\t$dist\t$strand_dist\t$relative_dist\t$real_dist\t$count\n";
}

close $in_fh;
close $out_fh;
exit;

