#!bin/perl

use strict;
use warnings;

open (hand1, "dmel-all-transcript-r6.07.fasta") or die $!;
my %mrna; my $id;

while (<hand1>)    {
  chomp;
  if (/^>(FBtr[0-9]*)\s/)    {
     $id = $1;   next }
  $mrna{$id} .= $_;   }

close hand1;

open (hand2, "dmel-all-five_prime_UTR-r6.07.fasta") or die $!;
my %utr5; 

while (<hand2>)     {
  chomp;
  if (/^>(FBtr[0-9]*)\s/)    {
     $id = $1;   next }
  $utr5{$id} .= $_;   }

close hand2;

open (hand3, "dmel-all-three_prime_UTR-r6.07.fasta") or die $!;
my %utr3; 

while (<hand3>)     {
  chomp;
  if (/^>(FBtr[0-9]*)\s/)    {
     $id = $1;   next }
  $utr3{$id} .= $_;   }

close hand3;


open (OUT, ">dmel-all-trans-100.txt");

foreach my $ge (sort keys %mrna)   {
  my $gene_size = length($mrna{$ge});
  my $utr5; my $utr3;

  if (exists $utr5{$ge})  {
    $utr5 = length($utr5{$ge});}
  else {next}  

  if (exists $utr3{$ge})    {
    $utr3 = length($utr3{$ge});}
  else {next}

  my $seq = substr($mrna{$ge}, $utr5-100, $gene_size-$utr5-$utr3+200);
  
  next if length($seq) < 100;

  next if substr($seq, 100, 3) ne 'ATG';
 
  print OUT ">$ge\n$seq\n";        }
   
 close OUT;
