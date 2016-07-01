use strict;
use warnings;

#set extend region 

my $ext=50;

#get the sequence of tRNA based on the annotation gtf

# step1, get the genome sequences
my $ge; my %genome=();
open (hand1, "NC10.fa") or die $!;
while (<hand1>)                {
$_ =~ s/\s+$//;
if (/>/)                   {
$_ =~ s/>//; $ge=$_; next; }
$_ =~ tr/agct/AGCT/;
$genome{$ge}.= $_;             }
close hand1;


&trna ('merged_asm');

sub trna    {

my ($sample)=@_;

open (hand1, "$sample/merged.gff3") or die $!;
my %trna=(); my $id=" "; my %ext5; my %ext3;
while (<hand1>)    {
chomp;
 next if /^#/;
 my @a=split /\t/;
 
 if ($a[2] eq 'transcript')          { 
   if ($a[8] =~ /ID=(TCONS_\d+);/) {
    $id=$1; print "$id\n";   }

   $ext5{$id."\t".$a[6]} = substr($genome{$a[0]},$a[3]-1-$ext,$ext) if exists $genome{$a[0]};
   $ext3{$id."\t".$a[6]} = substr($genome{$a[0]},$a[4]-1,$ext) if exists $genome{$a[0]};
   next                                 }

 if ($a[2] eq 'exon')      {
   $trna{$id."\t".$a[6]}.=substr($genome{$a[0]},$a[3]-1,$a[4]-$a[3]+1) if exists $genome{$a[0]};             
               }
                   }   

close hand1;


open (tRNA, ">putative_tRNA_sequences/$sample.tRNA.fa");
open (tRNAext, ">putative_tRNA_sequences/$sample.tRNA.ext$ext.fa");


my $gene=scalar keys %trna; print "$gene\n";

foreach my $i (sort keys %trna)  {
  my $seq = $trna{$i};
  my $ext_seq = $ext5{$i}.$trna{$i}.$ext3{$i} if exists $ext5{$i} and exists $ext3{$i};

  my ($tRNA,$strand)=split /\t/, $i;
  if ($strand eq '-')                       {
      $seq=reverse $seq; $seq =~ tr/ATCG/TAGC/;
      $ext_seq=reverse $ext_seq; $ext_seq =~ tr/ATCG/TAGC/;
      print tRNA ">$tRNA\n$seq\n";
      print tRNAext ">$tRNA\n$ext_seq\n";   }
  else {
      print tRNA ">$tRNA\n$seq\n";
      print tRNAext ">$tRNA\n$ext_seq\n";   } }

close tRNA;
close tRNAext;
      
       }
 

