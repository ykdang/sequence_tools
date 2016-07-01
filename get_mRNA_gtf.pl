use strict;
use warnings;

#set extend region 

my $ext=5;

#get the sequence of mRNA based on gtf
# step1, get the genome sequences
my $ge; my %genome=();
open (hand1, "dmel-all-chromosome-r5.45.fasta") or die $!;
while (<hand1>)                {
$_ =~ s/\s+$//;
if (/>(\S+)\stype/)                   {
$ge = "chr$1";
print "$ge\n"; next; }
$_ =~ tr/agct/AGCT/;
$genome{$ge}.= $_;             }
close hand1;


open (hand1, "MDv3_annotation.gtf") or die $!;
my %mrna=(); my $id=" "; 
while (<hand1>)    {
 chomp;
 my @a=split /\t/;
 next if $a[2] ne "CDS";
 my @b=split /"/, $a[8]; 
 if ($id eq $b[3])    {               #$b[3] is TCONS, this is to ensure that if two exons from the same gene(TCONS), then it will be connected. the offset of gtf file is 1, not 0 
     my $exon=substr($genome{$a[0]},$a[3]-1,$a[4]-$a[3]+1) if exists $genome{$a[0]};
     $mrna{$id."\t".$a[6]}.=$exon;  
          }
 else{  
     $id=$b[3];
     $mrna{$id."\t".$a[6]}=substr($genome{$a[0]},$a[3]-1,$a[4]-$a[3]+1) if exists $genome{$a[0]};
      }   }

close hand1;


open (mRNA, ">dmelr5.45.mRNA.largest.fa");

my $gene=scalar keys %mrna; print "$gene\n";
my %uni_seq;
foreach my $i (keys %mrna)  {
  my $size = length($mrna{$i}); 
  my $seq = $mrna{$i};
  my ($mRNA,$strand)=split /\t/, $i;
  $mRNA =~ s/\.\w+$//;
  $mRNA =~ s/-R\w+$//;
  #print "$mRNA\n";
  if ($strand eq '-')   {
      $seq=reverse $seq; $seq =~ tr/ATCG/TAGC/;
      }
  if (defined $uni_seq{$mRNA} && length($uni_seq{$mRNA}) <= $size)   {
           $uni_seq{$mRNA} = $seq;
            }
        else  {
           $uni_seq{$mRNA} = $seq;
            }
   }

foreach my $i (sort keys %uni_seq)   {
  print mRNA ">$i\n$uni_seq{$i}\n";    }
close mRNA;

          
        
       


