use strict;
use warnings;
use List::Util 'shuffle';

my $infile = $ARGV[0];

my $number_of_samples = $ARGV[1];

chomp $ARGV[1];

open(IN, $infile);

my @vector; #will store file lines
my @indices; #will store random positions to sample;
my $i = 0; #index

while (my $line = <IN>) {
  chomp $line;
  $vector[$i] = $line;
  $i++;
}

close IN;

#creating a vector of lines to be sampled;

my %seen ;
my $i2 = 0;
while( $i2 <= ($number_of_samples - 1)) {
  my $x = int rand($i-1);
  $seen{$x} and next ;
  $seen{$x} = 1;
  print $vector[$x]."\n";
  $i2++;
}

#print("@indices\n");

