#!/usr/bin/perl
use strict;
use warnings;

@ARGV or die "Usage: $0 <prokka.faa>";

my %ignore = ('hypothetical protein'=>1);

my @gene;
while (<ARGV>) {
  next unless m/^>(\S+)\s+(.+)$/;
  push @gene, [ $1, $2 ];
}
my $N = scalar(@gene);
print STDERR "Found $N genes.\n";

my $dir = $ARGV;
#my @file = split('/',$dir);
#my $filename = "$file[-1].prokka_pseudos.tsv";

open(FH, '>', "$dir.prokka_pseudos.tsv") or die $!;

print FH "Found $N genes.\n";
print FH "Start\tEnd\tGeneProduct\n";

my $P = 0;
if ($N > 1) {
  for my $i (1 .. $N) {
    my $prod = $gene[$i-1][1];
    if ( !$ignore{$prod} and $gene[$i][1] eq $prod ) {
      print FH "$gene[$i-1][0]\t$gene[$i][0]\t$prod\n";
      print "$gene[$i-1][0] & $gene[$i][0] => $prod\n";
      $P++;
    }
  }
}
print FH "Found $P potential pseudo-genes\n";
close(FH);
print STDERR "Found $P potential pseudo-genes\n";