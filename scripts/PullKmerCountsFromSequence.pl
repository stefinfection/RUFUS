#!/usr/bin/perl

use strict;
my $kgJhash = "/scratch/ucgd/lustre/work/marth/resources/RUFUS/1000G.RUFUSreference.min45.Jhash";
my $jellyfishPath="/uufs/chpc.utah.edu/common/HIPAA/u0746015/marth_software/RUFUS/stock/bin/externals/jellyfish/src/jellyfish_project/bin/jellyfish";
my $seq = $ARGV[0];
my $kmerLength = $ARGV[1];
my $hashFile = $ARGV[2];

print "sequence\t$ARGV[2]\t";
for (my $j = 3; $j < scalar(@ARGV); $j++)
	        {print "$ARGV[$j]\t";}
		print "1kg\n";

for ( my $i = 0; $i < length($seq) - $kmerLength; $i++)
{

    # print kmer
	my $kmer =  substr($seq, $i, $kmerLength);
	print "$kmer\t";

	# print counts from first hashFile
	my $first = `$jellyfishPath  query $hashFile $kmer`;
        chomp $first; 
        my @temp1 = split / /, $first;
        print "$temp1[1]\t";

    # print counts from other hashFiles
	for (my $j = 3; $j < scalar(@ARGV); $j++)
	{
		my $first = `$jellyfishPath  query $ARGV[$j] $kmer`;
		chomp $first; 
		my @temp1 = split / /, $first;
		print "$temp1[1]\t";
	}

    # print counts from 1kg hash
	my $first = `$jellyfishPath  query $kgJhash $kmer`;
	chomp $first;
	my @temp1 = split / /, $first;
	if ($temp1[1] eq 0){
		my $revcomp = reverse $kmer;
		$revcomp =~ tr/ATGCatgc/TACGtacg/;
		$first = `$jellyfishPath  query $kgJhash $revcomp`;
		chomp $first;
		@temp1 = split / /, $first;
	}
	print "$temp1[1]\t";
	print "\n"; 
}
print "\n";
