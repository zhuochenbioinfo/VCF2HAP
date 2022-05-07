# by Zhuo Chen, IGDB, CAS
# email1: zhuochen@genetics.ac.cn
# email2: zhuochenbioinfo@gmail.com

use strict;
use warnings;
#use Text::NSP::Measures::2D::Fisher::right;

my($hapfile,$sampleinfo,$outfile) = @ARGV;

my $usage = "USAGE:\nperl $0 <haplotype file> <sample info> <outfile>\n";
$usage .= "<haplotype file>: #TAG\tRANK\tHAPLOTYPE\tNOTE\tSAMPLES\n";
$usage .= "<sample info>: #SAMPLE\tGROUP\n";

die $usage unless(@ARGV == 3);

my %hash_sample;
my %hash_group;

open(IN,"<$sampleinfo") or die $!;
while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($sample,$subtype) = split/\t/;
	if($sample =~ /^Sample_/){
		$sample =~ s/Sample_//;
	}
	$hash_sample{$sample}{group} = $subtype;
	push @{$hash_group{$subtype}{samples}}, $sample;
}
close IN;

my @groups = sort keys %hash_group;
my $samplenum = keys %hash_sample;
foreach my $subtype(@groups){
	my $subtypecount = @{$hash_group{$subtype}{samples}};
	my $ratio = $subtypecount/$samplenum;
	$hash_group{$subtype}{ratio} = $ratio;
}

#push @groups, "unknown";

open(IN,"<$hapfile") or die $!;
open(OUT,">$outfile") or die $!;

print OUT "#TAG\tRANK\tHAPLOTYPE\tNOTE\tSUM\t".join("\t",@groups)."\n";

while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($tag,$rank,$haplotype,$note,$sum,$samples_join) = split/\t/;
	my @samples = split/,/,$samples_join;
	my @ratios = ();
	my %hash_tmp;
	foreach my $sample(@samples){
		if($sample =~ /^Sample_/){
			$sample =~ s/Sample_//;
		}
		unless(exists $hash_sample{$sample}){
			$sum--;
			next;
		}
		my $group = $hash_sample{$sample}{group};
		$hash_tmp{$group}{count} ++;
	}
	next if($sum == 0);
	foreach my $group(@groups){
		unless(exists $hash_tmp{$group}){
			$hash_tmp{$group}{count} = 0;
		}
		my $ratio = $hash_tmp{$group}{count}/@{$hash_group{$group}{samples}};
		$ratio = sprintf("%.3f", $ratio);
		push @ratios, $hash_tmp{$group}{count};
	}
	print OUT "$tag\t$rank\t$haplotype\t$note\t$sum\t".join("\t",@ratios)."\n";
}
close IN;
close OUT;

