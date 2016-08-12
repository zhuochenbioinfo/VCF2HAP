use strict;
use warnings;

my($queryfile,$targetfile,$outpath) = @ARGV;

my $usage = "USAGE:\nperl $0 <query haplotype file> <target haplotype group file> <outpath>\n";
$usage .= "<query haplotype file>: #TAG\tRANK\tHAPLOTYPE\tNOTE\tSUM\tSAMPLES\n";
$usage .= "<target haplotype group file>: #TAG\tRANK\tHAPLOTYPE\tNOTE\tSUM\tMAJORGROUPS\tGROUPS\n";

die $usage unless(@ARGV == 3);

my %hash_target;
my %hash_sample;

open(IN,"<$targetfile") or die $!;
while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($tag,$rank,$haplotype,$note,$sum,$majorgroups,@groups) = split/\t/;
	$hash_target{$tag}{$rank}{haplotype} = $haplotype;
	$hash_target{$tag}{$rank}{sum} = $sum;
	$hash_target{$tag}{$rank}{majorgroups} = $majorgroups;
}
close IN;

open(IN,"<$queryfile") or die $!;
while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($tag,$rank,$haplotype,$note,$sum,$samples_joint) = split/\t/;
	my @samples = split/,/,$samples_joint;
	my $target_rank = "-";
	my $majorgroups = "-";
	my %hash_tmp;
	foreach my $trank(sort {$a <=> $b} keys %{$hash_target{$tag}}){
		$hash_tmp{$trank} = compare_hap($hash_target{$tag}{$trank}{haplotype},$haplotype);
	}
	my @tranks = sort {$hash_tmp{$b} <=> $hash_tmp{$a} or $a <=> $b} keys %hash_tmp;
	if($hash_tmp{$tranks[0]} > 0.75){
		$target_rank = $tranks[0];
		$majorgroups = $hash_target{$tag}{$target_rank}{majorgroups};
	}
	foreach my $sample(@samples){
		$hash_sample{$sample}{$tag}{rank} = $target_rank;
		$hash_sample{$sample}{$tag}{majorgroups} = $majorgroups;
	}
}
close IN;

foreach my $sample(sort keys %hash_sample){
	open(OUT,">$outpath/$sample.hapclass");
	foreach my $tag(sort keys %{$hash_sample{$sample}}){
		print OUT "$tag\t$hash_sample{$sample}{$tag}{rank}\t$hash_sample{$sample}{$tag}{majorgroups}\n";
	}
	close OUT;
}


sub compare_hap{
	my($hapref,$hapin) = @_;
	my @arr_hapref = split/\|/,$hapref;
	my @arr_hapin = split/\|/,$hapin;
	unless(@arr_hapref == @arr_hapin){
		print "the two haplotypes for comparation must be the same length.\nref:$hapref\nin:$hapin\n";
		die;
	}
	my $length_compare = 0;
	my $length_same = 0;
	my $percentage = 0;
	for(my $i = 0;$i < @arr_hapref;$i++){
		#if($arr_hapin[$i] eq '-' or $arr_hapref[$i] eq '-'){
		if($arr_hapin[$i] eq '-'){
			next;
		}
		$length_compare++;
		if($arr_hapin[$i] eq $arr_hapref[$i]){
			$length_same++;
		}
	}
	if($length_compare > 0){
		$percentage = $length_same/$length_compare;
	}
	return $percentage;
}