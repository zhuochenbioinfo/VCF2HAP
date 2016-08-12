use warnings;
use strict;
use Cwd 'abs_path';

my($hapfile,$outfile) = @ARGV;
my $usage = "USAGE:\nperl $0 <haplotype file> <outfile>\n";

die $usage unless(@ARGV == 2);

$hapfile = abs_path($hapfile);

my %hash_haps;
my $gene_tmp = "";

open(IN,"<$hapfile") or die $!;
open(OUT,">$outfile");
print OUT "##hap=$hapfile\n";
print OUT "#TAG\tRANK\tHAPLOTYPE\tNOTE\tSAMPLENUM\tSAMPLES\n";

while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($gene,$rank,$hap,$anno,$num,$samples_join) = split/\t/;
	
	if($gene_tmp ne $gene and $gene_tmp ne ""){
		
		my @ranks = sort {$hash_haps{$gene_tmp}{$a}{lostcount} <=> $hash_haps{$gene_tmp}{$b}{lostcount} or $a <=> $b} keys %{$hash_haps{$gene_tmp}};
		for(my $i = 0; $i < @ranks; $i++){
			my $rank1 = $ranks[$i];
			my $hap1 = $hash_haps{$gene_tmp}{$rank1}{haplotype};
			if(exists $hash_haps{$gene_tmp}{$rank1}{superior}){
				splice(@ranks,$i,1);
				$i--;
				next;
			}
			for(my $j = $i + 1; $j < @ranks; $j++){
				my $rank2 = $ranks[$j];
				my $hap2 = $hash_haps{$gene_tmp}{$rank2}{haplotype};
				my $percentage = compare_hap($hap1,$hap2);
				if($percentage == 1){
					push @{$hash_haps{$gene_tmp}{$rank1}{subtype}}, $rank2;
					push @{$hash_haps{$gene_tmp}{$rank2}{superior}}, $rank1;
				}
				if(exists $hash_haps{$gene_tmp}{$rank2}{superior} and @{$hash_haps{$gene_tmp}{$rank2}{superior}} > 1){
					delete($hash_haps{$gene_tmp}{$rank2});
					splice(@ranks,$j,1);
					$j--;
					next;
				}
			}
		}
		
		@ranks = sort {$a <=> $b} @ranks;
		foreach my $rank1(@ranks){
			my $hap1 = $hash_haps{$gene_tmp}{$rank1}{haplotype};
			my $num1 = $hash_haps{$gene_tmp}{$rank1}{samplenum};
			my $anno1 = $hash_haps{$gene_tmp}{$rank1}{anno};
			my $allnum1 = $num1;
			my @samples1 = @{$hash_haps{$gene_tmp}{$rank1}{samples}};
			foreach my $rank2(@{$hash_haps{$gene_tmp}{$rank1}{subtype}}){
				next unless(exists $hash_haps{$gene_tmp}{$rank2});
				push @samples1, @{$hash_haps{$gene_tmp}{$rank2}{samples}};
				$allnum1 += $hash_haps{$gene_tmp}{$rank2}{samplenum};
			}
			print OUT "$gene_tmp\t$rank1\t$hap1\t$anno1\t$allnum1\t".join(",",@samples1)."\n";
		}
		
		delete($hash_haps{$gene_tmp});
	}
	$gene_tmp = $gene;
	
	my @bases = split/\|/,$hap;
	my $lostcount = 0;
	foreach(@bases){
		if($_ eq "-"){
			$lostcount ++;
		}
	}
	$hash_haps{$gene}{$rank}{lostcount} = $lostcount;
	$hash_haps{$gene}{$rank}{samplenum} = $num;
	$hash_haps{$gene}{$rank}{allsamplenum} = $num;
	$hash_haps{$gene}{$rank}{haplotype} = $hap;
	@{$hash_haps{$gene}{$rank}{samples}} = split/,/,$samples_join;
	$hash_haps{$gene}{$rank}{anno} = $anno;
}
close IN;

my @ranks = sort {$hash_haps{$gene_tmp}{$a}{lostcount} <=> $hash_haps{$gene_tmp}{$b}{lostcount} or $a <=> $b} keys %{$hash_haps{$gene_tmp}};

for(my $i = 0; $i < @ranks; $i++){
	my $rank1 = $ranks[$i];
	my $hap1 = $hash_haps{$gene_tmp}{$rank1}{haplotype};
	if(exists $hash_haps{$gene_tmp}{$rank1}{superior}){
		splice(@ranks,$i,1);
		$i--;
		next;
	}
	for(my $j = $i + 1; $j < @ranks; $j++){
		my $rank2 = $ranks[$j];
		my $hap2 = $hash_haps{$gene_tmp}{$rank2}{haplotype};
		my $percentage = compare_hap($hap1,$hap2);
		if($percentage == 1){
			push @{$hash_haps{$gene_tmp}{$rank1}{subtype}}, $rank2;
			push @{$hash_haps{$gene_tmp}{$rank2}{superior}}, $rank1;
		}
			if(exists $hash_haps{$gene_tmp}{$rank2}{superior} and @{$hash_haps{$gene_tmp}{$rank2}{superior}} > 1){
			delete($hash_haps{$gene_tmp}{$rank2});
			splice(@ranks,$j,1);
			$j--;
			next;
		}
	}
}
		
@ranks = sort {$a <=> $b} @ranks;
foreach my $rank1(@ranks){
	my $hap1 = $hash_haps{$gene_tmp}{$rank1}{haplotype};
	my $num1 = $hash_haps{$gene_tmp}{$rank1}{samplenum};
	my $anno1 = $hash_haps{$gene_tmp}{$rank1}{anno};
	my $allnum1 = $num1;
	my @samples1 = @{$hash_haps{$gene_tmp}{$rank1}{samples}};
	foreach my $rank2(@{$hash_haps{$gene_tmp}{$rank1}{subtype}}){
		next unless(exists $hash_haps{$gene_tmp}{$rank2});
		push @samples1, @{$hash_haps{$gene_tmp}{$rank2}{samples}};
		$allnum1 += $hash_haps{$gene_tmp}{$rank2}{samplenum};
	}
	print OUT "$gene_tmp\t$rank1\t$hap1\t$anno1\t$allnum1\t".join(",",@samples1)."\n";
}

close OUT;

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