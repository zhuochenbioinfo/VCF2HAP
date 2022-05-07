use strict;
use warnings;

my $hap_file = shift;
my $var_file = shift;
my $sample_info = shift;
my $out_file = shift;

my $usage = "USAGE:\nperl $0 <hap> <var> <sample info> <out prefix>\n";
$usage .= "InDels will be transfrom to SNP and may lead to false results.\n";

die $usage unless(defined $out_file);

my $ntax = 0;
my $nchar = 0;
my $ntrait = 0;
my $groups_joint;
my $hap_data;
my $group_data;

open(IN,"<$sample_info") or die $!;
my %hash_sample;
my %hash_group;
while(<IN>){
	chomp;
	my($sample,$group) = split/\t/;
	$hash_sample{$sample}{group} = $group;
	$hash_group{$group}{count} ++;
}
close IN;

my @groups = sort keys %hash_group;
$ntrait = @groups;
$groups_joint = join(" ",@groups);

my %hash_var = ();
my $rank = 0;
open(IN,"<$var_file") or die $!;
while(<IN>){
    chomp;
    next if($_ =~ /^#/);
    my($chr,$pos,$ref_base,$alt_joint,$cov_count,$alt_count,$snpeff) = split/\t/;
    my @alt_bases = split/,/,$alt_joint;
	my @bases = ($ref_base,@alt_bases);
	# transform indel to A and T
	my $is_indel = 0;
	foreach my $base(@bases){
		next if(length($base) == 1);
		$is_indel = 1;
	}
	if($is_indel == 1){
		$bases[0] = "A";
		for(my $i = 1; $i < @bases; $i++){
			$bases[$i] = "T";
		}
	}
    @{$hash_var{$rank}{bases}} = @bases;
    $rank++;
}
close IN;

$nchar = $rank;

open(IN,"<$hap_file") or die $!;

my %hash_hap;
my $id = 0;

while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($symbol,$hap_rank,$haplotype,$anno,$num,$samples_joint) = split/\t/;
	my @samples = split/,/, $samples_joint;
    my @codes = split/\|/,$haplotype;
    unless(@codes == $rank){
        die "ERROR: number of vars in the var file and hapltype file is note the same.\n";
    }
	my @bases = ();
    for(my $i = 0; $i < @codes; $i++){
        my $code = $codes[$i];
        my $base = $code;
        if($code ne "-" and $code ne "h"){
            $base = ${$hash_var{$i}{bases}}[$code];
        }else{
			$base = "-";
		}
        push @bases, $base;
    }
	my $seq = join("",@bases);
    foreach my $sample(@samples){
		next unless(exists $hash_sample{$sample});
		my $group = $hash_sample{$sample}{group};
		$hash_hap{$id}{$group}{count} ++;
	}
	$hap_data .= "Hap_$id ".join("",@bases)."\n";
	$id++;
}
close IN;

$ntax = $id;

foreach my $i(sort {$a <=> $b} keys %hash_hap){
	my @values = ();
	foreach my $group(@groups){
		my $value = 0;
		if(exists $hash_hap{$i}{$group}){
			$value = $hash_hap{$i}{$group}{count};
		}
		push @values, $value;
	}
	$group_data .= "Hap_$i ".join(",", @values)."\n";
}

open(OUT,">$out_file.nex");

print OUT <<SET;
#NEXUS

Begin Data;
Dimensions ntax=$ntax nchar=$nchar;
Format datatype=DNA missing=N gap=-;
Matrix
$hap_data
;
END;

Begin Traits;
Dimensions NTraits=$ntrait;
format labels=yes missing=? separator=Comma;
TraitLabels $groups_joint;
Matrix
$group_data
;
End;
SET
close OUT;

