use strict;
use warnings;

my $hap_file = shift;
my $var_file = shift;
my $sample_info = shift;
my $out_file = shift;

my $usage = "USAGE:\nperl $0 <hap> <var> <sample info> <out>\n";
$usage .= "InDels will be transfrom to SNP and may lead to false results.\n";

die $usage unless(defined $out_file);

open(IN,"<$sample_info") or die $!;
my %hash_sample;
while(<IN>){
	chomp;
	my($sample,$group) = split/\t/;
	$hash_sample{$sample}{group} = $group;
}
close IN;

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


open(IN,"<$hap_file") or die $!;
open(OUT,">$out_file");

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
		print OUT ">$group\n$seq\n";
	}
}
close IN;
close OUT;
