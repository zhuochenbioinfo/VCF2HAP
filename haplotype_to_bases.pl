use strict;
use warnings;

my $hap_file = shift;
my $var_file = shift;
my $out_file = shift;

my $usage = "USAGE:\nperl $0 <haplotype file> <vars file> <output table>\n";
die $usage unless(defined $out_file);

my %hash_var = ();
my $rank = 0;
open(IN,"<$var_file") or die $!;
while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($chr,$pos,$ref_base,$alt_joint,$cov_count,$alt_count,$snpeff) = split/\t/;
	my @alt_bases = split/,/,$alt_joint;
	@{$hash_var{$rank}{bases}} = ($ref_base,@alt_bases);
	$rank++;
}
close IN;

open(IN,"<$hap_file") or die $!;
open(OUT,">$out_file");
while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($symbol,$haprank,$haplotype,$anno,$num,$samples_joint) = split/\t/;
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
		}
		push @bases, $base;
	}
	print OUT "$symbol\t$haprank\t".join("\t",@bases)."\n";
}
close IN;
close OUT;
