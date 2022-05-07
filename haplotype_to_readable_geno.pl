use strict;
use warnings;
use Getopt::Long;

my($in,$out,$lof_ranks,$gof_ranks,$het_ranks,$wt_ranks,$high_ranks,$low_ranks,$novel_ranks);

my $usage = "USAGE:\nperl $0 --in <input haplotype file> --out <output geno file>\n";
$usage .= "--lof [loss of function hapranks], comma seperated\n";
$usage .= "--gof [gain of function hapranks], comma seperated\n";
$usage .= "--wt [wild type hapranks], comma seperated\n";
$usage .= "--low [low function hapranks], comma seperated\n";
$usage .= "--high [high function hapranks], comma seperated\n";
$usage .= "--novel [novel allele haprank], comma seperated\n";
$usage .= "--het [heterozygous haprank], comma seperated\n";

GetOptions(
	"in=s" => \$in,
	"out=s" => \$out,
	"lof=s" => \$lof_ranks,
	"gof=s" => \$gof_ranks,
	"wt=s" => \$wt_ranks,
	"high=s" => \$high_ranks,
	"low=s" => \$low_ranks,
	"novel=s" => \$novel_ranks,
	"het=s" => \$het_ranks,
) or die $usage;

die $usage unless(defined $in);

unless(defined $out){
	$out = $in;
	$out =~ s/\.haplotype/\.geno/;
	#die "$out";
}

# define hapranks
my %hash_rank;
if(defined $lof_ranks){
	my @ranks = split/,/,$lof_ranks;
	my %hash = map{$_ => "LoF"} @ranks;
	%hash_rank = (%hash_rank,%hash);
}
if(defined $gof_ranks){
	my @ranks = split/,/,$gof_ranks;
	my %hash = map{$_ => "GoF"} @ranks;
	%hash_rank = (%hash_rank,%hash);
}
if(defined $wt_ranks){
	my @ranks = split/,/,$wt_ranks;
	my %hash = map{$_ => "WT"} @ranks;
	%hash_rank = (%hash_rank,%hash);
}
if(defined $high_ranks){
	my @ranks = split/,/,$high_ranks;
	my %hash = map{$_ => "High"} @ranks;
	%hash_rank = (%hash_rank,%hash);
}
if(defined $low_ranks){
	my @ranks = split/,/,$low_ranks;
	my %hash = map{$_ => "Low"} @ranks;
	%hash_rank = (%hash_rank,%hash);
}
if(defined $novel_ranks){
	my @ranks = split/,/,$novel_ranks;
	my %hash = map{$_ => "Novel"} @ranks;
	%hash_rank = (%hash_rank,%hash);
}
if(defined $het_ranks){
	my @ranks = split/,/,$het_ranks;
	my %hash = map{$_ => "Het"} @ranks;
	%hash_rank = (%hash_rank,%hash);
}

=test
foreach my $haprank(sort keys %hash_rank){
	print "$haprank\t$hash_rank{$haprank}\n";
}
=cut

open(IN,"<$in") or die $!;
open(OUT,">$out") or die "# ERROR: unable to open: $out\n";
print OUT "#SAMPLE\trank\tanno\tfunc\n";

while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($locus,$haprank,$haplotype,$func,$num,$samples_join) = split/\t/;
	my @samples = split/,/,$samples_join;
	my $anno = "NA";
	if(defined $hash_rank{$haprank}){
		$anno = $hash_rank{$haprank};
	}
	foreach my $sample(@samples){
		print OUT "$sample\t$haprank\t$anno\t$func\n";
	}
}
close IN;
close OUT;
