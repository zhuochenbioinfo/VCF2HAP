# by Zhuo Chen, IGDB, CAS
# email1: zhuochen@genetics.ac.cn
# email2: zhuochenbioinfo@gmail.com

use strict;
use warnings;
use Getopt::Long;

my($vcf,$varlist,$nohet,$outpath);
my $usage = "USAGE:\nperl $0 --vcf <vcf> --var <variants list> --out <outpath>\n";
$usage .= "<vcf> is the input vcf file. [Necessary]\n";
$usage .= "<variants list> is a list of variants.[Necessary]\n";
$usage .= "\t#CHROM\tPOS\tREF\tALT\tTAGs\n";
$usage .= "<outpath> is the path to output haplotype files. [Necessary]\n";
$usage .= "Use --nohet to ignore hetero position.\n";

GetOptions(
	"vcf=s" => \$vcf,
	"var=s" => \$varlist,
	"out=s" => \$outpath,
	"nohet!" => \$nohet,
) or die $usage;

die $usage unless(defined $vcf and defined $varlist and defined $outpath);

my %hash_tag;
my %hash_sample;
my %hash_var;
my %hash_haplotype;
my %hash_chr; # This hash is used to check the end of the region list and stop reading vcf

print "# Reading variants list...\n";
open(IN,"<$varlist") or die $!;
while(<IN>){
	chomp;
	my($chr,$pos,$ref,$alt,$tag) = split/\t/;
	push @{$hash_var{$chr}{$pos}{$ref}{$alt}{tags}}, $tag;
	push @{$hash_tag{$tag}{vars}}, "$chr;$pos;$ref;$alt";
	unless(exists $hash_chr{$chr}){
		$hash_chr{$chr}{end} = $pos;
	}
	if($pos > $hash_chr{$chr}{end}){
		$hash_chr{$chr}{end} = $pos;
	}
}
close IN;

print "# Reading VCF...\n";

if($vcf =~ /gz$|gzip$/){
	open(VCF,"zcat $vcf|") or die $!;
}else{
	open(VCF,"<$vcf") or die $!;
}

my @keepsamples = ();
my @remained_chrs = ();
my $chr_tmp = "";
my @regions = ();
my $skip = 0;
my $chrend = 0;

while(<VCF>){
	chomp;
	
	# pick sample names
	if($_ =~ /^#CHROM/){
		my(undef,undef,undef,undef,undef,undef,undef,undef,undef,@samplenames) = split/\t/;
		for(my $i = 0; $i < @samplenames; $i++){
			$hash_sample{$samplenames[$i]}{rank} = $i;
			push @keepsamples, $samplenames[$i];
		}
		print "\t# sample number=".@keepsamples."\n";
		next;
	}
	
	next if($_ =~ /^#/);
	
	# split vcf line
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@datas) = split/\t/;
	
	# pick regions when entering a new chr
	if($chr_tmp ne $chr){
		$chr_tmp = $chr;
		if(exists $hash_chr{$chr}){
			$chrend = $hash_chr{$chr}{end};
			delete($hash_chr{$chr});
			print "\t# Reading chr:$chr\n";
			$skip = 0;
		}else{
			$chrend = 0;
			print "\t# Skipping chr:$chr\n";
			$skip = 1;
		}
		@remained_chrs = keys %hash_chr;
	}
	
	
	# check if all the regions have been read through
	if(@remained_chrs == 0 and $skip == 1){
		last;
	}
	if(@remained_chrs == 0 and $pos > $chrend){
		last;
	}
	
	# check if the position has reach the end of the chromosome region with variants
	if($skip == 1 or $pos > $chrend){
		next;
	}
	
	# check if the position is contained in the variants list
	next unless(exists $hash_var{$chr});
	next unless(exists $hash_var{$chr}{$pos});
	next unless(exists $hash_var{$chr}{$pos}{$ref});
	
	foreach my $alts_join(keys %{$hash_var{$chr}{$pos}{$ref}}){
		my %hash_tmp;
		my @alts = split/,/,$alts_join;
		my @alts_new = split/,/,$alt;
		my @bases_new = ($ref,@alts_new);
		my @bases = ($ref,@alts);
		for(my $i = 0; $i < @bases; $i++){
			$hash_tmp{$bases[$i]}{code} = $i;
		}
		foreach my $sample(@keepsamples){
			my $i = $hash_sample{$sample}{rank};
			my $spot = $datas[$i];
			my $gt = "-";
			if($spot =~ /^(\d+)\/(\d+)/){
				$gt = $1;
				my $base = $bases_new[$gt];
				if(exists $hash_tmp{$base}){
					$gt = $hash_tmp{$base}{code};
				}else{
					$gt = "-";
				}
				if($1 ne $2 and defined $nohet){
					$gt = "h";
				}
			}
			$hash_var{$chr}{$pos}{$ref}{$alts_join}{sample}{$sample} = $gt;
		}
	}
}
close VCF;

print "# Constructing haplotype...\n";

# construct haplotype

my @sorted_tags = sort keys %hash_tag;
for(my $i = 0; $i < @sorted_tags; $i++){
	my $tag = $sorted_tags[$i];
	my @vars = @{$hash_tag{$tag}{vars}};
	foreach my $var(@vars){
		my($chr,$pos,$ref,$alt) = split/;/,$var;
		foreach my $sample(@keepsamples){
			my $gt = 0;
			if(exists $hash_var{$chr}{$pos}{$ref}{$alt}{sample}{$sample}){
				$gt = $hash_var{$chr}{$pos}{$ref}{$alt}{sample}{$sample};
			}
			push @{$hash_sample{$sample}{tag}{$tag}}, $gt;
		}
	}
	
	my $varcount = 0;
	foreach my $sample(sort keys %hash_sample){
		if($varcount == 0){
			if(exists $hash_sample{$sample}{tag}{$tag}){
				$varcount = @{$hash_sample{$sample}{tag}{$tag}};
			}else{
				splice(@sorted_tags, $i, 1);
				$i--;
				last;
			}
			my $hap_ref = "0";
			$hap_ref .= "|0" x ($varcount - 1);
			$hash_haplotype{$tag}{$hap_ref}{count} = 0;
			$hash_haplotype{$tag}{$hap_ref}{rank} = 0;
			@{$hash_haplotype{$tag}{$hap_ref}{samples}} = ();
		}
		if($varcount != @{$hash_sample{$sample}{tag}{$tag}}){
			die "# ERROR: variant count not identical at gene: $tag\n";
		}
		my $hap = join("|", @{$hash_sample{$sample}{tag}{$tag}});
		push @{$hash_haplotype{$tag}{$hap}{samples}}, $sample;
		$hash_haplotype{$tag}{$hap}{count}++;
	}
}

# rank the haplotypes
foreach my $tag(@sorted_tags){
	my $count = 0;
	foreach my $hap(sort {$hash_haplotype{$tag}{$b}{count} <=> $hash_haplotype{$tag}{$a}{count}} keys %{$hash_haplotype{$tag}}){
		next if(exists $hash_haplotype{$tag}{$hap}{rank});
		$count++;
		$hash_haplotype{$tag}{$hap}{rank} = $count;
	}
}

print "# Outputing results...\n";

# output haplotypes
#open(OUT,">$outfile");

foreach my $tag(@sorted_tags){
	open(OUT,">$outpath/$tag.vars");
	print OUT "#samplenum=".@keepsamples."\n";
	print OUT "#tag=$tag\n";
	print OUT "#varnum=".@{$hash_tag{$tag}{vars}}."\n";
	print OUT "#CHROM\tPOS\tREF\tALT\tCOV\tALT\n";
	foreach(@{$hash_tag{$tag}{vars}}){
		my $vars = $_;
		$vars =~ s/;/\t/g;
		print OUT "$vars\n";
	}
	open(OUT,">$outpath/$tag.haplotype");
	foreach my $hap(sort {$hash_haplotype{$tag}{$a}{rank} <=> $hash_haplotype{$tag}{$b}{rank}} keys %{$hash_haplotype{$tag}}){
		my $count = $hash_haplotype{$tag}{$hap}{count};
		my $rank = $hash_haplotype{$tag}{$hap}{rank};
		print OUT "$tag\t$rank\t$hap\t-\t$count\t".join(",",@{$hash_haplotype{$tag}{$hap}{samples}})."\n";
	}
	close OUT;
}


sub uniq{
	my @items = @_;
	my %hash;
	my @outarr = ();
	foreach(@items){
		next if(exists $hash{$_});
		push @outarr, $_;
		$hash{$_} = "";
	}
	return @outarr;
}
