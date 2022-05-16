# by Zhuo Chen, IGDB, CAS
# email1: zhuochen@genetics.ac.cn
# email2: zhuochenbioinfo@gmail.com

use strict;
use warnings;
use Getopt::Long;

my($vcf,$regionlist,$keeplist,$noindel,$maf,$mincov,$outpath,$nohet,$pass);
my $usage = "USAGE:\nperl $0 --vcf <vcf> --reg <region list> --keep <sample list> --out <outpath> --cov <minium cover> --maf <maf> --noindel\n";
$usage .= "<vcf> is the input vcf file. [Necessary]\n";
$usage .= "<region list> is a list of regions such as genes. [Necessary]\n";
$usage .= "<outpath> is the path to output haplotype files. [Necessary]\n";
$usage .= "<sample list> is a list of samples. [Optional]\n";
$usage .= "<maf> minor allele frequency. [Optional]\n";
$usage .= "<minium cover> minium coverage ratio of a variant to be picked.\n";
$usage .= "Use --noindel to ignore indels. [Optional]\n";
$usage .= "Use --nohet to mark heterozygous genotype by \"h\".\n";
$usage .= "Use --pass to pick variants with ONLY 'PASS' or 'SnpCluster' tag.\n";

GetOptions(
	"vcf=s" => \$vcf,
	"reg=s" => \$regionlist,
	"keep=s" => \$keeplist,
	"out=s" => \$outpath,
	"maf=s" => \$maf,
	"cov=s" => \$mincov,
	"noindel!" => \$noindel,
	"nohet!" => \$nohet,
	"pass!" => \$pass,
) or die $usage;

die $usage unless(defined $vcf and defined $regionlist and defined $outpath);

my %hash_reg;
my %hash_tag;
my %hash_sample;
my %hash_haplotype;
my %hash_chr; # This hash is used to check the end of the region list and stop reading vcf

print "# Reading region list...\n";

open(IN,"<$regionlist") or die $!;

=cu
#EXAMPLE:
Chr1	1000	2500	gene1
Chr1	2760	2900	gene1
Chr1	3000	3200	gene2

=cut

while(<IN>){
	chomp;
	my($chr,$start,$end,$tag) = split/\t/;
	push @{$hash_reg{$chr}{$start}{$end}{tag}}, $tag;
	my($chrnum) = $chr =~ /Chr(\d+)/;
	if(exists $hash_tag{$tag}{chr}){
		if($hash_tag{$tag}{chr} ne $chr){
			die "# ERROR: the regions in one tag shall all be on the same chromosome!\n";
		}
	}
	$hash_tag{$tag}{chr} = $chr;
	$hash_tag{$tag}{chrnum} = $chrnum;
	unless(exists $hash_tag{$tag}{start}){
		$hash_tag{$tag}{start} = $start;
		$hash_tag{$tag}{end} = $end;
	}
	if($start < $hash_tag{$tag}{start}){
		$hash_tag{$tag}{start} = $start;
	}
	if($end > $hash_tag{$tag}{end}){
		$hash_tag{$tag}{end} = $end;
	}
	$hash_chr{$chr} = "";
}
close IN;

if(defined $keeplist){
	print "# Reading sample list...\n";
	open(IN,"<$keeplist") or die $!;
	while(<IN>){
		chomp;
		$hash_sample{$_}{rank} = "";
	}
	close IN;
}

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

while(<VCF>){
	chomp;
	
	# pick sample names
	if($_ =~ /^#CHROM/){
		my(undef,undef,undef,undef,undef,undef,undef,undef,undef,@samplenames) = split/\t/;
		if(defined $keeplist){
			for(my $i = 0; $i < @samplenames; $i++){
				next unless(exists $hash_sample{$samplenames[$i]});
				$hash_sample{$samplenames[$i]}{rank} = $i;
				push @keepsamples, $samplenames[$i];
			}
		}else{
			for(my $i = 0; $i < @samplenames; $i++){
				$hash_sample{$samplenames[$i]}{rank} = $i;
				push @keepsamples, $samplenames[$i];
			}
		}
		print "\t# sample number=".@keepsamples."\n";
		next;
	}
	
	next if($_ =~ /^#/);
	
	# split vcf line
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$datas_join) = split/\t/,$_,10;
	
	# check pass
	if(defined $pass){
		next unless($filter eq "PASS" or $filter eq "SnpCluster");
	}

	# pick regions when entering a new chr
	if($chr_tmp ne $chr){
		@regions = ();
		foreach my $start(sort {$a <=> $b} keys %{$hash_reg{$chr}}){
			foreach my $end(sort {$a <=> $b} keys %{$hash_reg{$chr}{$start}}){
				push @regions, "$start,$end";
			}
		}
		$chr_tmp = $chr;
		if(exists $hash_chr{$chr}){
			delete($hash_chr{$chr});
		}
		@remained_chrs = keys %hash_chr;
		print "\t# Reading chr:$chr\n";
		if(@regions == 0){
			print "\t# Skipping chr:$chr\n";
		}
	}
	
	# check if all the regions have been read through
	if(@remained_chrs == 0 and @regions == 0){
		last;
	}
	
	# check if the chr contains any regions
	if(@regions == 0){
		next;
	}
	
	# check if the position is within any region. If so, pick the tags(genes) for the position
	my @tags = ();
	
	for(my $i = 0; $i < @regions; $i++){
		my($start,$end) = split/,/,$regions[$i];
		if($pos < $start){
			last;
		}elsif($pos > $end){
			splice(@regions, $i, 1);
			$i--;
			next;
		}else{
			push @tags, @{$hash_reg{$chr}{$start}{$end}{tag}};
		}
	}
	next unless(@tags > 0);
	@tags = uniq(@tags);
	
	
	my @alts = split/,/,$alt;
	my @bases = ($ref, @alts);
	
	# check if the variant is an indel
	my $checkindel = 0;
	if(defined $noindel){
		foreach(@bases){
			if(length($_) > 1){
				$checkindel = 1;
			}
		}
	}
	next if($checkindel == 1);
	
	# check if the variant exists in the selected samples
	my $altnum = @alts;
	my $altcount = 0;
	my $covcount = 0;
	my $lostcount = 0;
	
	my @datas = split/\t/,$datas_join;

	foreach my $sample(@keepsamples){
		my $i = $hash_sample{$sample}{rank};
		my $spot = $datas[$i];
		my $gt = "-";
		if($spot =~ /^(\d+)\/(\d+)/){
			$gt = $1;
			if(defined $nohet and $1 ne $2){
				$gt = "h";
			}
		}
		if($gt ne '-' and $gt ne "h"){
			$covcount++;
		}else{
			$lostcount++;
		}
		if($gt =~ /\d/ and $gt > 0){
			$altcount++;
		}
	}
	
	next if($altcount == 0);
	
	my $allcount = $lostcount + $covcount;
	if(defined $mincov){
		next unless($lostcount/$allcount < (1 - $mincov));
	}
	if(defined $maf){
		next unless($altcount/$allcount > $maf and ($covcount - $altcount)/$allcount > $maf);
	}
	
	# push the variant into selected samples
	foreach my $tag(@tags){
		push @{$hash_tag{$tag}{vars}}, "$chr;$pos;$ref;$alt;$covcount;$altcount";
		foreach my $sample(@keepsamples){
			my $i = $hash_sample{$sample}{rank};
			my $spot = $datas[$i];
			my $gt = "-";
			if($spot =~ /^(\d+)\/(\d+)/){
				$gt = $1;
				if(defined $nohet and $1 ne $2){
					$gt = "h";
				}
			}
			push @{$hash_sample{$sample}{tag}{$tag}}, $gt;
		}
	}
}
close VCF;

print "# Constructing haplotype...\n";

# construct haplotype

my @sorted_tags = sort {$hash_tag{$a}{chrnum} <=> $hash_tag{$b}{chrnum} or $hash_tag{$a}{start} <=> $hash_tag{$b}{start}} keys %hash_tag;
for(my $i = 0; $i < @sorted_tags; $i++){
	my $tag = $sorted_tags[$i];
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
	close OUT;
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
