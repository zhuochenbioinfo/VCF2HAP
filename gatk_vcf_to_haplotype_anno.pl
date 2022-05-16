# by Zhuo Chen, IGDB, CAS
# email1: zhuochen@genetics.ac.cn
# email2: zhuochenbioinfo@gmail.com

use strict;
use warnings;
use Getopt::Long;

my($vcf,$regionlist,$keeplist,$noindel,$nointron,$nohet,$maf,$mincov,$outpath,$pass,$preferAlt);
my $usage = "USAGE:\nperl $0 --vcf <vcf> --reg <region list> --keep <sample list> --out <outpath> --maf <maf> --cov <minium cover> --noindel\n";
$usage .= "<vcf> is the input vcf file. [Necessary]\n";
$usage .= "<region list> is a list of regions such as genes. [Necessary]\n";
$usage .= "<outpath> is the path to output haplotype files. [Necessary]\n";
$usage .= "<sample list> is a list of samples. [Optional]\n";
$usage .= "<maf> minor allele frequency. [Optional]\n";
$usage .= "<minium cover> minium coverage ratio of a variant to be picked.\n";
$usage .= "Use --noindel to ignore indels. [Optional]\n";
$usage .= "Use --nointron to ignore variants in intron region without any predictable effect.\n";
$usage .= "Use --pass to pick variants with ONLY 'PASS' or 'SnpCluster' tag.\n";
$usage .= "Use --alt to select the second type in a hetero position.\n";
$usage .= "Use --nohet to ignore hetero position.\n";
$usage .= "Note: variants must be annotated by snpEff version SnpEff 3.3 (build 2013-05-30).\n";

GetOptions(
	"vcf=s" => \$vcf,
	"reg=s" => \$regionlist,
	"keep=s" => \$keeplist,
	"out=s" => \$outpath,
	"maf=s" => \$maf,
	"cov=s" => \$mincov,
	"noindel!" => \$noindel,
	"nointron!" => \$nointron,
	"pass!" => \$pass,
	"alt!" => \$preferAlt,
	"nohet!" => \$nohet,
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
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@datas) = split/\t/;
	
	# filter PASS
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
	
	foreach my $sample(@keepsamples){
		my $i = $hash_sample{$sample}{rank};
		my $spot = $datas[$i];
		my $gt = "-";
		if($spot =~ /^(\d+)\/(\d+)/){
			$gt = $1;
			if(defined $preferAlt){
				$gt = $2;
			}
			if($1 ne $2 and defined $nohet){
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
		my @effs = extract_anno($info,$tag);
		
		# throw a SNP/InDel if it locate in intron with no effect
		if(defined $nointron){
			my $pick = 0;
			foreach(@effs){
				next if($_ =~ /^INTRON/);
				$pick++;
			}
			next if($pick == 0);
		}
		
		my @alleffs = ("REF:REF",@effs);
		push @{$hash_tag{$tag}{vars}}, "$chr;$pos;$ref;$alt;$covcount;$altcount;".join(",",@effs);
		foreach my $sample(@keepsamples){
			my $eff = "NA";
			my $i = $hash_sample{$sample}{rank};
			my $spot = $datas[$i];
			my $gt = "-";
			if($spot =~ /^(\d+)\/(\d+)/){
				$gt = $1;
				if(defined $preferAlt){
					$gt = $2;
				}
				if($1 ne $2 and defined $nohet){
					$gt = "h";
				}
			}
			push @{$hash_sample{$sample}{tag}{$tag}{gt}}, $gt;
			if($gt ne "-" and $gt ne "h" and defined $alleffs[$gt]){
				($eff) = $alleffs[$gt] =~ /(\S+)\:/;
			}
			$hash_sample{$sample}{tag}{$tag}{eff}{$eff}++;
		}
	}
}
close VCF;

print "# Constructing haplotype...\n";

# construct haplotype

my @sorted_tags = sort {$hash_tag{$a}{chrnum} <=> $hash_tag{$b}{chrnum} or $hash_tag{$a}{start} <=> $hash_tag{$b}{start}} keys %hash_tag;
my @impacts	= qw(HIGH MODERATE LOW MODIFIER REF NA);

for(my $i = 0; $i < @sorted_tags; $i++){
	my $tag = $sorted_tags[$i];
	
	unless(exists $hash_tag{$tag}{vars}){
		print "\t# No variants in tag:$tag.\n";
		splice(@sorted_tags,$i,1);
		$i--;
		next;
	}
	
	my $varcount = @{$hash_tag{$tag}{vars}};
	my $hap_ref = "0";
	$hap_ref .= "|0" x ($varcount - 1);
	$hash_haplotype{$tag}{$hap_ref}{count} = 0;
	$hash_haplotype{$tag}{$hap_ref}{rank} = 0;
	@{$hash_haplotype{$tag}{$hap_ref}{samples}} = ();
	$hash_haplotype{$tag}{$hap_ref}{eff} = "HIGH:0;MODERATE:0;LOW:0;MODIFIER:0;REF:$varcount;NA:0";
	
	foreach my $sample(sort keys %hash_sample){
		if($varcount != @{$hash_sample{$sample}{tag}{$tag}{gt}}){
			die "# ERROR: variant count not identical at gene: $tag\n";
		}
		my $hap = join("|", @{$hash_sample{$sample}{tag}{$tag}{gt}});
		push @{$hash_haplotype{$tag}{$hap}{samples}}, $sample;
		$hash_haplotype{$tag}{$hap}{count}++;
		unless(exists $hash_haplotype{$tag}{$hap}{eff}){
			my @effcounts = ();
			foreach my $eff(@impacts){
				my $effcount = 0;
				if(exists $hash_sample{$sample}{tag}{$tag}{eff}{$eff}){
					$effcount = $hash_sample{$sample}{tag}{$tag}{eff}{$eff};
				}
				push @effcounts, "$eff:$effcount";
			}
			$hash_haplotype{$tag}{$hap}{eff} = join(";",@effcounts);
		}
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
	print OUT "#CHROM\tPOS\tREF\tALT\tCOV\tALT\tIMPACT\n";
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
		my $hapeff = $hash_haplotype{$tag}{$hap}{eff};
		print OUT "$tag\t$rank\t$hap\t$hapeff\t$count\t".join(",",@{$hash_haplotype{$tag}{$hap}{samples}})."\n";
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

sub extract_anno{
	my($info,$tag) = @_;
	$info =~ /EFF=(\S+)/;
	my $eff = $1;
	my @EFF = split/,/,$eff;
	
	my $effect = "NA";
	my @effs = ();
	
	my %hash_impact;
	# HIGH, MODERATE, LOW, MODIFIER
	$hash_impact{HIGH} = 3;
	$hash_impact{MODERATE} = 2;
	$hash_impact{LOW} = 1;
	$hash_impact{MODIFIER} = 0;
	$hash_impact{INTRON} = -1;
	
	foreach(@EFF){
		$_ =~ /(\S+?)\((\S+)\)/;
		my $effect_tmp = $1;
		my $desc_tmp = $2;
		my @descs = split/\|/,$desc_tmp;
		
		my $impact = $descs[0];
		my $locus = $descs[5];
		my $transcript = $descs[8];
		my $altrank = $descs[10] - 1;
		
		next if($locus ne $tag);
		
		if(defined $effs[$altrank]){
			my($efftmp) = $effs[$altrank] =~ /(\S+?)\:/;
			next unless($hash_impact{$impact} > $hash_impact{$efftmp});
		}
		$effs[$altrank] = "$impact:$_";
		if($effect_tmp eq "INTRON" and $impact eq "MODIFIER"){
			$effs[$altrank] = "INTRON:$_";
		}
		#print "$effect\n$desc\n";
		
		# UPSTREAM(MODIFIER|||||Gene_SSIV1-Kasalath.t01|||SSIV1-Kasalath.t01||1)
		# NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|gTc/gCc|V469A||Gene_SSIV1-Kasalath.t01|||SSIV1-Kasalath.t01|4|1)
		# UPSTREAM(MODIFIER|||||LOC_Os12g39830|||LOC_Os12g39830.1||1)
		# INTRON(MODIFIER|||||Gene_SSIV1-Kasalath.t01|||SSIV1-Kasalath.t01|4|1)
		# SPLICE_SITE_DONOR(HIGH|||||LOC_Os06g04200|||LOC_Os06g04200.1|1|1)
		# UTR_5_PRIME(MODIFIER|||||LOC_Os06g04200|||LOC_Os06g04200.3|1|1),UTR_5_PRIME(MODIFIER|||||LOC_Os06g04200|||LOC_Os06g04200.3|1|2)
	}
	return @effs;
}
