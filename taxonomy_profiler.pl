#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use experimental 'smartmatch';
use Getopt::Std;
use Encode;

my %options=();
getopts("hm:e:c:do:t:fa", \%options);

#DB files
my $taxdump_dir=$options{t} || "/home/manan/NCL/kaiju_decodes/ncbi_tax_dmp/new_taxdump/";
my $fullnamelineagefile=$taxdump_dir."/taxidlineage.dmp";
my $nodesfile=$taxdump_dir."/nodes.dmp";
my $namesfile=$taxdump_dir."/names.dmp";

#input file names
my $metadata_file=$options{m} ||"/home/ncim/ganga_sampling/EPI2ME_CSV/barcodes_to_keep.tsv";
my $extract_list=$options{e} ||"/home/ncim/ganga_sampling/EPI2ME_CSV/extract.txt";
my $csv_dir=$options{c} ||"/home/ncim/ganga_sampling/EPI2ME_CSV/csv/";

#input_parameters
my $demultiplex_req= $options{d} ? 1:0;
my $keywordextract_req=$options{e} ? 1:0;
my $amr_req=$options{a} ? 1:0;
if ($keywordextract_req==0){$extract_list="extract not required"}
my $outdir=$options{o} ||"/home/ncim/ganga_sampling/test/";


#output Files names
my $lineage_read=$outdir."/read_lineage.tsv";
my $level1_outfile=$outdir."/Domain.tsv";
my $level2_outfile=$outdir."/phylum.tsv";
my $level3_outfile=$outdir."/class.tsv";
my $level4_outfile=$outdir."/order.tsv";
my $level5_outfile=$outdir."/family.tsv";
my $level6_outfile=$outdir."/genus.tsv";
my $classified_outfile=$outdir."/classified.tsv";
my @created_dirs;
my @created_files;

#variable definitions
my %lineage;
my %profile_l1;
my %profile_l2;
my %profile_l3;
my %profile_l4;
my %profile_l5;
my %profile_l6;
my @barcodes_list;
my %classified_count;
my %unclassified_count;
my %metadata;
my $epi2me_uploaded="";
my @metadata_headers;
my %keycheck;
my @samples_list;
my @extract_list;
my $state;
my %nodes;
my %names;

print"
Input Data:
Taxonomy DMP folder from NCBI:$taxdump_dir
Metadata File: $metadata_file
CSV folder:$csv_dir
Demultiplexing:$demultiplex_req
Extract req:$keywordextract_req
Extract keywords list file:$extract_list
Carryout_AMR against MegaresDB:$amr_req
Output Directory:$outdir

Starting Analysis\n\n";



print "______________________Loading Metadata file__________________\n\n";
open(my $samples_metadata,"<",$metadata_file);
binmode $samples_metadata, ":utf8";
foreach my $line (<$samples_metadata>){
	$line=decode_utf8($line);
	$line=~s/\r\n//;
	$line=~s/\n//;
	if($line=~/^#/){
		@metadata_headers=split /\t/,$line;
		$metadata_headers[0]=~s /^#//;
	}else{
		my @entry=split /\t/,$line;
		if (not defined $entry[2]){}else{
			my $keycheck_key=$entry[2].";".$entry[4];
			#print "$keycheck_key\t $entry[0]\n";
			push(@samples_list,$entry[0]);
			$keycheck{$keycheck_key}=$entry[0];
			for (my $i=0;$i<scalar(@entry);$i++){
				chomp($entry[$i]);
				my $value=$entry[$i];
				my $header=$metadata_headers[$i];
				$metadata{$keycheck_key}{$header}=$value;
				#print "$keycheck_key:$header:$value\n";
			}
			#print "______\n"
		}
	}
}

my $barcode_count=scalar(keys %metadata);
@samples_list=uniq @samples_list;
my $sample_count=scalar(@samples_list);
print "loaded metadata: @metadata_headers\nfor $sample_count samples from $barcode_count barcodes\n\n";

if ($keywordextract_req==1){
	print "______________________Loading Extract List__________________\n\n";
	my $extractlist=`cat $extract_list`;
	@extract_list=split /\n/, $extractlist;
	print "loaded extract list into memory: @extract_list\n";
}
print "______________________Loading Taxonomy file__________________\n\n";

#lineage
open (my $linfile,"<",$fullnamelineagefile);
foreach my $line(<$linfile>){
	my @entry=split /\|/,$line;
	$entry[0]=~s/\s+$//g;
	$entry[1]=~s/\s+$//g;
	$entry[1]=~s/^\s+//g;
	$lineage{$entry[0]}=$entry[1];
}
close $linfile;
my $db_count=scalar(keys %lineage);
print "Loaded $db_count lineage entries in memory\n\n";


#names
$db_count=0;
open (my $namfile,"<",$namesfile);
foreach my $line(<$namfile>){
	my @entry=split /\|/,$line;
	if ($entry[3]=~/scientific name/){
		$entry[0]=~s/\s+$//g;
		$entry[1]=~s/\s+$//g;
		$entry[1]=~s/^\s+//g;
		$names{$entry[0]}=$entry[1];
	}
}
close $namfile;
$db_count=scalar(keys %names);
print "Loaded $db_count names entries in memory\n\n";

#nodes
$db_count=0;
open (my $nodfile,"<",$nodesfile);
foreach my $line(<$nodfile>){
	my @entry=split /\|/,$line;
	$entry[0]=~s/\s+$//g;
	$entry[2]=~s/\s+$//g;
	$entry[2]=~s/^\s+//g;
	$nodes{$entry[0]}=$entry[2];
}
close $nodfile;
$db_count=scalar(keys %nodes);
print "Loaded $db_count nodes entries in memory\n\n";
#foreach (keys %nodes){
#	print "$_\t$nodes{$_}\n";
#}


print "_____________________Processing epi2me csv files___________________________";
if (exists $options{f}){
	print "\nare you sure you want to delete previously created folder (press CTRL+C to exit)";
	#chomp(my $r=<>);
	my $rem=`rm -r $outdir`;
	print "\n removed previously created folder";
}
`mkdir $outdir`;
`mkdir $outdir/segregated`;

opendir my $dir, $csv_dir or die "Cannot open directory: $csv_dir $!";
my @in_files = readdir $dir;
closedir $dir;

open (my $out,">", $lineage_read);
foreach my $csvfile (sort @in_files){
	if ($csvfile eq ".." || $csvfile eq "."){next;}
	print "\n-->Starting $csvfile";
	open (my $csv,"<", $csv_dir.$csvfile);
	binmode $csv, ":utf8";
	my $dummy=<$csv>;
	while (my $line=<$csv>){
		$line=decode_utf8($line);
		my $lineage ="";
		my $state="";
		my @entry=split /,/,$line;
		my $keycheck_key=$entry[2].";".$entry[4];
		#print "\n$keycheck_key";
		if (exists $keycheck{$keycheck_key}){
			my $barcode=$keycheck_key;
			#print "--get\n";
			if(@barcodes_list~~$barcode){}else{push(@barcodes_list,$barcode);}
			if($entry[3]=~/Classified/){

				my $taxid=$entry[5];
				if ($taxid == 1){unclassified_counter($barcode);next;}
				if(exists $lineage{$taxid}){$lineage=$lineage{$taxid}." $taxid"}else{$lineage=check_merged($taxid)." $taxid";}
				$lineage=~s/\D/\;/g;
				$lineage=~s/131567//;
				$lineage=~s/^\;//;
				if ($lineage eq ""){unclassified_counter($barcode);next;}else{
					#print "$lineage\t";
					$lineage=taxonomy_refiner($lineage);
					print $out "$entry[1]\t$barcode\t$lineage\n";
					add_to_profiles($barcode,$lineage);
					classified_counter($barcode);
				}
			}else{
				#print "-----";
				unclassified_counter($barcode);
			}
			if($demultiplex_req==1){
				my $readid=$entry[1];
				my $demultiplexed_dir=$outdir."/demultiplexed/";
				my $sample_dir=$demultiplexed_dir.$keycheck{$barcode}."/";
				my $fastq_sample=$demultiplexed_dir.$keycheck{$barcode}.".fastq";
				my $fastqfilename=$entry[0];
				my $fastqfile=$metadata{$barcode}{"Epi2me_uploads_folders"}."/".$fastqfilename;
				#print "$readid $fastqfile\n";
				my $fastqentry=`grep -m 1 -A3 $readid $fastqfile`;
				if ($demultiplexed_dir~~@created_dirs){}else{`mkdir $demultiplexed_dir`;push(@created_dirs,$demultiplexed_dir);}
				#if ($sample_dir~~@created_dirs){}else{`mkdir $sample_dir`;push(@created_dirs,$sample_dir);}
				open (my $sf,">>",$fastq_sample);
				print $sf $fastqentry;
				close $sf;
			}
			if($keywordextract_req==1){
				my $edir=$outdir."/extracted";
				mkdir $edir;
				foreach my $eterm (@extract_list){
					my $eedir=$edir."/".$eterm;
					if ($eedir~~@created_dirs){}else{`mkdir $eedir`;push(@created_dirs,$eedir);}
					if($lineage=~/$eterm/){
						my $orgname=$entry[6];
						$orgname=~s/\s/_/g;
						$orgname=~s/\(/_/g;
						$orgname=~s/\)/_/g;
						my $org_dir=$eedir."/".$orgname;
						if ($org_dir~~@created_dirs){}else{`mkdir $org_dir`;push(@created_dirs,$org_dir);}
						my $readid=$entry[1];
						my $outfasta=$org_dir."/".$orgname."_".$keycheck{$barcode}.".fastq";
						my $fastqfilename=$entry[0];
						my $fastqfile=$metadata{$barcode}{"Epi2me_uploads_folders"}."*";
						#my @metadata_feilds=keys $metadata{$barcode};
						#print "barcode: $barcode\t fastqfile is $fastqfile\n";
						`grep -A3 $readid $fastqfile >>$outfasta`;
						#print "--`grep -A3 $readid $fastqfile >>$outfasta`\n";
					}
				}
			}
		}
	}
	close $csv;
	print "--done\n";
}
close $out;

open (my $classi, ">", $classified_outfile);
@barcodes_list= uniq @barcodes_list;

print $classi "barcode\t";
foreach my $bc(sort @barcodes_list){
	print $classi "$bc\t";
}
print $classi "\nsampleID\t";
foreach my $bc(sort @barcodes_list){
	print $classi "$keycheck{$bc}\t";
}
print $classi "\nClassified\t";
foreach my $bc(sort @barcodes_list){
	print $classi "$classified_count{$bc}\t";
}
print $classi "\nUnclassified\t";
foreach my $bc(sort @barcodes_list){
	print $classi "$unclassified_count{$bc}\t";
}
close $classi;



print "\n\n______________________writingclassification_______________\n-->L1";
my $header_all;
$header_all="sample\t";
open (my $out1,">", $level1_outfile);
print $out1 "sample\t";
foreach my $bc(sort @samples_list){
	print $out1 "$bc\t";
	$header_all=$header_all."$bc\t";
}
print $out1 "\n";
foreach my $org (sort keys %profile_l1){
	print $out1 "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l1_score;
		if (exists$profile_l1{$org}{$bc}){$profile_l1_score=$profile_l1{$org}{$bc}}else{$profile_l1_score=0}
		print $out1 "$profile_l1_score\t";
	}
	print $out1 "\n";
}
close $out1;

print "--done\n-->L2";
open (my $out2,">", $level2_outfile);
print $out2 "sample\t";
foreach my $bc(sort @samples_list){
	print $out2 "$bc\t";
}
print $out2 "\n";
foreach my $org (sort keys %profile_l2){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l2_score;
		if (exists$profile_l2{$org}{$bc}){$profile_l2_score=$profile_l2{$org}{$bc}}else{$profile_l2_score=0}
		$line=$line."$profile_l2_score\t";
	}
	level_wise_writer("level_2",$line);
	print $out2 "$line\n";
}
close $out2;

print "--done\n-->L3";
open (my $out3,">", $level3_outfile);
print $out3 "sample\t";
foreach my $bc(sort @samples_list){
	print $out3 "$bc\t";
}
print $out3 "\n";
foreach my $org (sort keys %profile_l3){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l3_score;
		if (exists$profile_l3{$org}{$bc}){$profile_l3_score=$profile_l3{$org}{$bc}}else{$profile_l3_score=0}
		$line=$line."$profile_l3_score\t";
	}
	level_wise_writer("level_3",$line);
	print $out3 "$line\n";
}
close $out3;

print "--done\n-->L4";
open (my $out4,">", $level4_outfile);
print $out4 "sample\t";
foreach my $bc(sort @samples_list){
	print $out4 "$bc\t";
}
print $out4 "\n";
foreach my $org (sort keys %profile_l4){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l4_score;
		if (exists$profile_l4{$org}{$bc}){$profile_l4_score=$profile_l4{$org}{$bc}}else{$profile_l4_score=0}
		$line=$line."$profile_l4_score\t";
	}
	level_wise_writer("level_4",$line);
	print $out4 "$line\n";
}
close $out4;

print "--done\n-->L5";
open (my $out5,">", $level5_outfile);
print $out5 "sample\t";
foreach my $bc(sort @samples_list){
	print $out5 "$bc\t";
}
print $out5 "\n";
foreach my $org (sort keys %profile_l5){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l5_score;
		if (exists$profile_l5{$org}{$bc}){$profile_l5_score=$profile_l5{$org}{$bc}}else{$profile_l5_score=0}
		$line=$line."$profile_l5_score\t";
	}
	level_wise_writer("level_5",$line);
	print $out5 "$line\n";
}
close $out5;

print "--done\n-->L6";
open (my $out6,">", $level6_outfile);
print $out6 "sample\t";
foreach my $bc(sort @samples_list){
	print $out6 "$bc\t";
}
print $out6 "\n";
foreach my $org (sort keys %profile_l6){
	my $line= "$org\t";
	foreach my $bc(sort @samples_list){
		my $profile_l6_score;
		if (exists$profile_l6{$org}{$bc}){$profile_l6_score=$profile_l6{$org}{$bc}}else{$profile_l6_score=0}
		$line=$line."$profile_l6_score\t";
	}
	level_wise_writer("level_6",$line);
	print $out6 "$line\n";
}
close $out6;
print "done\n";

#==========Finding location of running code to use in next step===========
use Cwd 'abs_path';
my $code_loc=abs_path($0);
my $run_loc;
if($code_loc=~/(\/.*)(\/.*\.pl$)/){$run_loc=$1;}
#========================================================================

if ($amr_req==1){
	my $demultiplexed_dir=$outdir."/demultiplexed/";
	opendir my $dir, $demultiplexed_dir or die "demultiplexing required for amr to run rerun with -d (it will take exponentially more time)\n major error  $demultiplexed_dir not found $!";
	my @amrfiles = readdir $dir;
	closedir $dir;
	my $amrout_folder=$outdir."/amr_out/";
	my $dem_fasta=$outdir."/dem_fasta/";
	`mkdir $amrout_folder`;
	`mkdir $dem_fasta`;
	foreach my $amrfile(@amrfiles){
		if ($amrfile eq "."||$amrfile eq ".."){next;}
		my $inputfile_fq=$demultiplexed_dir."/".$amrfile;
		$amrfile=~s/\.fastq//;
		my $inputfile_fa=$dem_fasta."/".$amrfile.".fasta";
		my $bp_outfile=$amrout_folder.$amrfile."/";
		`awk '/^@/{gsub(/^@/,">");print;getline;print}' $inputfile_fq >$inputfile_fa`;
		my @args=($inputfile_fa,$bp_outfile,"133.33",$amrfile);
		system($run_loc."/blast_parse_amr.pl",@args)
	}
	my $freqfile=$amrout_folder."/frequency.txt";
	my $amrprofile=$amrout_folder."/amr_profile.tsv";
	my @args2=($freqfile,$amrprofile);
	system($run_loc."/parse_freq.pl",@args2)
}
print "created dirs: @created_dirs";

print "Fininshed Processing\n Good Bye\n";




#Sub Routine Definitions

sub taxonomy_refiner{
	my $lineage=$_[0];
	my @lineage=split /\;/,$lineage;
	my $newlineage;
	if ($names{$lineage[0]} eq "Bacteria"){
		my %lineage;
		my @keys=("superkingdom","phylum","class","order","family","genus","species");
		foreach my $x (@keys){
			$lineage{$x}="-";
		}
		foreach my $i(@lineage){
			if ($nodes{$i}~~@keys){
				$lineage{$nodes{$i}}=$names{$i};
			}
		}
		for (my $y=0;$y<scalar(@keys);$y++){
			if ($lineage{$keys[$y]} eq "-"){$lineage{$keys[$y]}="other_".$lineage{$keys[$y-1]}}
			$newlineage=$newlineage."$lineage{$keys[$y]}; ";
		}
		$newlineage=~s/; $//;
		#print "$newlineage\n\n";
	}
	if ($names{$lineage[0]} eq "Viruses"){
		my %lineage;
		my @keys=("superkingdom","order","family","genus","species");
		foreach my $x (@keys){
			$lineage{$x}="-";
		}
		foreach my $i(@lineage){
			if ($nodes{$i}~~@keys){
				$lineage{$nodes{$i}}=$names{$i};
			}
		}
		for (my $y=0;$y<scalar(@keys);$y++){
			if ($lineage{$keys[$y]} eq "-"){$lineage{$keys[$y]}="other_".$lineage{$keys[$y-1]}}
			$newlineage=$newlineage."$lineage{$keys[$y]}; ";
		}
		$newlineage=~s/; $//;
		#print "$newlineage\n\n";
	}
	if ($names{$lineage[0]} eq "Eukaryota"){
		my %lineage;
		my @keys=("superkingdom","kingdom","phylum","class","order","family","species");
		foreach my $x (@keys){
			$lineage{$x}="-";
		}
		foreach my $i(@lineage){
			if ($nodes{$i}~~@keys){
				$lineage{$nodes{$i}}=$names{$i};
			}
		}
		for (my $y=0;$y<scalar(@keys);$y++){
			if ($lineage{$keys[$y]} eq "-"){$lineage{$keys[$y]}="other_".$lineage{$keys[$y-1]}}
			$newlineage=$newlineage."$lineage{$keys[$y]}; ";
		}
		$newlineage=~s/; $//;
		#print "$newlineage\n\n";
	}
	if ($names{$lineage[0]} eq "Archaea"){
		my %lineage;
		my @keys=("superkingdom","phylum","class","order","family","genus","species");
		foreach my $x (@keys){
			$lineage{$x}="-";
		}
		foreach my $i(@lineage){
			if ($nodes{$i}~~@keys){
				$lineage{$nodes{$i}}=$names{$i};
			}
		}
		for (my $y=0;$y<scalar(@keys);$y++){
			if ($lineage{$keys[$y]} eq "-"){$lineage{$keys[$y]}="other_".$lineage{$keys[$y-1]}}
			$newlineage=$newlineage."$lineage{$keys[$y]}; ";
		}
		$newlineage=~s/; $//;
		#print "$newlineage\n\n";
	}

	return $newlineage;
}


sub level_wise_writer{
	my $level=$_[0];
	my $line=$_[1];
	#print "$line\n\n";
	if ($line=~/^(\w+);.*/){
		my $out_folder_name=$outdir."/segregated/".$1."/";
		#print "\n$out_folder_name--";
		if ($out_folder_name~~@created_dirs){}else{
			`mkdir $out_folder_name`;
			push (@created_dirs,$out_folder_name);
			@created_dirs=uniq @created_dirs;
		}
		my $outfilename=$out_folder_name.$level.".tsv";
		if($outfilename~~@created_files){}else{
			open (my $file,">",$out_folder_name.$level.".tsv");
			print $file "$header_all\n";
			close $file;
			push(@created_files,$outfilename);
		}
		open (my $file,">>",$out_folder_name.$level.".tsv");
		print $file "$line\n";
		close $file;
	}
}
#print "@created_dirs\n\n\n";

sub classified_counter{
	my $barcode=$_[0];
	if (exists $classified_count{$barcode}){
		my $tmp= $classified_count{$barcode};
		$tmp++;
		$classified_count{$barcode}=$tmp;
	}else{
		$classified_count{$barcode}=1;
	}
}


sub unclassified_counter{
	my $barcode=$_[0];
	if (exists $unclassified_count{$barcode}){
		my $tmp= $unclassified_count{$barcode};
		$tmp++;
		$unclassified_count{$barcode}=$tmp;
	}else{
		$unclassified_count{$barcode}=1;
	}
}

sub check_merged{
	my $taxid=$_[0];
	my $mergedfile=$taxdump_dir."/merged.dmp";
	my $merline= `grep "^$taxid\t" $mergedfile`;
	my $lineage="";
	if ($merline=~/\S+\s+\|\s+(\S+)\s+\|/){$lineage=$lineage{$1}}
	return $lineage;
}

sub add_to_profiles{
	my $barcode=$_[0];
	my $lineage=$_[1];
	my @lineage=split/\;/,$lineage;
	if (exists $lineage[1]){}else{$lineage[1]="other_".$lineage[0]}
	if (exists $lineage[2]){}else{$lineage[2]="other_".$lineage[1]}
	if (exists $lineage[3]){}else{$lineage[3]="other_".$lineage[2]}
	if (exists $lineage[4]){}else{$lineage[4]="other_".$lineage[3]}
	if (exists $lineage[5]){}else{$lineage[5]="other_".$lineage[4]}
	my $level_1=$lineage[0];
	my $level_2=$lineage[0].";".$lineage[1];
	my $level_3=$lineage[0].";".$lineage[1].";".$lineage[2];
	my $level_4=$lineage[0].";".$lineage[1].";".$lineage[2].";".$lineage[3];
	my $level_5=$lineage[0].";".$lineage[1].";".$lineage[2].";".$lineage[3].";".$lineage[4];
	my $level_6=$lineage[0].";".$lineage[1].";".$lineage[2].";".$lineage[3].";".$lineage[4].";".$lineage[5];
	if (exists $profile_l6{$level_6}{$keycheck{$barcode}}){
		my $tmp=$profile_l6{$level_6}{$keycheck{$barcode}};
		$tmp++;
		$profile_l6{$level_6}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l6{$level_6}{$keycheck{$barcode}}=1;
	}

	if (exists $profile_l5{$level_5}{$keycheck{$barcode}}){
		my $tmp=$profile_l5{$level_5}{$keycheck{$barcode}};
		$tmp++;
		$profile_l5{$level_5}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l5{$level_5}{$keycheck{$barcode}}=1;
	}


	if (exists $profile_l4{$level_4}{$keycheck{$barcode}}){
		my $tmp=$profile_l4{$level_4}{$keycheck{$barcode}};
		$tmp++;
		$profile_l4{$level_4}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l4{$level_4}{$keycheck{$barcode}}=1;
	}



	if (exists $profile_l3{$level_3}{$keycheck{$barcode}}){
		my $tmp=$profile_l3{$level_3}{$keycheck{$barcode}};
		$tmp++;
		$profile_l3{$level_3}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l3{$level_3}{$keycheck{$barcode}}=1;
	}


	if (exists $profile_l2{$level_2}{$keycheck{$barcode}}){
		my $tmp=$profile_l2{$level_2}{$keycheck{$barcode}};
		$tmp++;
		$profile_l2{$level_2}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l2{$level_2}{$keycheck{$barcode}}=1;
	}


	if (exists $profile_l1{$level_1}{$keycheck{$barcode}}){
		my $tmp=$profile_l1{$level_1}{$keycheck{$barcode}};
		$tmp++;
		$profile_l1{$level_1}{$keycheck{$barcode}}=$tmp;
	}else{
		$profile_l1{$level_1}{$keycheck{$barcode}}=1;
	}
}