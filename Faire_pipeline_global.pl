#!/usr/bin/perl
#fastx_toolkit_global.pl
use warnings;
use strict;
use YAML::Tiny;
use File::Basename;
use Data::Dumper;
use change_format_fork;
#my @files;

use macs_job;
use csar_job;

#fill @files using yaml file
my $config_file = YAML::Tiny->new;
# Open the config
$config_file= YAML::Tiny->read("conect_script_global_try2.yml") or die "Cannot read global config file\n";

#==============Manage Fastx-toolkit and save file group features=============================================
my $count = 0;
my %file_CorT;
unlink 'output_files.log' if (-e 'output_files.log');
foreach my $yml_file_ref (@{$config_file->[0]->{'fastx_toolkit_files'}}) {
	foreach my $yml_file (keys %{$yml_file_ref}) {
		warn "Fastx-toolkit\tUSING FILE: $yml_file...\n";
		!system "perl fastx-toolkit.pl $yml_file" or die "Cannot execute fast-toolkit.pl: $!";
		my $type_file = $config_file->[0]->{'fastx_toolkit_files'}->[$count]->{$yml_file};
		if ($type_file) {
			$file_CorT{$yml_file} = $type_file;
		}
		$count += 1;
	}
}
#print "file_CorT hash\n";
#print Dumper %file_CorT;
#==============Manage bowtie build process of multiple genomes ====================
foreach my $yml_file (@{$config_file->[0]->{'bowtie-build'}}) {
	if ($yml_file) {
		warn "bowtie-build\tUSING FILE: $yml_file...\n";
		#!system "perl bowtie_build.pl $yml_file" or die "Cannot execute bowtie-build.pl: $!";
	}
	else {
		warn "bowtie-build execution canceled : no fasta file available for this part";
	}
}

#==============
my $fastx_config_file = YAML::Tiny->new;
my %yml_infile;
foreach (keys %file_CorT) {
	# Open the fastx-toolkit config
	$fastx_config_file= YAML::Tiny->read("$_");
	my $infile = $fastx_config_file->[0]->{'in_file'};
	$yml_infile{$_} = $infile;
}

#==========Get raw file => clean file relation==============
#===== using output_files.log========================
my %rawfile_cleanfile;
open LOG_FILES,'<','output_files.log' or die "Cannot read output_files.log: $!";
while (<LOG_FILES>) {
	chomp;
	my (undef,$rawfile,$cleanfile) = split ('\t',$_);
	$rawfile_cleanfile{$rawfile} = $cleanfile;
}
my %cleanfile_rawfile = reverse %rawfile_cleanfile;
#print "rawfile_cleanfile hash\n";
#print Dumper %rawfile_cleanfile;

#==================
$count = 0;
my %CorT_organelleBowtie;
foreach my $organelle_bowtie_ref (@{$config_file->[0]->{'bowtie-organelle'}}) {
	foreach my $organelle_bowtie_yml (keys %{$organelle_bowtie_ref}) {
		my $CorT_organelle_bowtie = $config_file->[0]->{'bowtie-organelle'}->[$count]->{$organelle_bowtie_yml};
		$CorT_organelleBowtie{$CorT_organelle_bowtie}= $organelle_bowtie_yml;
	}
	$count += 1;
}
$count = 0;
my %CorT_nuclearBowtie;
foreach my $nuclear_bowtie_ref (@{$config_file->[0]->{'bowtie-nuclear'}}) {
	foreach my $nuclear_bowtie_yml (keys %{$nuclear_bowtie_ref}) {
		my $CorT_nuclear_bowtie = $config_file->[0]->{'bowtie-nuclear'}->[$count]->{$nuclear_bowtie_yml};
		$CorT_nuclearBowtie{$CorT_nuclear_bowtie}= $nuclear_bowtie_yml;
	}
	$count += 1;
}
my %organelleBowtie_CorT = reverse %CorT_organelleBowtie;
my %nuclearBowtie_CorT = reverse %CorT_nuclearBowtie;
#print "organelleBowtie_CorT hash\n";
#print Dumper %organelleBowtie_CorT;
#=====================

my $organelle_bowtie_config_file = YAML::Tiny->new;
my $nuclear_bowtie_config_file = YAML::Tiny->new;

# Open the config
#$organelle_bowtie_config_file_control = YAML::Tiny->read("$organelle_bowtie_file");
my %CorT_file = reverse %file_CorT;
#========= 

#============Bowtie Mapping to organelle reference Genomes manage I/O options ===============
foreach my $organelle_yml_file (keys %organelleBowtie_CorT) {
	if($organelle_yml_file) {
		$organelle_bowtie_config_file = YAML::Tiny->read($organelle_yml_file);
		#print "organelle_bowtie_config_file:  $organelle_bowtie_config_file\n";
		#print "organelle_yml_file: $organelle_yml_file\n";
		foreach my $fastx_yml_file (keys %file_CorT) {
			my $fastx_group_str = $file_CorT{$fastx_yml_file};
			my ($fastx_group) = $fastx_group_str =~ /^([A-Z]+)_[0-9]+$/i;
			#print "$fastx_group | $organelleBowtie_CorT{$organelle_yml_file}\n";
			if (lc($fastx_group) eq lc($organelleBowtie_CorT{$organelle_yml_file})) {
				#print "ENTER\n";
				
				for (1..4) {
					#print "1-4: $_\n";
					
					foreach my $main_cat (keys %{$organelle_bowtie_config_file->[0]->{'1_MainArguments'}->{'main'}->[$_]}) {
						my $count_fastqc_files = 0;
						#print "main_cat: $main_cat\n";
						foreach my $file_label (@{$organelle_bowtie_config_file->[0]->{'1_MainArguments'}->{'main'}->[$_]->{$main_cat}}) {
							if ($file_label and $file_label eq $fastx_group_str) {
								#print "$file_label:SUCCES!!!!\n";
								my $fastx_yml_file = $CorT_file{$file_label};
								my $raw_file = $yml_infile{$fastx_yml_file};
								my $clean_file = $rawfile_cleanfile{$raw_file};
								#print "$fastx_yml_file | $raw_file | $clean_file\n";
								$organelle_bowtie_config_file->[0]->{'1_MainArguments'}->{'main'}->[$_]->{$main_cat}->[$count_fastqc_files]= "$clean_file";
								#print "COUNT: $count_fastqc_files\n";
								#$organelle_bowtie_config_file->write($organelle_yml_file);
							}
							$count_fastqc_files += 1;
						}
						#print "COUNT OUT: $count_fastqc_files\n";
					}
				}
			}
		}
		print "\n";
		my $base_name  = $organelle_bowtie_config_file->[0]->{'0_General_system'}->{'sys'}->[5]->{'outfile_base_name'};
		$organelle_bowtie_config_file->[0]->{'5_Output'}->{'values'}->[2]->{'un'} = $base_name."_unaligned";
		$organelle_bowtie_config_file->write($organelle_yml_file);
		!system "perl config_YAML_file_spaces.pl $organelle_yml_file" or die "Cannot fortmat YAML bowtie organelle files: $!";
		warn "bowtie mapping for $organelle_yml_file starts...\n";
		!system "perl bowtie_job.pl $organelle_yml_file" or die "Cannot execute bowtie_job: $!";
		warn "\n";
		foreach my $nuclear_yml_file (keys %nuclearBowtie_CorT) {
			if ($nuclear_yml_file and $nuclearBowtie_CorT{$nuclear_yml_file} eq $organelleBowtie_CorT{$organelle_yml_file}) {
				$nuclear_bowtie_config_file = YAML::Tiny->read($nuclear_yml_file);
				foreach my $unaligned_file (glob "$base_name\_unaligned*") {
					if ($unaligned_file ~~ /_1[.fq]?$/) {
						$nuclear_bowtie_config_file->[0]->{'1_MainArguments'}->{'main'}->[1]->{'m1'}->[0]= "$unaligned_file";
					}
					elsif ($unaligned_file ~~ /_2[.fq]?$/) {
						$nuclear_bowtie_config_file->[0]->{'1_MainArguments'}->{'main'}->[2]->{'m2'}->[0]= "$unaligned_file";
					}
					elsif ($unaligned_file ~~ /[^_1-2]$/) {
						$nuclear_bowtie_config_file->[0]->{'1_MainArguments'}->{'main'}->[4]->{'s'}->[0]= "$unaligned_file";
					}
				}
				$nuclear_bowtie_config_file->write($nuclear_yml_file);
				!system "perl config_YAML_file_spaces.pl $nuclear_yml_file" or die "Cannot fortmat YAML bowtie nuclear files: $!";
				warn "bowtie mapping for $nuclear_yml_file starts...\n";
				!system "perl bowtie_job.pl $nuclear_yml_file" or die "Cannot execute bowtie_job: $!";
				warn "\n";
			}
		}
	}
}

#==============================================================================================================

#=============Remove organelle mapping files===========
if ($config_file->[0]->{'format-change-processing'}->[0]->{'erase_intermediate_mapped_files'}) {
	foreach my $bowtie_organelle_config_hash (@{$config_file->[0]->{'bowtie-organelle'}}) {
		foreach my $bowtie_organelle_config (keys %{$bowtie_organelle_config_hash}) {
			#print $bowtie_organelle_config ."\n";
			$organelle_bowtie_config_file = YAML::Tiny->read($bowtie_organelle_config);
			my $outfile_basename = $organelle_bowtie_config_file->[0]->{'0_General_system'}->{'sys'}->[5]->{'outfile_base_name'};
			#print "HERE:$outfile_basename\n";
			foreach my $erase_file (glob "$outfile_basename*[!.report]") {
				#print "FILE:$erase_file\n";
				if ($erase_file ~~ /.+\.sam(\.gz)?$/) {
					unlink $erase_file;
				}
				else {
					print "Could not erase file $erase_file becouse suffix was not .sam or .sam.gz\n";
				}
			}
		}
	}
}
#======================================================


#============Bowtie Mapping to nuclear reference Genomes manage I/O options ===============
unless (%organelleBowtie_CorT) {
	 print "bypass to nuclear bowtie!!!\n";
	 foreach my $nuclear_yml_file (keys %nuclearBowtie_CorT) {
		 if ($nuclear_yml_file) {
			$nuclear_bowtie_config_file = YAML::Tiny->read($nuclear_yml_file);
			foreach my $fastx_yml_file (keys %file_CorT) {
				my $fastx_group_str = $file_CorT{$fastx_yml_file};
				my ($fastx_group) = $fastx_group_str =~ /^([A-Z]+)_[0-9]+$/i;
				#print "$fastx_group | $nuclearBowtie_CorT{$nuclear_yml_file}\n";
				if (lc($fastx_group) eq lc($nuclearBowtie_CorT{$nuclear_yml_file})) {
					#print "ENTER\n";
					
					for (1..4) {
						foreach my $main_cat (keys %{$nuclear_bowtie_config_file->[0]->{'1_MainArguments'}->{'main'}->[$_]}) {
							my $count_fastqc_files = 0;
							foreach my $file_label (@{$nuclear_bowtie_config_file->[0]->{'1_MainArguments'}->{'main'}->[$_]->{$main_cat}}) {
								if ($file_label and $file_label eq $fastx_group_str) {
									#print "$file_label:SUCCES!!!!\n";
									my $fastx_yml_file = $CorT_file{$file_label};
									my $raw_file = $yml_infile{$fastx_yml_file};
									my $clean_file = $rawfile_cleanfile{$raw_file};
									#print "$fastx_yml_file | $raw_file | $clean_file\n";
									$nuclear_bowtie_config_file->[0]->{'1_MainArguments'}->{'main'}->[$_]->{$main_cat}->[$count_fastqc_files]= "$clean_file";
									#print "COUNT: $count_fastqc_files\n";
								}
								$count_fastqc_files += 0;
							}
						}
					}
				}
			}
			print "\n";
			$nuclear_bowtie_config_file->write($nuclear_yml_file);
			!system "perl config_YAML_file_spaces.pl $nuclear_yml_file" or die "Cannot fortmat YAML bowtie nuclear files: $!";
			warn "bowtie mapping for $nuclear_yml_file starts...\n";
			!system "perl bowtie_job.pl $nuclear_yml_file" or die "Cannot execute bowtie_job: $!";
			warn "\n";
		}
	}
}
#================================================================================================


#=============Remove organelle mapping files===========
if ($config_file->[0]->{'format-change-processing'}->[0]->{'erase_intermediate_mapped_files'}) {
	foreach my $bowtie_organelle_config_hash (@{$config_file->[0]->{'bowtie-organelle'}}) {
		foreach my $bowtie_organelle_config (keys %{$bowtie_organelle_config_hash}) {
			#print $bowtie_organelle_config ."\n";
			$organelle_bowtie_config_file = YAML::Tiny->read($bowtie_organelle_config);
			my $outfile_un_basename = $organelle_bowtie_config_file->[0]->{'5_Output'}->{'values'}->[2]->{'un'};
			foreach my $erase_file (glob "$outfile_un_basename*[!.report]") {
				#print "ERASE_FILE:$erase_file.\n";
				unlink $erase_file;
			}
		}	
	}
}
#======================================================

my $erase_sam_files_flag = $config_file->[0]->{'format-change-processing'}->[1]->{'erase_sam_files'};
print "ERASE_FILES:$erase_sam_files_flag\n";
my @samfiles = (glob "*sam.gz");
print "Data format sam -> bed starts...\n";
change_format_fork::sam2bed_fork($erase_sam_files_flag,@samfiles);
print "Data format sam -> bed ends...\n";

print Dumper %nuclearBowtie_CorT;

#=========Specified input files for peak_callers======================================

my @peak_caller_ymls;
foreach my $program_yml_tag (keys %{$config_file->[0]->{'peak_callers'}}) {

	my $program_config_file = $config_file->[0]->{'peak_callers'}->{"$program_yml_tag"};
	
	open YML_OUT,'>',$program_config_file.".mod" or die "Cannot write config mod file\n";

	open YML,'<',$program_config_file or die "Cannot open config file $program_config_file\n";
	push (@peak_caller_ymls,$program_config_file);

	while (<YML>) {
		chomp;
		my $line = $_;
		foreach my $nuclear_bowtie_yml (keys %nuclearBowtie_CorT) {
			#print "HERE:$nuclear_bowtie_yml\t$nuclearBowtie_CorT{$nuclear_bowtie_yml}\n";
			my $bowtie_build_config_file = YAML::Tiny->read($nuclear_bowtie_yml);
			my $out_bowtie_basemane = $bowtie_build_config_file->[0]->{'0_General_system'}->{'sys'}->[5]->{'outfile_base_name'};
			my $bed_file = $out_bowtie_basemane."_sorted.bed";
			my $bam_file = $out_bowtie_basemane."_sorted.bam";

			if ($nuclearBowtie_CorT{$nuclear_bowtie_yml} =~ /treatment/i) {
				if ($line =~ /^\s+treatment:/) {
					print "T_SUCCES\t$line\n";
					if ($program_yml_tag ~~ /csar_config_file/) {
						$line =~ s/(treatment:).+$/$1 $bam_file/;
					} 
					else {
						$line =~ s/(treatment:).+$/$1 $bed_file/;
					}
					print "$line\n";
				}
			}
			elsif ($nuclearBowtie_CorT{$nuclear_bowtie_yml} =~ /control/i) {
				if ($line =~ /^\s+control:/) {
					print "C_SUCCES\t$line\n";
					if ($program_yml_tag ~~ /csar_config_file/) {
						$line =~ s/(control:).+$/control: $bam_file/;
					}
					else {
						$line =~ s/(control:).+$/control: $bed_file/;
					}
					print "$line\n";
				}
			}
			
		}
		print YML_OUT $line."\n";

	}
	close YML;
	close YML_OUT;
}

foreach my $peak_caller_name (@peak_caller_ymls) {
	unlink $peak_caller_name or warn "Cannot erase peak-caller config file : $peak_caller_name: $!";
	rename ($peak_caller_name.'.mod',$peak_caller_name) or warn "Cannot rename $peak_caller_name.mod: $!";
}
#===================================================================================


my $macs_peaks_bed_fl = macs_job::Macs_Job($config_file->[0]->{'peak_callers'}->{'macs_config_file'});
print "PEAKS_FILES:($macs_peaks_bed_fl)";
#my $csar_peaks_bed_fl = csar_job::Csar_Job($config_file->[0]->{'peak_callers'}->{'csar_config_file'});
#TODO: Mappability.pl integration
!system "perl MOSAICS.pl $config_file->[0]->{'peak_callers'}->{'mosaics_config_file'}" or die "Cannort execute MOSAICS.pl: $!";

#print "PEAKS_FILES:($macs_peaks_bed_fl,$csar_peaks_bed_fl)\n";


