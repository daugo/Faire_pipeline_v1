#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use YAML::Tiny;
use Data::Dumper;

my $bowtie_config_file = shift @ARGV; #config file , the structure of this config file could be generates by parse_bowtie_arg.pl
die "bowtie_config_file ($bowtie_config_file) doesn't exists: $!" unless (-e $bowtie_config_file); #check if config_file exists
die "bowtie_config_file ($bowtie_config_file) doesn't end in .yml" unless ($bowtie_config_file =~ /\.yml$/);

open IN,'<','mapping_bowtie_yaml_args.txt' or die "Cannot open mapping_bowtie_yaml_args.txt:$!"; #mappings file this file have exists in the same folder
my %yamlName_arg; #mapping between yaml names and command line arguments
while (<IN>) {
	chomp;
	my ($yaml_name, $arg) = $_ =~ /^'(.+)'\s=>\s'(.+)'$/;
	$yamlName_arg{$yaml_name} = $arg;
}
#=====Open config_file===========
my $config_file = YAML::Tiny->new;
# Open the config
$config_file= YAML::Tiny->read("$bowtie_config_file") or die "Cannot read YAML config file $bowtie_config_file: $!";

#====Building bowtie command======
my $bowtie_cmd = "bowtie ";
#This block extract most of the information contained in the structure of the YAML config file=========

foreach my $cmd_category (sort {$b cmp $a} keys %{$config_file->[0]}) {
	#print $cmd_category."\n";
	foreach my $cmd_cat_type (keys %{$config_file->[0]->{$cmd_category}}) {
		#print "\t$cmd_cat_type\n";
		my $count = 0;
		foreach my $array_cat_type (@{$config_file->[0]->{$cmd_category}->{$cmd_cat_type}}) {
			#print "\t\t$array_cat_type\n";
			foreach my $key_final (keys %{$array_cat_type}) {
				#print "\t\t$key_final\n";
				my $arg = $config_file->[0]->{$cmd_category}->{$cmd_cat_type}->[$count]->{$key_final};
				#print "\t\t\tbuuu:$arg\n" if $arg;
				
				if ($cmd_cat_type ~~ /^flags$/) {
					$bowtie_cmd = $bowtie_cmd."$yamlName_arg{$key_final} " if $arg; 
				}
				elsif ($cmd_cat_type ~~ /^values$/) {
					$bowtie_cmd = $bowtie_cmd."$yamlName_arg{$key_final} $arg " if $arg;
				}
				elsif ($cmd_cat_type ~~ /^main$/) {
					#bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]
					$bowtie_cmd = $bowtie_cmd."$arg " if ($arg and ($key_final ~~ /^(ebwt|hit)$/));
					
					if ($key_final ~~ /^(m1)|(m2)|(r)|(s)$/) {
						my @keyfinalfiles;
						foreach my $files (@{$arg}) {
							if ($files) {
								push (@keyfinalfiles,$files);
							}
						}
						if (@keyfinalfiles) {
							my $keyfinalfiles = join (',',@keyfinalfiles);
							$bowtie_cmd = $bowtie_cmd."-1 $keyfinalfiles " if ($key_final ~~ /^m1$/);
							$bowtie_cmd = $bowtie_cmd."-2 $keyfinalfiles " if ($key_final ~~ /^m2$/);
							$bowtie_cmd = $bowtie_cmd."--12 $keyfinalfiles " if ($key_final ~~ /^r$/);
							$bowtie_cmd = $bowtie_cmd."$keyfinalfiles " if ($key_final ~~ /^s$/);
						}
					}
				} 
				elsif ($arg and $cmd_cat_type ~~ /^sys$/ and $key_final =~ /bowtie_path/ ) {
					$arg =~ s/\/$//;
					$bowtie_cmd = $arg."/".$bowtie_cmd;
				}
			 
				#print "count: $count\n";
				$count += 1;
			}
		}
	}
}
#my $try = $config_file->[0]->{Input}->{values}->[5]->{trim5};
#print "Try: $try\n";
#print "bowtie_cmd: $bowtie_cmd\n"; 
#my $bowtie_path = $config_file->[0]->{'0_General_system'}->{'sys'}->[0]->{'bowtie_path'}

#===Save the arguments in the YAML config file related with general system configurations=====
my $use_pigz = $config_file->[0]->{'0_General_system'}->{'sys'}->[1]->{'use_pigz'};
my $pigz_installation_path = $config_file->[0]->{'0_General_system'}->{'sys'}->[2]->{'pigz_installation_path'};
my $use_gzip = $config_file->[0]->{'0_General_system'}->{'sys'}->[3]->{'use_gzip'};
my $pigz_threads = $config_file->[0]->{'0_General_system'}->{'sys'}->[4]->{'pigz_threads'};


#====Construct the compress command pigz, gzip or none according to YAML config_file======
my $cmd_compress;
if ($use_pigz) {
	my $pigz_cmd = 'pigz ';
	$pigz_installation_path =~ s/\/$//;
	$pigz_cmd = $pigz_installation_path.'/'.$pigz_cmd if $pigz_installation_path;
	if ($pigz_threads) {
		$pigz_cmd = $pigz_cmd."-p $pigz_threads ";
	}
	else {
		$pigz_cmd = $pigz_cmd.'-p 2';
	}
	$cmd_compress = '| '.$pigz_cmd;
}
elsif ($use_gzip) {
	my $gzip_cmd = '| gzip ';
	$cmd_compress = $gzip_cmd;
}

#=====Output files=========
my $outfile_base_name = $config_file->[0]->{'0_General_system'}->{'sys'}->[5]->{'outfile_base_name'};
my $cmd_report_part = "2> $outfile_base_name.report ";
my $cmd_out = "> $outfile_base_name.sam";

#====Concatenate compress and output options with bowtie main command===========
if ($cmd_compress) {
	$bowtie_cmd = $bowtie_cmd.$cmd_report_part.$cmd_compress.$cmd_out.'.gz';
}
else {
	$bowtie_cmd = $bowtie_cmd.$cmd_report_part.$cmd_out;
}

#=====Finally!!, execute the command======================
warn  "BOWTIE CMD:$bowtie_cmd\n";
system $bowtie_cmd;
open REPORT,'<',"$outfile_base_name.report" or die "Cannot open bowtie report file : $!";
my $check_succes;
while (<REPORT>) {
	chomp;
	if ($_ ~~ /^#\sreads\sprocessed:\s[0-9]+$/) {
		$check_succes = 1;
	}
}
close REPORT;
die "Cannot execute bowtie command: $!" unless $check_succes;
