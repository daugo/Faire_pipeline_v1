#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;
use YAML::Tiny;
my $build_config_file = shift @ARGV;

unless ($build_config_file and -e $build_config_file and $build_config_file ~~ /[.]yml$/) {
	die "usage: Give a valid YAML config file for bowtie-build: $!";
	}
#Mappings YAML labels => command line argument
my $arg_cmd = {
'noauto_index' => '--noauto',
'packed' => '--packed',
'no_rf' => '--noref',
'justref' => '--justref',
'ntoa' => '--ntoa',
'bmax' => '--bmax',
'bmaxdivn' => '--bmaxdivn',
'no_dc' => '--nodc',
'off_rate' => '--offrate', 
'ftabchars' => '--ftabchars',
'seed' => '--seed',
'cutoff' => '--cutoff',
'color_index' => '--color'
};

my %hash_arg_cmd = %$arg_cmd;

my $config_file = YAML::Tiny->new;
my $build_cmd = "bowtie-build ";

# Open the config
$config_file= YAML::Tiny->read($build_config_file) or die "Cannot read $build_config_file: $!";

foreach my $type_arg (sort {$b cmp $a} keys %{$config_file->[0]}) {
	#print "$type_arg\n";
	foreach my $cat_arg (keys %{$config_file->[0]{$type_arg}}) {
		#print "\t$cat_arg\n";
		my $count = 0;
		foreach my $arg_hash (@{$config_file->[0]->{$type_arg}->{$cat_arg}}) {
			
			foreach my $arg (keys %{$arg_hash}) {
				#print "\t\t$count\t";
				#print "$arg\n";
				my $opt = $config_file->[0]->{$type_arg}->{$cat_arg}->[$count]->{$arg};
				if ($cat_arg ~~ /^flags$/ ) {
					#print "$opt\n" if $opt;
					my $values_cmd = $hash_arg_cmd{$arg} if ($opt  && $hash_arg_cmd{$arg});
					$build_cmd = $build_cmd."$values_cmd " if $values_cmd;
				}

				elsif ($cat_arg ~~ /^values$/ && $type_arg !~ /^1_main_arguments$/) {
					#print "$opt\n" if $opt;
					my $values_cmd = $hash_arg_cmd{$arg}.' '."$opt" if ($opt && $hash_arg_cmd{$arg});
					$build_cmd = $build_cmd."$values_cmd " if $values_cmd;
				}
				
				elsif ($type_arg ~~ /^1_main_arguments$/) {
					#print "$opt\n" if $opt;
					$build_cmd = $build_cmd."$opt " if ($opt and $arg ~~ /(reference_in)|(ebwt_base)/) ;
					if ($opt and $arg ~~ /bowtie-path/) {
						$opt =~ s/\/$//;
						$build_cmd = $opt.'/'.$build_cmd;
					}
				}
			}
			
			$count += 1;
		}
	}
}
print $build_cmd."\n";
!system $build_cmd or die "bowtie-build fails: $!";

