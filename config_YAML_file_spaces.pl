#!/usr/bin/perl
use warnings;
use strict;

my $config_file_name = shift @ARGV;

open CONFIG,'<',$config_file_name or die "Cannot open $config_file_name: $!";
my $file =  join('', <CONFIG>);
close CONFIG;
$file =~ s/-\n\s{5}/-/g;
$file =~ s/(\s{6})\s{2}-/$1-/g;
#print $file;
unlink $config_file_name;
open OUT,'>',$config_file_name or die "Cannot create $config_file_name: $!";
print OUT "#YAML\n#$config_file_name\n";
print OUT $file;
close OUT;
