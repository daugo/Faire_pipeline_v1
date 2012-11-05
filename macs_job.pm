package macs_job;
#macs_job.pm
use warnings;
use strict;
use YAML::Tiny;
use File::Basename;


#my $file = shift @ARGV;
#Macs_Job($file);

sub Macs_Job {
	my $macs_yml_file = shift;
	my $config_file = YAML::Tiny->new;
	$config_file= YAML::Tiny->read($macs_yml_file) or die "Cannot read config file: $!\n";

	#========Get Parameters====================
	#Input Args
	my $treatment = $config_file->[0]->{'input'}->{'treatment'};
	my $control = $config_file->[0]->{'input'}->{'control'};
	my $format = $config_file->[0]->{'input'}->{'format'};
	#Output Args
	my $output_name = $config_file->[0]->{'output'}->{'output_name'};
	my $wig = $config_file->[0]->{'output'}->{'wig'};
	my $bed_graph= $config_file->[0]->{'output'}->{'bed_graph'};
	#Common Options
	my $pet_dist = $config_file->[0]->{'common_options'}->{'pet_dist'};
	my $gsize = $config_file->[0]->{'common_options'}->{'gsize'};
	my $tag_size = $config_file->[0]->{'common_options'}->{'tag_size'};
	my $band_width = $config_file->[0]->{'common_options'}->{'band_width'};
	my $p_value = $config_file->[0]->{'common_options'}->{'p_value'};
	my $m_fold = $config_file->[0]->{'common_options'}->{'m_fold'};
	#more_model_options
	my $nolambda =$config_file->[0]->{'more_model_options'}->{'nolambda'};
	my $s_local = $config_file->[0]->{'more_model_options'}->{'s_local'};
	my $l_local = $config_file->[0]->{'more_model_options'}->{'l_local'};
	my $off_auto = $config_file->[0]->{'more_model_options'}->{'off_auto'};
	my $no_model = $config_file->[0]->{'more_model_options'}->{'no_model'};
	my $shift_size = $config_file->[0]->{'more_model_options'}->{'shift_size'};
	#aditional_options
	my $kepp_dup = $config_file->[0]->{'more_model_options'}->{'keep_duplicates'};
	my $to_small = $config_file->[0]->{'more_model_options'}->{'to_small'};
	my $space_wig = $config_file->[0]->{'more_model_options'}->{'space_wig'};
	my $call_subpeaks = $config_file->[0]->{'more_model_options'}->{'call_subpeaks'};
	my $diag = $config_file->[0]->{'more_model_options'}->{'diag'};
	my $fe_min = $config_file->[0]->{'more_model_options'}->{'fe_min'};
	my $fe_max = $config_file->[0]->{'more_model_options'}->{'fe_max'};
	my $fe_step = $config_file->[0]->{'more_model_options'}->{'fe_step'};

	my $macs_command = "macs14 ";
	if ($treatment) {
		$macs_command .= "-t $treatment ";
	}
	else {die "Treatment file $treatment is MANDATORY"}
	$macs_command .= "-c $control " if $control;
	$macs_command .= "-f $format " if $format;
	$macs_command .= "-n $output_name " if $output_name;
	$macs_command .= "-w " if ($wig and $wig !~ /false/i);
	$macs_command .= "-B " if ($bed_graph and $bed_graph !~ /false/i);
	$macs_command .= "--petdist $pet_dist " if $pet_dist;
	$macs_command .= "-g $gsize " if $gsize;
	$macs_command .= "-s $tag_size " if $tag_size;
	$macs_command .= "--bw $band_width " if $band_width;
	$macs_command .= "-p $p_value " if $p_value;
	$macs_command .= "-m $m_fold " if $m_fold;
	$macs_command .= "--nolambda " if ($nolambda and $nolambda !~ /false/i);
	$macs_command .= "--slocal $s_local " if $s_local;
	$macs_command .= "--llocal $l_local " if $l_local;
	$macs_command .= "--off-auto " if ($off_auto and $off_auto !~ /false/i);
	$macs_command .= "--nomodel " if ($no_model and $no_model !~ /false/i);
	$macs_command .= "--shiftsize $shift_size " if $shift_size;
	$macs_command .= "--keep-dup $kepp_dup" if $kepp_dup;
	$macs_command .= "--to-small " if ($to_small and $to_small !~ /false/i);
	$macs_command .= "--space $space_wig " if $space_wig;
	$macs_command .= "--call-subpeaks " if ($call_subpeaks and $call_subpeaks !~ /false/i);
	$macs_command .= "--diag " if ($diag and $diag !~ /false/i);
	$macs_command .= "-fe-min $fe_min " if $fe_min;
	$macs_command .= "-fe-max $fe_max " if $fe_max;
	$macs_command .= "-fe-step $fe_step " if $fe_step;

	print "MACS CMD: $macs_command\n";
	!system $macs_command or die "Cannot: error executing macs";
	return "$output_name.bed"
}

1;


