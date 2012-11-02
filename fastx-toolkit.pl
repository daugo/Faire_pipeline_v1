#!/usr/bin/perl
#fastx-toolkit.pl;
use warnings;
use strict;
use File::Basename;
use YAML::Tiny;
use Data::Dumper;
use File::Copy;

my $yaml_file = shift (@ARGV);
unless ($yaml_file and -e $yaml_file and $yaml_file ~~ /.yml/) {
	die "usage:You must provided a valid YAML config_file with .yml extension.
	Please check the config_file example (fastx-toolkit_config_file.yml) and follow the structure:
		Just fill in the fields with your specific data info!!! 
		($!)"; 
}

my $config_file = YAML::Tiny->new;
# Open the config
$config_file= YAML::Tiny->read("$yaml_file");
#Get YAML parameters
#general options
my $fastx_path = $config_file->[0]->{fastx_path};
$fastx_path =~ s/\/$//;
$fastx_path = $fastx_path.'/';
my $in_file = $config_file->[0]->{in_file};
my @in_files;
push ( @in_files,$in_file);
#print "$_\n" for keys %{ $config_file->[0] };
#stats options
my $stats_perform = $config_file->[0]->{stats_graphs}->{perform};
#clipper options
my $clipper_perform = $config_file->[0]->{clipper}->{perform};
my $adapter = $config_file->[0]->{clipper}->{adapter_seq};
my $discard_len = $config_file ->[0]->{clipper}->{discard_length};
my $keep_adp_lenplus = $config_file->[0]->{clipper}->{keep_adapt_lenplus};
my $discard_non_clip= $config_file->[0]->{clipper}->{discard_non_clip};
my $discard_clip = $config_file->[0]->{clipper}->{discard_clip};
my $adapter_only = $config_file->[0]->{clipper}->{adapter_only};
my $keep_n_seq = $config_file->[0]->{clipper}->{keep_n_seq};
my $clipper_save = $config_file->[0]->{clipper}->{save_previous_file};
#my $histogram_perform = $config_file->[0]->{clipper}->{histogram};

#trimmer options
my $trimmer_perform = $config_file->[0]->{trimmer}->{perform};
my $first_base = $config_file->[0]->{trimmer}->{first_base};
my $last_base = $config_file->[0]->{trimmer}->{last_base};
my $trimmer_save = $config_file->[0]->{trimmer}->{save_previous_file};
#renamer options
my $renamer_perform = $config_file->[0]->{renamer}->{perform};
my $renamer_type = $config_file->[0]->{renamer}->{rename_type};
my $renamer_save = $config_file->[0]->{renamer}->{save_previous_file};
#collapser options
my $collapser_perform = $config_file->[0]->{collapser}->{perform};
my $collapser_save = $config_file->[0]->{collapser}->{save_previous_file};
#artifacts_filter options
my $artifacts_filter_perform = $config_file->[0]->{artifacts_filter}->{perform};
my $artifacts_save = $config_file->[0]->{artifacts_filter}->{save_previous_file};
#quality_filter options
my $quality_filter_perform = $config_file->[0]->{quality_filter}->{perform};
my $minimum_quality = $config_file->[0]->{quality_filter}->{minimum_quality};
my $minimum_percentage = $config_file->[0]->{quality_filter}->{minimum_percentage};
my $quality_filter_save = $config_file->[0]->{quality_filter}->{save_previous_file};

#reverse_complement
my $reverse_perform = $config_file->[0]->{reverse_complement}->{perform};
my $reverse_save = $config_file->[0]->{reverse_complement}->{save_previous_file};

#quality_trimmer
my $quality_trimmer_perform = $config_file->[0]->{quality_trimmer}->{perform};
my $quality_trimmer_thr = $config_file->[0]->{quality_trimmer}->{quality_thr};
my $quality_trimmer_min_len = $config_file->[0]->{quality_trimmer}->{min_length};
my $quality_trimmer_save = $config_file->[0]->{quality_trimmer}->{save_previous_file};
#quality_converter
my $quality_converter_perform = $config_file->[0]->{quality_converter}->{perform};
my $quality_converter_scores = $config_file->[0]->{quality_converter}->{scores};
my $quality_converter_save = $config_file->[0]->{quality_converter}->{save_previous_file};
#uncollapser
my $uncollapser_perform  = $config_file->[0]->{uncollapser}->{perform};
my $uncollapser_in_format = $config_file->[0]->{uncollapser}->{tab_format};
my $uncollapser_tab_col = $config_file->[0]->{uncollapser}->{tab_col};
my $uncollapser_save = $config_file->[0]->{uncollapser}->{save_previous_file};
#masker
my $masker_perform = $config_file->[0]->{masker}->{perform};
my $masker_quality_thr = $config_file->[0]->{masker}->{quality_thr};
my $masker_save = $config_file->[0]->{masker}->{save_previous_file};


#=======================================================================================

#ORDER

#print "$_\n" for keys %{ $config_file->[0] };
my %order_tool;
foreach my $tool (keys %{$config_file->[0]}) {
	if ($tool ne 'fastx_path' and $tool ne 'in_file') {
		my $order = $config_file->[0]{$tool}->{order};
		
		#print "$tool => $order\n";
		if ($order =~ /,/) {
			my @order_mul = split (',',$order);
			foreach (@order_mul) {
				$order_tool{$_} = $tool; 
			}
		}
		else {
			$order_tool{$order} = $tool;
		}
	}
}
#print Dumper %order_tool;
#==========================================================================================
foreach my $order (sort {$a<=>$b} keys(%order_tool)) {
	my $tool = $order_tool{$order};
	#print "$tool\n";
	stats_graphs() if $tool ~~ /^stats_graphs$/;
	clipper() if $tool ~~ /^clipper$/;
	trimmer() if $tool ~~ /^trimmer$/;
	renamer() if $tool ~~ /^renamer$/;
	collapser() if $tool ~~ /^collapser$/;
	artifacts() if $tool ~~ /^artifacts_filter$/;
	quality() if $tool ~~ /^quality_filter$/;
	reverse_complement() if $tool ~~ /^reverse_complement$/;
	quality_trimmer() if $tool ~~ /^quality_trimmer$/;
	quality_converter() if $tool ~~ /^quality_converter$/;
	uncollapser() if $tool ~~ /^uncollapser$/;
	masker() if $tool ~~ /^masker$/;

}
#stats_graphs
sub stats_graphs {
	
	if ($stats_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastq_quality_stats starts ...\n";
		!system $fastx_path."fastx_quality_stats -Q33 -i $in_files[-1] -o $filename.quality_stats" 
		or die "fastx_quality_stats fails : $!";
		
		warn "fastq_quality_boxplot_graph.sh starts ...\n";
		!system $fastx_path."fastq_quality_boxplot_graph.sh -i $filename.quality_stats -o $filename\_quality_boxplot_graph.png"
		or die "fastq_quality_boxplot_graph.sh fails: $!";
		
		warn "fastx_nucleotide_distribution_graph.sh starts ...\n";
		!system $fastx_path."fastx_nucleotide_distribution_graph.sh -i $filename.quality_stats -o $filename\_nucleotide_distribution_graph.png"
		or die "fastx_nucleotide_distribution_graph.sh fails: $!";
		
		#print "fastx_nucleotide_distribution_line_graph.sh starts ...\n";
		#!system $fastx_path."fastx_nucleotide_distribution_line_graph.sh -i $filename.quality_stats -o $filename\_nucleotide_distribution_line_graph.png"
		#or die "fastx_nucleotide_distribution_line_graph.sh fails: $!";
	}
}


#clipper
sub clipper {
	#print "@in_files\n";
	if ($clipper_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastx_clipper starts ..\n";
		my $clipper_cmd = $fastx_path."fastx_clipper -Q33";
		$clipper_cmd = $clipper_cmd."-a $adapter " if $adapter;
		$clipper_cmd = $clipper_cmd."-l $discard_len " if $discard_len;
		$clipper_cmd = $clipper_cmd."-d $keep_adp_lenplus " if $keep_adp_lenplus;
		$clipper_cmd = $clipper_cmd."-c " if $discard_non_clip;
		$clipper_cmd = $clipper_cmd."-C " if $discard_clip;
		$clipper_cmd = $clipper_cmd."-k " if $adapter_only;
		$clipper_cmd = $clipper_cmd."-n " if $keep_n_seq;
		$clipper_cmd = $clipper_cmd."-i $in_files[-1] -o $filename\_fastx_clipper$suffix";
		!system $clipper_cmd or die "fastx_clipper fails: $!";
		push (@in_files,"$filename\_fastx_clipper$suffix");
		discard_pre_file(@in_files) unless $clipper_save;
		
=c histogram need GD::Graph perl module to be installed. Too much effort for just get a graph
		
		if($histogram_perform) {
			!system $fastx_path."fastq_to_fasta -Q33 -i $filename\_fastx_clipper$suffix -o $filename\_fastx_clipper.fasta" or die "fastq_to_fasta fails: $!";
			!system $fastx_path."fasta_clipping_histogram.pl $filename\_fastx_clipper.fasta $filename\_fastx_clipper.png" or die "fasta_clipping_histogram.pl fails: $!";
		}
=cut 	
	}
}	

#trimmer
sub trimmer {
	
	if ($trimmer_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastx_trimmer starts ...\n";
		my $trimmer_cmd = $fastx_path."fastx_trimmer -Q33 ";
		$trimmer_cmd = $trimmer_cmd."-f $first_base " if $first_base;
		$trimmer_cmd = $trimmer_cmd."-l $last_base " if $last_base;
		$trimmer_cmd = $trimmer_cmd."-i $in_files[-1] -o $filename\_fastx_trimmer$suffix";
		#print "trimmer_cmd: $trimmer_cmd\n";
		!system $trimmer_cmd or die "fastx_trimmer fails: $!";
		push (@in_files,"$filename\_fastx_trimmer$suffix");
		discard_pre_file(@in_files) unless $trimmer_save;
	
	}
}

#renamer
sub renamer {
	
	if ($renamer_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastx_renamer starts ...\n";
		my $renamer_cmd = $fastx_path."fastx_renamer -Q33 ";
		$renamer_cmd = $renamer_cmd."-n $renamer_type " if $renamer_type;
		$renamer_cmd = $renamer_cmd."-i $in_files[-1] -o $filename\_fastx_rename$suffix";
		!system $renamer_cmd or die "fastx_rename fails: $!";
		push (@in_files,"$filename\_fastx_rename$suffix");
		discard_pre_file(@in_files) unless $renamer_save;
	}
}

#collapser
sub collapser {
	
	if ($collapser_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastx_collapser starts ...\n";
		my $collapser_cmd = $fastx_path."fastx_collapser -Q33 ";
		$collapser_cmd = $collapser_cmd."-i $in_files[-1] -o $filename\_fastx_collapser$suffix";
		!system  $collapser_cmd or die "fastx_collapser fails : $!";
		push (@in_files,"$filename\_fastx_collapser$suffix");
		discard_pre_file(@in_files) unless $collapser_save;
	}
}

#artifacts
sub artifacts {
	
	if ($artifacts_filter_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastx_artifacts starts ...\n";
		my $artifacts_cmd = $fastx_path."fastx_artifacts_filter -Q33 ";
		$artifacts_cmd= $artifacts_cmd."-i $in_files[-1] -o $filename\_fastx_artifacts$suffix";
		#print "$artifacts_cmd\n";
		!system $artifacts_cmd or die "fastx_artifacts_filter fails: $!";
		push (@in_files,"$filename\_fastx_artifacts$suffix");
		discard_pre_file(@in_files) unless $artifacts_save;
	}
}

#quality
sub quality {
	
	if ($quality_filter_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastq_quality_filter starts ...\n";
		my $quality_filter_cmd = $fastx_path."fastq_quality_filter -Q33 ";
		$quality_filter_cmd = $quality_filter_cmd."-q $minimum_quality " if $minimum_quality;
		$quality_filter_cmd = $quality_filter_cmd."-p $minimum_percentage " if $minimum_percentage;
		$quality_filter_cmd = $quality_filter_cmd."-i $in_files[-1] -o $filename\_fastq_quality_filter$suffix";
		!system $quality_filter_cmd or die "fastq_quality_filter fails: $!";
		push(@in_files,"$filename\_fastq_quality_filter$suffix");
		discard_pre_file(@in_files) unless $quality_filter_save;
		
	}
}

#reverse
sub reverse_complement {
	
	if ($reverse_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastx_reverse_complement starts ...\n";
		my $reverse_cmd = $fastx_path."fastx_reverse_complement -Q33 ";
		$reverse_cmd = $reverse_cmd."-i $in_files[-1] -o $filename\_fastx_reverse_complement$suffix";
		!system $reverse_cmd or die "fastx_reverse_complement fails: $!";
		push(@in_files,"$filename\_fastx_reverse_complement$suffix");
		discard_pre_file(@in_files) unless $reverse_save;
	}
}

#quality_trimmer
sub quality_trimmer {
	
	if ($quality_trimmer_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastq_quality_trimmer starts ...\n";
		my $qtrimmer_cmd = $fastx_path."fastq_quality_trimmer -Q33 ";
		$qtrimmer_cmd = $qtrimmer_cmd."-t $quality_trimmer_thr " if $quality_trimmer_thr;
		$qtrimmer_cmd = $qtrimmer_cmd."-l $quality_trimmer_min_len " if $quality_trimmer_min_len;
		$qtrimmer_cmd = $qtrimmer_cmd."-i $in_files[-1] -o $filename\_quality_trimmer$suffix";
		!system $qtrimmer_cmd or die "fastq_quality_trimmer fails: $!";
		push(@in_files,"$filename\_quality_trimmer$suffix");
		discard_pre_file(@in_files) unless $quality_trimmer_save;
	}
}

#quality_converter
sub quality_converter {
	
	if ($quality_converter_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastq_quality_converter starts ...\n";
		my $qconverter_cmd = $fastx_path."fastq_quality_converter -Q33 ";
		#print "quality_converter_scores: $quality_converter_scores\n";
		$qconverter_cmd = $qconverter_cmd."-n " if ($quality_converter_scores ~~ /^numeric$/i);
		$qconverter_cmd = $qconverter_cmd."-a " if ($quality_converter_scores ~~ /^ascii$/i);
		$qconverter_cmd = $qconverter_cmd."-i $in_files[-1] -o $filename\_quality_converter$suffix";
		!system $qconverter_cmd or die "fastq_quality_converter fails: $!";
		push(@in_files,"$filename\_quality_converter$suffix");
		discard_pre_file(@in_files) unless $quality_trimmer_save;
	}
}
#uncollapser
sub uncollapser {
	if ($uncollapser_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastx_uncollapser starts ...\n";
		my $uncollapser_cmd = $fastx_path."fastx_uncollapser -Q33 ";
		$uncollapser_cmd = $uncollapser_cmd."-c $uncollapser_tab_col " if $uncollapser_in_format;
		$uncollapser_cmd = $uncollapser_cmd."-i $in_files[-1] -o $filename\_uncollapser$suffix";
		!system $uncollapser_cmd or die "fastx_uncollapser fails: $!";
		push(@in_files,"$filename\_uncollapser$suffix");
		discard_pre_file(@in_files) unless $collapser_save;
	}
}

#masker
sub masker {
	
	if ($masker_perform) {
		my ($filename,$directories,$suffix) = fileparse($in_files[-1],qr/\.[^.]*/);
		warn "fastq_masker starts ...\n";
		my $masker_cmd = $fastx_path."fastq_masker -Q33 ";
		$masker_cmd = $masker_cmd."-q $masker_quality_thr " if $masker_quality_thr;
		$masker_cmd = $masker_cmd."-i $in_files[-1] -o $filename\_fastq_masker$suffix";
		!system $masker_cmd or die "fastq_masker fails: $!";
		push(@in_files,"$filename\_fastq_masker$suffix");
		discard_pre_file(@in_files) unless $masker_save;
	}
}

sub discard_pre_file {
	my $previous_file = splice (@_,-2,1) if ($_[-2] ne $_[0]);
	#print "previous_file: $previous_file\n";
	unlink $previous_file if $previous_file;
	
}

print "last file = $in_files[-1]\n";
my ($filename,$directories,$suffix) = fileparse($in_file,qr/\.[^.]*/);
my $folder_name = $filename;
unless (-d $folder_name) {
	mkdir "$folder_name" or die "Cannot create $folder_name: $!";
}
print "\n";
foreach ( glob "$filename*") {
	unless (-d $_ or $_ eq $in_files[0]) {
		move($_,"$folder_name/$_") or die "Cannot move $_ output files to $folder_name: $!";
	}
}

open OUT ,'>>','output_files.log' or die "Cannot create output_files.log :$!";
print OUT localtime()."\t$in_file\t$filename\/$in_files[-1]\n";
close OUT;
