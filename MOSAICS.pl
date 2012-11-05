#!/usr/bin/perl
#MOSAICS.pl
use warnings;
use strict;
use File::Basename;
use YAML::Tiny;
use File::Copy;

#=============
#Parameters

my $mosaics_yml_file = shift @ARGV;
my $config_file = YAML::Tiny->new;
$config_file= YAML::Tiny->read($mosaics_yml_file) or die "Cannot read config file: $!\n";


#CONSTRUCT_BIN PARAMETERS
my $bed_sample = $config_file->[0]->{'input'}->{'treatment'};
my $bed_control = $config_file->[0]->{'input'}->{'control'};
my $out_file_location = $ENV{'PWD'};
my $frag_len = $config_file->[0]->{'construct_bin_parameters'}->{'fragment_length'};
my $bin_size = $config_file->[0]->{'construct_bin_parameters'}->{'bin_size'};
my $capping = $config_file->[0]->{'construct_bin_parameters'}->{'capping'};

#check fragment length, bin size and read length are numbers
unless (($frag_len && $bin_size) ~~ /^[0-9]+$/) {
	die "average_read_length, fragment length and bin size must be integers";
}

#parse filename, dir path and suffix from the input files

my ($fn_S,$dir_S,$suf_S) = fileparse($bed_sample,qr/\.[^.]*/);
my ($fn_C,$dir_C,$suf_C) = fileparse($bed_control,qr/\.[^.]*/);

#========================

#BIN_DATA PARAMETERS
my $bin_suffix = "_fragL$frag_len\_bin$bin_size.txt";
my $bin_file_S = $fn_S.$suf_S.$bin_suffix;
my $bin_file_C = $fn_C.$suf_C.$bin_suffix;
my $M_file = $config_file->[0]->{'bin_files'}->{'mappability_file'};
my $GC_file = $config_file->[0]->{'bin_files'}->{'gc_file'};
my $N_file = $config_file->[0]->{'bin_files'}->{'n_file'};
my $data_type = $config_file->[0]->{'other_parameters'}->{'data_type'};

my $parallel_flag = 'FALSE';
$parallel_flag = 'TRUE' if ($config_file->[0]->{'other_parameters'}->{'parallel_flag'} and $config_file->[0]->{'other_parameters'}->{'parallel_flag'} !~ /false/i) ;
my $cores = $config_file->[0]->{'other_parameters'}->{'cores'};
#=======================
#FIT PARAMETERS
my $analysis_type = $config_file->[0]->{'analysis'}->{'analysis_type'};
#PEAK CALL PARAMETERS
my $FDR = $config_file->[0]->{'analysis'}->{'desired_fdr'};
my $maxgap = $config_file->[0]->{'analysis'}->{'maxgap'};
my $minsize = $config_file->[0]->{'analysis'}->{'minsize'};
my $thres = $config_file->[0]->{'analysis'}->{'thr_score'};

#======================


#== Write Rscript for constructBins ====
construct_bin_Rscript($bed_sample,$fn_S);
construct_bin_Rscript($bed_control,$fn_C);

open R_Script, '>', "mosaics_analysis_$fn_S\_$fn_C\_fragLen$frag_len\_binsize$bin_size.R" or die "Cannot write mosaics_analysis_$fn_S\_$fn_C\_fragL$frag_len\_binsize$bin_size.R: $!";
print R_Script 'library(mosaics)'."\n".'library(multicore)'."\n";
my $Bin_Data_cmd = 
"BinData<-readBins(type=c('chip','input','M','GC','N'), fileName=c('$bin_file_S','$bin_file_C','$M_file','$GC_file','$N_file'), dataType= '$data_type', parallel=$parallel_flag, nCore= $cores)\nBinData\n";
my $plot_cmds = "
pdf(file= '$fn_S\_$frag_len\_$bin_size.pdf')
plot(BinData)
plot(BinData,plotType='input')
plot(BinData, plotType='M')
plot(BinData, plotType='GC')
plot( BinData, plotType='M|input' )
plot( BinData, plotType='GC|input' )\n";

#TS have to be able to change
my $fit_cmd = "
Fit <- mosaicsFit(object=BinData,analysisType='$analysis_type',bgEst='automatic')\nFit\n";

my $plot_fit = "
plot(Fit)\ndev.off()\n";
my $peak_cmd = "
Peak <- mosaicsPeak( Fit, signalModel='2S', FDR=$FDR, maxgap=$maxgap, minsize=$minsize, thres=$thres)\nPeak\n";
my $export_cmd = "
export(Peak, type='txt', filename= '$fn_S\_$frag_len\_$bin_size.txt')
export(Peak, type='bed', filename='$fn_S\_$frag_len\_$bin_size.bed')
export(Peak, type='gff', filename='$fn_S\_$frag_len\_$bin_size.gff')";

print R_Script $Bin_Data_cmd.$plot_cmds.$fit_cmd.$plot_fit.$peak_cmd.$export_cmd;
close R_Script;
!system 'Rscript',"mosaics_analysis_$fn_S\_$fn_C\_fragLen$frag_len\_binsize$bin_size.R" or die "Cannot run Rscript mosaics_analysis_$fn_S\_$fn_C\_fragLen$frag_len\_binsize$bin_size.R: $!";

sub construct_bin_Rscript {
	my $file = shift;
	my $fn = shift;
	open R_Script,'>',"constructBins_$fn.R" or die "Cannot write constructBins_$fn_S.R : $!";
	print R_Script 'library(mosaics)'. "\n";
	my $contruct_bins_cmd = 
	"constructBins(infile='$file',fileFormat='bed',outfileLoc='$out_file_location', byChr=FALSE, fragLen=$frag_len, binSize=$bin_size,capping=$capping)";
	print R_Script "$contruct_bins_cmd";
	close R_Script;
	!system 'Rscript', "constructBins_$fn.R" or die "Cannot run Rscript: $!";
}
