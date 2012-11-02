#!/usr/bin/perl
#macs_job.pl
use warnings;
use strict;
use YAML::Tiny;
use File::Basename;

my $file = shift @ARGV;
Macs_Job($file);

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


sub CSAR_JOB {
	my $csar_yml_file = shift;
	my $config_file = YAML::Tiny->new;
	$config_file= YAML::Tiny->read($csar_yml_file) or die "Cannot read config file: $!\n";

	#========Get Parameters====================
	my $treatment_bam = $config_file->[0]->{'input'}->{''};
	my $control_bam = $config_file->[0]->{'input'}->{''};
	my $genome_file = $config_file->[0]->{'input'}->{''};
	
	my $output_name =  $config_file->[0]->{'output'}->{''};

	my $norm = $config_file->[0]->{'parameters'}->{''};
	my $fragment_length = $config_file->[0]->{'parameters'}->{''};
	my $min_score = $config_file->[0]->{'parameters'}->{''}; 
	my $max_gap = $config_file->[0]->{'parameters'}->{''};
	my $test = $config_file->[0]->{'parameters'}->{''};
	my $no_permutations = $config_file->[0]->{'parameters'}->{''};
	my $desired_FDR = $config_file->[0]->{'parameters'}->{''};

	my $no_digits = $config_file->[0]->{'other_options'}->{''};
	my $times_ram_call = $config_file->[0]->{'other_options'}->{''};
	my $times_ram_per = $config_file->[0]->{'other_options'}->{''};
	#============================================
	
	#==========Get Chr ID_lengths========================
	my @Seq_IDs; #array where the different identifiers of the mutifasta are going to be saved
	my @lengths;
	my %ID_length; #where lengths of each seqID are going to be saved
	my $seqio = Bio::SeqIO->new(-file => "$fasta_file", '-format' => 'Fasta'); #read the fasta file


	#foreach sequence in the fasta file  -> get sequence lengths and stored
	while ( my $seq = $seqio -> next_seq) {
	        my $id = $seq->id; #get id (ex. Chr)
	        push (@Seq_IDs,$id); #save id
	        my $length = $seq->length; #get sequence length
	        push (@lengths,$length);
	        $ID_length {$id} = $length;#save length indexed by seqid
	}
	my $length_list = 
	
	my $chrm_list =
	

	 
	
	

	my $CSAR_R_script = "
#=====libraries======
library(Rsamtools)
library(CSAR)
library(rtracklayer)
#==========================================
source("mappedReads2Nhits_mod.R")
#============================================
source("permutatedWinScores_mod.R")
#===========================================

#============Parameters======================

sample<-'$treatment_bam'
control<-'$control_bam'
length_Chr <- c($length_list) #30427671
norm <- $norm #300*10^8
Chr_names <- c($chrm_list)
w <- $fragment_length
digits <- $no_digits
t<- $min_score
g<-$max_gap
times_call <- $times_ram_call
times_per <- $times_ram_per

test -> $test
number_permutations<- $no_permutations
FDR <- $desired_FDR
output_name <- $output_name

#============================================

what <- c('qwidth', 'rname', 'strand', 'pos')
#which <- RangesList(Chr1=IRanges(1,length_Chr))
param <- ScanBamParam(what=what,simpleCigar=TRUE,reverseComplement=TRUE)

bam_S <- scanBam(file=sample,param=param)

bam_C <- scanBam(file=control,param=param)

data_frame_baminfo <- function (bam) {
  Nhits<-rep(1,length(bam[[1]]$pos))
  data_bam <-data.frame(Nhits=Nhits,lengthRead=bam[[1]]$qwidth,strand=bam[[1]]$strand,chr=bam[[1]]$rname,pos=bam[[1]]$pos)
  return (data_bam)
}

data_bam_S <- data_frame_baminfo(bam_S)
summary(data_bam_S)
data_bam_C <- data_frame_baminfo(bam_C)
summary(data_bam_C)


nhits <- function (data_bam,index) {
 
  nhits<-mappedReads2Nhits_mod(data_bam,file=paste('data_bam_',index),chr=Chr_names,chrL=c(length_Chr))
  return(nhits)
}

nhitsS <- nhits(data_bam_S,'S')
nhitsC <- nhits(data_bam_C,'C')

nhitsS$filenames
nhitsC$filenames

bam_scores<-ChIPseqScore(control=nhitsC,sample=nhitsS,file='bam_scores',times=times_call,norm=norm,digits=digits,test= test,backg=1)
bam_scores$filenames

#Calculate regions of read-enrichment
enriched_regions <-sigWin(experiment=bam_scores,t=t,g=g)
export(object=enriched_regions,paste(output_name,'.bed',sep=''))

#Save the read-enrichment scores at each nucleotide position in a .wig file format
score2wig(experiment=bam_scores,file=paste(output_name,'.wig',sep=''),times=times_call)

for (i in 1:number_permutations) {
  permutatedWinScores_mod(nn=i,sample=data_bam_S,control=data_bam_C,fileOutput='bam_permutations',chr=Chr_names,chrL=length_Chr,w=w,norm=norm,backg=1,times=times_per,digits=digits,t=t,g=g)
}

nulldist<-getPermutatedWinScores(file='bam_permutations',nn=1:number_permutations)

getThreshold(winscores=values(enriched_regions)$score,permutatedScores=nulldist,FDR=FDR)
"


}

sub MOSAICS_JOB {

}

=c
#==========Get Chr ID_lengths========================
my @Seq_IDs; #array where the different identifiers of the mutifasta are going to be saved
my %ID_length; #where lengths of each seqID are going to be saved
my $seqio = Bio::SeqIO->new(-file => "$fasta_file", '-format' => 'Fasta'); #read the fasta file


#foreach sequence in the fasta file  -> get sequence lengths and stored
while ( my $seq = $seqio -> next_seq) {
        my $id = $seq->id; #get id (ex. Chr)
        push (@Seq_IDs,$id); #save id
        my $length = $seq->length; #get sequence length
        $ID_length {$id} = $length;#save length indexed by seqid
}
=cut 
