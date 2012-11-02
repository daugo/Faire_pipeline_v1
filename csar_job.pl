#!/usr/bin/perl
#csar_job.pl
use warnings;
use strict;
use YAML::Tiny;
use File::Basename;
use Bio::SeqIO;

my $file = shift @ARGV;
Csar_Job($file);

sub Csar_Job {
	my $csar_yml_file = shift;
	my $config_file = YAML::Tiny->new;
	$config_file= YAML::Tiny->read($csar_yml_file) or die "Cannot read config file: $!\n";

	#========Get Parameters====================
	my $treatment_bam = $config_file->[0]->{'input'}->{'treatment'};
	my $control_bam = $config_file->[0]->{'input'}->{'control'};
	my $genome_file = $config_file->[0]->{'input'}->{'genome_fasta'};
	
	my $output_name =  $config_file->[0]->{'output'}->{'output_basename'};

	my $norm = $config_file->[0]->{'parameters'}->{'norm'};
	my $fragment_length = $config_file->[0]->{'parameters'}->{'fragment_length'};
	my $min_score = $config_file->[0]->{'parameters'}->{'min_score'};
	my $max_gap = $config_file->[0]->{'parameters'}->{'max_gap'};
	my $test = $config_file->[0]->{'parameters'}->{'test'};
	my $no_permutations = $config_file->[0]->{'parameters'}->{'no_permutations'};
	my $desired_FDR = $config_file->[0]->{'parameters'}->{'desired_FDR'};

	my $no_digits = $config_file->[0]->{'other_options'}->{'no_digits'};
	my $times_ram_call = $config_file->[0]->{'other_options'}->{'times_ram_call'};
	my $times_ram_per = $config_file->[0]->{'other_options'}->{'times_ram_per'};
	#============================================
	my ($fn,$dir,$suf) = fileparse($treatment_bam,qr/\.[^.]*/);

	#==========Get Chr ID_lengths========================
	my @Seq_IDs; #array where the different identifiers of the mutifasta are going to be saved
	my @lengths;
	my %ID_length; #where lengths of each seqID are going to be saved
	my $seqio = Bio::SeqIO->new(-file => "$genome_file", '-format' => 'Fasta'); #read the fasta file


	#foreach sequence in the fasta file  -> get sequence lengths and stored
	while ( my $seq = $seqio -> next_seq) {
	        my $id = $seq->id; #get id (ex. Chr)
	        $id = "'$id'";
	        push (@Seq_IDs,$id); #save id
	        my $length = $seq->length; #get sequence length
	        push (@lengths,$length);
	        $ID_length {$id} = $length;#save length indexed by seqid
	}
	my $length_list = join(', ',@lengths);

	my $chrm_list = join(', ',@Seq_IDs);
	
	open R_SCRIPT,'>',$fn."_".$output_name.".R" or die "Cannot write CSAR R script\n";
	my $CSAR_R_script = "
#=====libraries======
library(Rsamtools)
library(CSAR)
library(rtracklayer)
#==========================================
source('mappedReads2Nhits_mod.R')
#============================================
source('permutatedWinScores_mod.R')
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

test_csar <- '$test'
number_permutations<- $no_permutations
FDR <- $desired_FDR
output_name <- '$output_name'

#============================================

what <- c('qwidth', 'rname', 'strand', 'pos')
#which <- RangesList(Chr1=IRanges(1,length_Chr))
param <- ScanBamParam(what=what,simpleCigar=TRUE,reverseComplement=TRUE)

bam_S <- scanBam(file=sample,param=param)

bam_C <- scanBam(file=control,param=param)

data_frame_baminfo <- function (bam) {
  Nhits<-rep(1,length(bam[[1]]\$pos))
  data_bam <-data.frame(Nhits=Nhits,lengthRead=bam[[1]]\$qwidth,strand=bam[[1]]\$strand,chr=bam[[1]]\$rname,pos=bam[[1]]\$pos)
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

nhitsS\$filenames
nhitsC\$filenames

bam_scores<-ChIPseqScore(control=nhitsC,sample=nhitsS,file='bam_scores',times=times_call,norm=norm,digits=digits,test=test_csar,backg=1)
bam_scores\$filenames

#Calculate regions of read-enrichment
enriched_regions <-sigWin(experiment=bam_scores,t=t,g=g)
export(object=enriched_regions,paste(output_name,'.bed',sep=''))

#Save the read-enrichment scores at each nucleotide position in a .wig file format
score2wig(experiment=bam_scores,file=paste(output_name,'.wig',sep=''),times=times_call)

for (i in 1:number_permutations) {
  permutatedWinScores_mod(nn=i,sample=data_bam_S,control=data_bam_C,fileOutput='bam_permutations',chr=Chr_names,chrL=length_Chr,w=w,norm=norm,backg=1,times=times_per,digits=digits,t=t,g=g)
}

nulldist<-getPermutatedWinScores(file='bam_permutations',nn=1:number_permutations)

getThreshold(winscores=values(enriched_regions)\$score,permutatedScores=nulldist,FDR=FDR)

";
	print R_SCRIPT $CSAR_R_script;
	close R_SCRIPT;

	!system "Rscript $fn\_$output_name.R" or die "Cannot execute CSAR R script";
	return "$output_name.bed";
	

}