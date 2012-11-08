Faire_pipeline_v1
=================

##Overview

The pipeline allows the user to perform a complete analysis for FAIRE-seq datasets. It integrates different analysis tools, starting from common but time consuming NGS data analysis procedures to functional approximations as over-representation of Gene Ontology categories and motif discovery.
Configuration files, in YAML markup language (e.g., MotifAnalysis options.yml), provides user interaction and traceability. The YAML files that users need to configure at each stage are labeled in this [Figure](http://bce-user-004.uniandes.edu.co/FAIRE_ARAB/pipeline.png). An example of each configuration file needed in each stage is provided. Users must provide in the YAML configuration files, information specific for their experiment. Users could also change many of the parameters for the different tools in the pipeline or let the default values unchanged.
￼

The structure of the pipeline, the order in which each process is carried out, the input data that must be supplied, the output format of each process and the configuration files required are summarized in this [Figure](http://bce-user-004.uniandes.edu.co/FAIRE_ARAB/pipeline.png). The pipeline was written following a modular approach. Users could just exploit some parts of the pipeline. The parts that can be independently requested, are: pre-processing of raw reads, mapping to reference genome using Bowtie, processing of alignment files up to BED formatted files, peak-calling, functional analysis (annotation parsing,Peak annotation, Term enrichment of GO and Motif discovery) and intersection analysis between peaks datasets. Advanced users can easily add new available tools and modify the interaction with the ones currently included.
All provided results, from identification of the enriched regions, to functional approximations are stored in a relational database powered by MySQL ([Figure DB](http://bce-user-004.uniandes.edu.co/FAIRE_ARAB/DB_squema.png)). This will aid in downstream analysis, publication of results and design of further experiments. Exploratory graphics are automatically generated using R and can help in the interpretation of results.
The pipeline is written in Perl, thus it can be used with any operative system where the interpreter runs. However, it has only been tested on Linux.

##Configuration Files Descriptions

Templates for each part of the pipeline are provided in ([ConfigurationFiles folder](https://github.com/daugo/Faire_pipeline_v1/tree/master/ConfigurationFiles)). Inside each configuration file template there is a description of the fields and which of them are mandatory. For non-mandatory fields, the indicator '~'  is requiered as  null/undefined value. Comments are supported by the '#' indicator.

##Availability and requirements
* Project name: FAIRE-seq pipeline (I will think in a better name)
* Project home page: github site
* Operating system(s): Platform independent (Tested on Ubuntu Desktop 11.10 64 bit) • Programming language: Perl (v5.12.4)
* Other requirements: Programs: FASTQC 0.10.1, FASTX Toolkit 0.0.13, Bowtie 0.12.8, samtools 0.1.17, BEDTools Suite v2.16.2, MACS (macs14 1.4.1 20110622), R-2.15 with Bioconductor packages (CSAR, MOSAiCS, multicore, Rsamtools, TopGO), MySQL 14.14, MEME-CHIP.
* Perl Additional Modules: YAML::Tiny 1.51, BioPerl. 
* License: GNU GPL
