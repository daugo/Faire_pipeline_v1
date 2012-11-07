-- MySQL dump 10.13  Distrib 5.1.63, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: GFF3_store
-- ------------------------------------------------------
-- Server version	5.1.63-0ubuntu0.11.10.1

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `BED_Peaks_Info`
--

DROP TABLE IF EXISTS `BED_Peaks_Info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `BED_Peaks_Info` (
  `Peak_ID` int(11) NOT NULL AUTO_INCREMENT,
  `file_path` varchar(255) NOT NULL,
  `filename` varchar(255) NOT NULL,
  `suffix` char(6) NOT NULL,
  `chr_seqID` char(20) NOT NULL,
  `start` bigint(20) NOT NULL,
  `end` bigint(20) NOT NULL,
  `program` char(20) NOT NULL,
  `score` float DEFAULT NULL,
  `Calculate_Summit` bigint(20) DEFAULT NULL,
  PRIMARY KEY (`Peak_ID`),
  UNIQUE KEY `unique_same` (`chr_seqID`,`Calculate_Summit`,`filename`),
  CONSTRAINT `BED_Peaks_Info_ibfk_1` FOREIGN KEY (`chr_seqID`) REFERENCES `SEQID_LANDMARK` (`Name`)
) ENGINE=InnoDB AUTO_INCREMENT=43359 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Closest_Analysis`
--

DROP TABLE IF EXISTS `Closest_Analysis`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Closest_Analysis` (
  `Closest_ID` int(11) NOT NULL AUTO_INCREMENT,
  `Peak_ID` int(11) NOT NULL,
  `RNA_ID` varchar(50) NOT NULL,
  `Distance` int(11) NOT NULL,
  `Feature` varchar(50) NOT NULL,
  PRIMARY KEY (`Closest_ID`),
  UNIQUE KEY `closest_unique_entry` (`Peak_ID`,`RNA_ID`,`Distance`,`Feature`),
  KEY `RNA_ID` (`RNA_ID`),
  CONSTRAINT `Closest_Analysis_ibfk_2` FOREIGN KEY (`RNA_ID`) REFERENCES `Parent_RNA` (`RNA_ID`),
  CONSTRAINT `Closest_Analysis_ibfk_1` FOREIGN KEY (`Peak_ID`) REFERENCES `BED_Peaks_Info` (`Peak_ID`)
) ENGINE=InnoDB AUTO_INCREMENT=78265 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Features`
--

DROP TABLE IF EXISTS `Features`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Features` (
  `ID` int(11) NOT NULL AUTO_INCREMENT,
  `seqid` varchar(20) NOT NULL,
  `Source` varchar(20) NOT NULL,
  `Feature` varchar(50) NOT NULL,
  `Start` bigint(20) NOT NULL,
  `End` bigint(20) NOT NULL,
  `Strand` enum('+','-','?','.') NOT NULL,
  `Phase` enum('0','1','2') DEFAULT NULL,
  `Description_GFF` varchar(255) DEFAULT NULL,
  `Parent_ID` varchar(50) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE KEY `seqid` (`seqid`,`Feature`,`Start`,`End`,`Parent_ID`),
  KEY `Parent_ID` (`Parent_ID`),
  CONSTRAINT `Features_ibfk_1` FOREIGN KEY (`seqid`) REFERENCES `SEQID_LANDMARK` (`Name`),
  CONSTRAINT `Features_ibfk_2` FOREIGN KEY (`Parent_ID`) REFERENCES `Parent_RNA` (`RNA_ID`)
) ENGINE=InnoDB AUTO_INCREMENT=479599 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `GO_Annotation`
--

DROP TABLE IF EXISTS `GO_Annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `GO_Annotation` (
  `GO_Ann_ID` int(11) NOT NULL AUTO_INCREMENT,
  `locus_name` varchar(50) NOT NULL,
  `db_specific_id` varchar(50) DEFAULT NULL,
  `object_name` varchar(50) DEFAULT NULL,
  `relationship_type` varchar(100) DEFAULT NULL,
  `GO_term` varchar(255) DEFAULT NULL,
  `GO_ID` varchar(50) NOT NULL,
  `keyword_ID` varchar(50) DEFAULT NULL,
  `aspect_ID` char(1) DEFAULT NULL,
  `GOslim_term` varchar(255) DEFAULT NULL,
  `evidence_code` char(3) DEFAULT NULL,
  `evidence_description` varchar(255) DEFAULT NULL,
  `evidence_with` varchar(50) DEFAULT NULL,
  `reference` varchar(100) DEFAULT NULL,
  `annotator` varchar(20) DEFAULT NULL,
  `date_annotated` date DEFAULT NULL,
  PRIMARY KEY (`GO_Ann_ID`),
  KEY `locus_name` (`locus_name`),
  KEY `GO_ID` (`GO_ID`),
  CONSTRAINT `GO_Annotation_ibfk_2` FOREIGN KEY (`GO_ID`) REFERENCES `GO_IDs` (`GO_ID`),
  CONSTRAINT `GO_Annotation_ibfk_1` FOREIGN KEY (`locus_name`) REFERENCES `Gene` (`Gene_ID`)
) ENGINE=InnoDB AUTO_INCREMENT=688759 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `GO_Enrichment`
--

DROP TABLE IF EXISTS `GO_Enrichment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `GO_Enrichment` (
  `Rank` int(11) NOT NULL,
  `GO_ID` varchar(50) NOT NULL,
  `Annotated` int(11) NOT NULL,
  `Significant` int(11) NOT NULL,
  `Expected` float DEFAULT NULL,
  `Rank_Classic` int(11) DEFAULT NULL,
  `classic` float DEFAULT NULL,
  `weight` float DEFAULT NULL,
  `list_peaks` char(20) NOT NULL,
  PRIMARY KEY (`GO_ID`,`list_peaks`),
  KEY `list_peaks` (`list_peaks`),
  CONSTRAINT `GO_Enrichment_ibfk_1` FOREIGN KEY (`GO_ID`) REFERENCES `GO_IDs` (`GO_ID`),
  CONSTRAINT `GO_Enrichment_ibfk_2` FOREIGN KEY (`list_peaks`) REFERENCES `Peak_Combinations` (`list_peaks`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `GO_IDs`
--

DROP TABLE IF EXISTS `GO_IDs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `GO_IDs` (
  `GO_ID` varchar(50) NOT NULL,
  PRIMARY KEY (`GO_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Gene`
--

DROP TABLE IF EXISTS `Gene`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Gene` (
  `seqid` varchar(20) NOT NULL,
  `Source` varchar(20) NOT NULL,
  `Type` varchar(50) NOT NULL,
  `Start` bigint(20) NOT NULL,
  `End` bigint(20) NOT NULL,
  `Strand` enum('+','-','?','.') NOT NULL,
  `Gene_ID` varchar(50) NOT NULL,
  `Description_GFF` varchar(255) NOT NULL,
  PRIMARY KEY (`Gene_ID`),
  KEY `seqid` (`seqid`),
  CONSTRAINT `Gene_ibfk_1` FOREIGN KEY (`seqid`) REFERENCES `SEQID_LANDMARK` (`Name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Intersect_Programs_Info`
--

DROP TABLE IF EXISTS `Intersect_Programs_Info`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Intersect_Programs_Info` (
  `chr_seqID` char(20) NOT NULL,
  `start` bigint(20) NOT NULL,
  `end` bigint(20) NOT NULL,
  `num_peaks` int(11) NOT NULL,
  `list_peaks` char(20) NOT NULL,
  `group_index` int(11) NOT NULL,
  `Peak_ID` int(11) NOT NULL,
  PRIMARY KEY (`chr_seqID`,`start`,`Peak_ID`),
  KEY `Peak_ID` (`Peak_ID`),
  KEY `list_peaks` (`list_peaks`),
  CONSTRAINT `Intersect_Programs_Info_ibfk_2` FOREIGN KEY (`list_peaks`) REFERENCES `Peak_Combinations` (`list_peaks`),
  CONSTRAINT `Intersect_Programs_Info_ibfk_1` FOREIGN KEY (`Peak_ID`) REFERENCES `BED_Peaks_Info` (`Peak_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Meme_Mast`
--

DROP TABLE IF EXISTS `Meme_Mast`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Meme_Mast` (
  `Peak_ID` int(11) NOT NULL,
  `p_value` float NOT NULL,
  `e_value` float NOT NULL,
  PRIMARY KEY (`Peak_ID`),
  CONSTRAINT `Meme_Mast_ibfk_1` FOREIGN KEY (`Peak_ID`) REFERENCES `Meme_Sites` (`Peak_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Meme_Sites`
--

DROP TABLE IF EXISTS `Meme_Sites`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Meme_Sites` (
  `motifID` char(20) NOT NULL,
  `seq_meme_id` char(20) NOT NULL,
  `Peak_ID` int(11) NOT NULL,
  `position` int(11) NOT NULL,
  `strand` enum('+','-') NOT NULL,
  `p_value` float NOT NULL,
  `left_flank` varchar(255) DEFAULT NULL,
  `site` varchar(255) NOT NULL,
  `rigth_flank` varchar(255) DEFAULT NULL,
  `filename` varchar(255) NOT NULL,
  PRIMARY KEY (`motifID`,`Peak_ID`),
  KEY `Peak_ID` (`Peak_ID`),
  KEY `motifID` (`motifID`,`filename`),
  CONSTRAINT `Meme_Sites_ibfk_2` FOREIGN KEY (`motifID`, `filename`) REFERENCES `Meme_motifs` (`motif_ID`, `filename`),
  CONSTRAINT `Meme_Sites_ibfk_1` FOREIGN KEY (`Peak_ID`) REFERENCES `BED_Peaks_Info` (`Peak_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Meme_motifs`
--

DROP TABLE IF EXISTS `Meme_motifs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Meme_motifs` (
  `motif_ID` char(20) NOT NULL,
  `width` int(11) NOT NULL,
  `sites` int(11) NOT NULL,
  `ic` float NOT NULL,
  `re` float NOT NULL,
  `llr` int(11) NOT NULL,
  `e_value` float NOT NULL,
  `bayes_thr` float NOT NULL,
  `elapsed_time` float NOT NULL,
  `Reg_ex` text,
  `filename` varchar(255) NOT NULL,
  PRIMARY KEY (`motif_ID`,`filename`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Overlap_Analysis`
--

DROP TABLE IF EXISTS `Overlap_Analysis`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Overlap_Analysis` (
  `Overlap_ID` int(11) NOT NULL AUTO_INCREMENT,
  `Peak_ID` int(11) NOT NULL,
  `RNA_ID` varchar(50) NOT NULL,
  `Feature` varchar(50) NOT NULL,
  `Feature_No` int(11) DEFAULT NULL,
  PRIMARY KEY (`Overlap_ID`),
  KEY `Peak_ID` (`Peak_ID`),
  KEY `RNA_ID` (`RNA_ID`),
  CONSTRAINT `Overlap_Analysis_ibfk_1` FOREIGN KEY (`Peak_ID`) REFERENCES `BED_Peaks_Info` (`Peak_ID`),
  CONSTRAINT `Overlap_Analysis_ibfk_2` FOREIGN KEY (`RNA_ID`) REFERENCES `Parent_RNA` (`RNA_ID`)
) ENGINE=InnoDB AUTO_INCREMENT=14285 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Parent_RNA`
--

DROP TABLE IF EXISTS `Parent_RNA`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Parent_RNA` (
  `seqid` varchar(20) NOT NULL,
  `Source` varchar(20) NOT NULL,
  `Type` varchar(50) NOT NULL,
  `Start` bigint(20) NOT NULL,
  `End` bigint(20) NOT NULL,
  `Strand` enum('+','-','?','.') NOT NULL,
  `RNA_ID` varchar(50) NOT NULL,
  `Description_GFF` varchar(255) NOT NULL,
  `Parent_ID` varchar(50) NOT NULL,
  PRIMARY KEY (`RNA_ID`),
  KEY `seqid` (`seqid`),
  KEY `Parent_ID` (`Parent_ID`),
  CONSTRAINT `Parent_RNA_ibfk_1` FOREIGN KEY (`seqid`) REFERENCES `SEQID_LANDMARK` (`Name`),
  CONSTRAINT `Parent_RNA_ibfk_2` FOREIGN KEY (`Parent_ID`) REFERENCES `Gene` (`Gene_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `Peak_Combinations`
--

DROP TABLE IF EXISTS `Peak_Combinations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `Peak_Combinations` (
  `list_peaks` char(20) NOT NULL,
  PRIMARY KEY (`list_peaks`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `SEQID_LANDMARK`
--

DROP TABLE IF EXISTS `SEQID_LANDMARK`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `SEQID_LANDMARK` (
  `Name` varchar(20) NOT NULL,
  `Length` bigint(20) DEFAULT NULL,
  PRIMARY KEY (`Name`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2012-11-06 12:09:52
