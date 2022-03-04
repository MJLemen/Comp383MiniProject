# Comp383MiniProject
The miniproject is a python wrapper using a combination of SPAdes, GeneMark, and Blast.

GOAL: The goal of the project was to look at the strain of K-12 E. coli and determine whether the strain had mutated and developed new genes by comparing the refseq number of genes and the genemark amount of coding sequences identified.

CODE EXPLANATION:
This code takes in an SRA number (format is SRR#), grabs the .sra file using the prefetch command in SRA tools, and then generates the .fastq file using the faster-q command. Then, the fastq file has the genome assebled using spades to creae a contigs.fasta file in SR#####_assembly directory. The file is then read in and sorted for all contigs over length of 1000 bases. Those contigs are then wrote out to a separate .fasta file. The .fasta file is then fed into genemark to determine the genes of the contigs that code for amino acid sequences, and writes that out to the predicted_cds.txt file. An ecoli database is then generated using an ecoli databise file from ncbi. Blast is then used to compare the predicted amino acid sequences in prerdicted_cds.txt and provides the top gene hit for each of the sequences and outputs to the file predicted_functionality.csv. Each line is for each of the predicted genes present in the assembled genome. Then, the number of predicted genes is compared against the number of genes present in the RefSeq genome to determine if RefSeq or GeneMarkS2 found more genes.

HOW THE CODE WAS TESTED:
The code was tested using the SRA number: SRR8185310 for the E. coli strain K-12

The Ecoli dataset used a multifasta file of genes from Ecoli from Prokka, which is located in the github directory under Ecoli.fasta.


PACKAGE DOCUMENTATION:
All the following Packages need to be installed in the home directory from which the code was executed. (Ex. if running the mp1.py in the /home/user_name directory, the package would need to be installed in /home/user_name directory as well.

PACKAGE: NCBI SRA TOOLKIT
Version:  v3.0.0
Directory: sratoolkit.2.11.2-ubuntu64 
Guide to Installation: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
code to install: 'wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz'

What SRA ToolKit Was Used For:
prefetch command to retrieve the .sra file, faster-q command to generate the .fastq file, Blastp to compare the fasta files to the db, makeblastdb to create the ecoli database

Configuration for SRA TOOLKIT:
While configuring the package (command is vdb-config -i), make sure the cache output is being sent to the results folder in the home directory you are running the code! (EX: home/wheree_running_code/results)

PACKAGE: SPAdes genome assembler
Version: SPAdes genome assembler v3.15.4
website link: https://cab.spbu.ru/files/release3.15.4/manual.html

Code to install: 
First: 'wget http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Linux.tar.gz' Second:'tar -xzf SPAdes-3.15.4-Linux.tar.gz'

What SPAdes genome assembler was used for:
This package was used to assemble the genome of the SRA number given at  the beginning of the code.

PACKAGE: GENEMARKS2
Version: GeneMarkS-2 version 1.14_1.25
website link: 'http://exon.gatech.edu/GeneMark/license_download.cgi'

What Genemarks2 was used for:
Predict the coding sequences of contigs and translate into amino acids.

PACKAGE: BLAST+
version: v2.2.31
installation link: 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/'

What Blast+ was used for:
Blast+ was used for the makeblastdb to create the ecoli database with the multi-fasta file from Prokka.
Blast+ was also used for Blastp to compare the contigs coding regions against the database to determine the genes in the assembled genome.
