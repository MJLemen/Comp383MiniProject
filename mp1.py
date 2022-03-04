#this is the miniproject code

###PROBLEM 1###

import os
#import os to write commands to prompt line

parent_path = os.getcwd() #get the home directory
directory = 'results' #create the filename to write out to
path = os.path.join(parent_path, directory) #get the entire path name

if not os.path.isdir(path): #if results folder does not exist, create it
    os.mkdir(path)

fileName = path + '/miniproject.log' #pathname to miniproject.log file

def read_fasta(fileName): #read in fasta file
	fastaDict= {} #create dictionary to hold sequence names and sequences
	with open(fileName, 'r') as inp: #open file
		line = inp.readline().rstrip() #read in first line 
		key = '' #initialize key and sequence
		sequence = ''
		while line: #run loop while there are still lines
			if ">" in line: #if there is a carrot, is the key
				key = line
                
			else:
				sequence += line #else is the sequence
			line = inp.readline().rstrip() #grab the next line
			if ">" in line: #if the carrot is in the next line,
				fastaDict[key] = sequence #write out key and seq to dict
				sequence = '' #reset sequence
		fastaDict[key] = sequence #write out last key and sequence
	return fastaDict #reutrn the dictionary

with open(fileName, 'w') as log: #open the file

	sra_Number = 'SRR8185310' #the number information that is downloaded
	
	prefetch = 'sratoolkit.2.11.2-ubuntu64/bin/prefetch ' + sra_Number
	#create the prefetch command
	os.system(prefetch) #call on shell to run the prefetch command#will download file to the /sra_output/sra directory in the results folder

	sra_file = "results/sra_output/sra/" + sra_Number + ".sra"
	#give the path the the .sra file
		
	fastq_dump = "sratoolkit.2.11.2-ubuntu64/bin/fasterq-dump " + sra_file + ' -o ' + path + '/' + sra_Number + '.fastq'
	os.system(fastq_dump)
	#generate the fastq file above using the sra file in the results folder
	

###PROBLEM 2###
	
	
	spades_command = 'python3 SPAdes-3.15.4-Linux/bin/spades.py  -t 2 -s ' + path +'/'+ sra_Number + '.fastq -o ' + path  + '/' + sra_Number +'_assembly/' 
	#the spades command, -t for two bytes, -s for single reads, -o for specified output directory

	log.write(spades_command+ '\n') #write out command to log file
	os.system(spades_command) #execute command in terminal
	
###PROBLEM 3###
	sra_assembly_path = path +'/'+sra_Number+'_assembly/' #the sra_assembly path
	fasta_dict = read_fasta(sra_assembly_path+'contigs.fasta') #find the contigs file and read the fasta file
	fasta_dict_keys = fasta_dict.keys() #get the keys of the dict
	contigsAbove1000 = 0 #initialize counter for # of contigs above thousand
	total_bp_contigs_over_1000 = 0 #counter for total bp
	contigdict1000 = {} #dictionary for contigs over 1000 bp	

	for i in fasta_dict_keys: #run through all keys of contig dict
		headerList = i.split('_') #take header b/c has length value and split into list
		length = int(headerList[3]) #4th item is always the length
		if length > 1000: #if greater than 1000
			contigsAbove1000 +=1 #adds to contig counter
			total_bp_contigs_over_1000 += length #adds to total bases counter
			contigdict1000[i] = fasta_dict[i] #adds sequence and key to another dictionary

	with open(path+'/1000_contigs.fasta','w') as out: #write out the larger sequences to a separate file to be used for genemark
		for i in contigdict1000.keys():
			out.write( i + '\n')
			out.write(str(contigdict1000[i])+'\n')
	contig_path = path + '/1000_contigs.fasta' #make file path variable
	
	log.write('There are ' + str(contigsAbove1000) + ' contigs > 1000 in the assembly.\n')
	#write out to log file amount of contigs above 1000
###Problem 4###
	log.write('There are ' + str(total_bp_contigs_over_1000) + ' bp in the assembly.\n')
	#write out to log file total bp for all contigs over 1000 bps
###Problem 5###
	predicted_cds_file = path + '/predicted_cds.txt'
	#create a file path variable for the predicted_cds.txt file

	genemark = 'perl '+ parent_path + '/gms2_linux_64/gms2.pl --seq '+ contig_path + ' --genome-type bacteria --format gff --output ' + path + '/genemarkS2.txt --faa ' + predicted_cds_file
	#genemarks2 command

	os.system(genemark)
	#write the command to the prompt
	#write results to results folder in the home directory
	#first writes out the seqeunces in nucl. form to genemarkS2.txt
	#second writes out predicted AA seqs to predicted_cds.txt

###PROBLEM 6###

	query_file = predicted_cds_file #create query file variable
	outfmt = '"10 qseqid sseqid pident qcovs"' #out format options
	db_file = 'Ecoli.fasta' #database file being used
	output_file = path + '/predicted_functionality.csv' #output file
	makedb = 'makeblastdb -in ' + db_file  + '-out ' + path + '/Ecoli -title ' + path + '/Ecoli -dbtype prot' #make db command
	db_path = path + '/Ecoli' #db path for blast command
	#os.system(makedb) #make db call
	blast_command='blastp -query '+ query_file + ' -db ' + db_path + ' -max_target_seqs 1 -out ' + output_file + ' -outfmt ' + outfmt 
	#os.system(blast_command) #blast command and call
	#blast the contigs file from genemark output against the Ecoli database and output to predicted_functionality.csv file in results folder

###PROBLEM 7###	
	with open(output_file, 'r') as blast:
		blast_cds = blast.readlines() #blast comparisons
		num_cds = len(blast_cds)
	#open the output file from the blast command
	#and get the total amount of cds (each line is a separate CDS)
	refseq_Ecoli_cds = 4140 #total refseq CDS for the SRR number
	discrepancy = 0 #variable to calculate
	if num_cds > refseq_Ecoli_cds: #write out the discrepancy b/w genmark and refseq
		discrepancy = num_cds - refseq_Ecoli_cds
		log.write('GeneMarkS found ' + str(discrepancy) + ' additional CDS than the RefSeq.\n')
	elif num_cds < refseq_Ecoli_cds:
		discrepancy = refseq_Ecoli_cds - num_cds
		log.write('RefSeq found ' + str(discrepancy) + ' additional CDS than GeneMarkS.\n')
	else:
		log.write('GeneMarkS found the same number of CDS (' + refseq_Ecoli_cds + ') as RefSeq.\n')
	
