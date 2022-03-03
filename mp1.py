#this is the miniproject code

###PROBLEM 1###
'Retrieve the Illumina reads for the resequencing of K-12 project: https://www.ncbi.nlm.nih.gov/sra/SRX5005282. These are single-end Illumina reads. Your code will retrieve this file; i.e., you cannot just retrieve it.'

import os

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
	'''
	prefetch = 'sratoolkit.2.11.2-ubuntu64/bin/prefetch ' + sra_Number
	#create the prefetch command
	os.system(prefetch) #call on shell to run the prefetch command#will download file to /home/mlemenager/sra-output/sra/SRR8185310.sra

	sra_file = "sra-output/sra/SRR8185310.sra"
	print("Generating fastq for: " + sra_Number)
	fastq_dump = "sratoolkit.2.11.2-ubuntu64/bin/fasterq-dump " + sra_file
	os.system(fastq_dump)
	#generate the fastq file above
	'''

###PROBLEM 2###
	'''
	
	spades_command = 'python3 SPAdes-3.15.4-Linux/bin/spades.py  -t 2 -s ' + sra_Number + '.fastq -o ' + path  + '/' + sra_Number +'_assembly/' 
	#the spades command, -t for two bytes, -s for single reads, -o for specified output directory

	log.write(spades_command+ '\n') #write out command to log file
	os.system(spades_command) #execute command in terminal
	'''
###PROBLEM 3###
	sra_assembly_path = path +'/'+sra_Number+'_assembly/' #the sra_assembly path
	fasta_dict = read_fasta(sra_assembly_path+'contigs.fasta') #find the contigs file and read the fasta file
	fasta_dict_keys = fasta_dict.keys() #get the keys of the dict
	contigsAbove1000 = 0 #initialize counter for # of contigs above thousand
	contigList = [] #list of keys for contigs above 1000
	total_bp_contigs_over_1000 = 0 #counter for total bp
	contigdict1000 = {}	

	for i in fasta_dict_keys:
		headerList = i.split('_')
		length = int(headerList[3])
		if length > 1000:
			contigsAbove1000 +=1
			contigList.append(i)
			total_bp_contigs_over_1000 += length
			contigdict1000[i] = fasta_dict[i]

	with open(path+'/1000_contigs.fasta','w') as out:
		for i in contigdict1000.keys():
			out.write( i + '\n')
			out.write(str(contigdict1000[i])+'\n')
	contig_path = path + '/1000_contigs.fasta'			
	
	log.write('There are ' + str(contigsAbove1000) + ' contigs > 1000 in the assembly.\n')

###Problem 4###
	log.write('There are ' + str(total_bp_contigs_over_1000) + ' bp in the assembly.\n')

###Problem 5###
	predicted_cds_file = path + '/predicted_cds.txt'
	'''
	genemark = 'perl '+ parent_path + '/gms2_linux_64/gms2.pl --seq '+ contig_path + ' --genome-type bacteria --format gff --output ' + path + '/genemarkS2.txt --faa ' + predicted_cds_file
	
	os.system(genemark)
	'''
	#write the command to the prompt
	#write results to results folder in the home directory
	#first writes out the seqeunces in nucl. form to genemarkS2.txt
	#second writes out predicted AA seqs to predicted_cds.txt

###PROBLEM 6###

	query_file = predicted_cds_file
	outfmt = '"10 qseqid sseqid pident qcovs"'
	db_file = 'Ecoli.fasta'
	output_file = path + '/predicted_functionality.csv'
	makedb = 'makeblastdb -in ' + db_file  + '-out ' + path + '/Ecoli -title ' + path + '/Ecoli -dbtype prot'
	db_path = path + '/Ecoli'
	os.system(makedb)
	blast_command='blastp -query '+ query_file + ' -db ' + db_path + ' -max_target_seqs 1 -out ' + output_file + ' -outfmt ' + outfmt 
	print(blast_command)
	os.system(blast_command)
	'''
	with open(output_file, 'r') as blast:
		blast_comps = blast.readlines() #blast comparisons
		print(len(blast_comps))
	'''
