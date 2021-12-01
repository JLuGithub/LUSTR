# LUSTR
LUSTR is designed to genotype germline or somatic alleles by their repeat sizes at each user pointed short tandem repeat (STR) locus. LUSTR requires pre-installation of samtools to process .sam/.bam files. LUSTR is written in perl script and can be run directly without installation.

A typical pipeline to run LUSTR:

(0) Prepare the reference genome/sequences/contigs etc. For each pair-end short read sequenced sample to be genotyped, prepare the raw read pairs as .fastq format. In the situation where raw reads were unavailable, users may use LUSTR_Extractor.pl or other existing tools such as samtools to extract raw read pairs from .sam/.bam files.

(1.1) Prepare a file containing basic information of the STRs of interest on the reference. The file can be generated either manually or translated from results processed by existing tools such as TandemRepeatFinder. The file shall contain lines for each STR with following columns separated by space:
	  Column 1: ID of the STR (no space included);
    Column 2: Chromosome or segment name (identical to the names in the reference) where the STR is located;
    Column 3 & 4: A region indicated by starting and ending position to initiate STR searching. Noting that the region can be imperfect and will be used to search for STR seeds and initiate repeat extension;
	  Column 5: Sequence of the repetitive unit of the STR.
	  Sample: STR_sample1 chr23 150 200 ACG
(1.2) Use LUSTR_Finder.pl to process the STR information file from Step (1.1) against the reference, to obtain standarized sequences of the repeat and flanking regions of each STR of interest as well as other information for subsequent steps. The output file will be plain text containing groups of eight lines for each STR, as follows:
	  Line 1: Information of the STR
	  Line 2 & 3: Sequence of the repeat region of the STR (from the positive strand according to the reference)
	  Line 4 & 5: Sequence of the upstream flanking region of the STR (from the positive strand as well)
	  Line 6 & 7: Sequence of the downstream flanking region of the STR (from the positive strand as well)
	  Line 8: A blank line to separate

(2) Use LUSTR_RefCreater.pl to generate fasta format file from the STR sequences obtained in Step (1.2). The fasta file will serve as the reference to remap raw read pairs to both flanking and repeat sequences. The remapping can be done by existing tools such as bowtie or bwa, but an output in .sam or .bam format without location sorting is required. The parameters for the remapping shall work for pair-end reads and be optimized for STRs. For example, the output may want to contain all mapped hits, without excluding unmapped reads, prefer clips to INDELs, and so forth. Giving bwa mem as an example, besides all other default settings, -a is necessary to export all hits, -D 0.01 -P is strongly recommend to have best results for overlapping STRs, and -T may need to be adjusted when read length is shorter than 100nt.

(3) According to the STR sequences obtained in Step (1.2), use LUSTR_Realigner.pl to process the remapped .sam/.bam files obtained in Step (2) to realign the read pairs to each STR locus. The output file will be plain text containing two header lines for general information and one blank header line to separate, followed by lines for each STR. The lines for each STR locus include one line for the information of the STR, one line indicating the pair number being realigned, several lines for the realignment details of each pair, and one blank line to separate. Each realignment line contains the ID of the pair, and paired "{}" for read 1 and read 2 for their realigned lengths to the upstream flank, the repeat, and the downstream flank, separated by commas ",". In case when the pair can be realigned in multiple ways, multiple paired "{}" will be shown, separated by spaces.

(4) Use LUSTR_Caller.pl to process the realignment results obtained in Step (3). The output file will be plain text containing a brief format instruction, followed by estimations of each STR. For each STR it will contain a header line for the information, followed by lines of all estimated alleles showing size, reliability, and allele fraction. For allele type, "D" indicates a size shorter than read length and can be directly called, "I" indicates a size longer than read length estimated by pairs overlapping with the flanking sequences, "I+" indicates an alternative estimation of "I" but includes pairs with only repeat sequences. An warning log file will also be generated to indicate detection of evidences to potential offtargets, mutations close to STR boundaries, and duplicates filteration. Noting that if the library is composed by reads with largely varied lengths, a warning message will be given and the results can be affected.

Usage:

(0) Extractor
	  Instruction:
		    The Extractor module extracts raw read pairs from given .sam/.bam files.
		    It can process both sorted or unsorted .sam/.bam files.
		    It outputs a single .fastq file written by read pairs.
		    In case the .sam/.bam is incomplete due to filteration of unmapped reads, a warning message will be given.
	  Usage Sample:
		    perl LUSTR_Extractor.pl -i <file_name.bam> -o <file_name.fastq>
    Options:
		    --help/-h		          Print usage instructions
	    	--input/-i <file>	    Input .sam/.bam file
		    --output/-o <file>	  Output .fastq file of raw read pairs
		    --max/-m <value>	    Maximum pair number of intermediate files, which will be automatically deleted when finish
				                      This setting affects memory usage. Setting of 10,000,000 usually costs 5~10G memory
					                    (Default = 10000000)

(1) Finder
	  Instruction:
		    The Finder module processes the basic information of the STRs provided by user, and retrieves standarized STR sequences by user defined parameters.
		    The parameters and other information retrieved will be inherited in subsequent steps.
		    See Step (1.1) and (1.2) for more details.
	  Usage Sample:
		    perl LUSTR_Finder.pl -r <reference.fa> --refindex <reference.fa.fai> -i <STRinfo.txt> -o <STRsequences.txt> -e 100
	  Options:
		    --help/-h		          Print usage instructions
		    --ref/-r <file> 		  The reference file to search for STR sequences
		    --refindex <file>	    Index .fai file of the reference that can be generated by samtools
					                    If not appointed, LUSTR will search for file with ".fai" extension after the file name of the reference
		    --input/-i <file>	    The file for basic information of STRs. See format details in Step (1.1)
		    --output/-o <file>	  Output file for STR sequences. See format details in Step (1.2)
		    --match/-m <value>	  Match score for repeat extension (Default = 2)
		    --mis/-x <value>  	  Mismatch penalty for repeat extension (Default = -5)
		    --gap/-g <value>	    Gap penalty for repeat extension (Default = -7)
		    --stop/-s <value> 	  Threshold to stop repeat extension (Default = -30)
		    --dynamic/-d	    	  Switch on dynamic adjustment of match score for long repetitive units (Default = off) 
		    --comple/-c     		  Switch on searching for the STR sequence on reference negative strand when positive strand searching fails (Default = off)
		    --extra/-e <value>	  This setting helps avoid long time running for certain STR cases when searching for extra repeat extensions allowing mismatches
					                    (Default = no limit, but a setting = 100 is recommended)
		    --flank/-f <value>	  Length to export flanking sequences
					                    This length determines whether the reads in pair are close enough to be considered as from the same STR locus
					                    Please make the best setting accordingly to the sequencing library preparation procedure
					                    (Default = 1000)

(2) RefCreater
	  Instruction:
		    The RefCreater module processes the repeat and flanking sequences of STRs to generate both real and artificial references to remap raw reads to.
		    It outputs a fasta format file with 60nt lines. Sequences with too many Ns can be filtered according to user setting.
	  Usage Sample:
		    perl LUSTR_RefCreater.pl -i <STRsequences.txt> -o <STRreferences.fa>
	  Options:
		    --help/-h	         	  Print usage instructions
		    --input/-i <file> 	  The STR sequences obtained by the Finder module. See format details in Step (1.2)
		    --output/-o <file>	  Output file in fasta format. If not appointed, LUSTR will generate a file with ".fa" extension after the input file name
		    --length/-l <value>	  Minimum length when generating the references for repeats. A setting with 50nt longer than read length is recommended
					                    (Default = 200)
		    --ratio/-r <value>	  A value between 0 and 1 as the threshold for the ratio of Ns. Sequences with Ns beyond this ratio will be discarded
					                    (Default = 0.3) 

(3) Realigner
	  Instruction:
		    The Realigner module realigns the remapped pairs to each STR locus, and generates a plain text file indicating how each pair is realigned.
		    The running speed is affected by both library size and number of STRs. For heavy duty work, the Realigner module allows multi-processing mode.
		    See Step (1.2) and (2) for input requirements, and Step (3) for more details of the output.
	  Usage Sample (single process, read length = ~ 150bp):
		    perl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt>
	  Usage Sample (single process, read length = ~ 100bp):
		    perl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt> --repeatm 0.15 --repeate 0.4 -d 100
	  Usage Sample (multiple parallel processes):
		    perl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt> -m 3 -n 1
		    perl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt> -m 3 -n 2
		    perl LUSTR_Realigner.pl -s <STRsequences.txt> -i <remapped.bam> -o <realign.txt> -m 3 -n 3
	  Options:
		    --help/-h		          Print usage instructions
		    --str/-s <file>	  	  The STR sequences obtained by the Finder module. See format details in Step (1.2)
		    --input/-i <file>	    Unsorted .sam/.bam obtained in the remapping step. See details in Step (2)
		    --output/-o <file>	  Output file for realignment results. See format details in Step (3)
		    --fr/--rr/--rf/--ff	  Indicate how reads are paired if presented (forward-reversely/reverse-reversely/reverse-forwardly/forward-forwardly)
					                    (Default = --fr)
		    --flank/-f <value>	  Maximum length to search for flankings
					                    A setting of the minimum mapping length used in previous step plus 20nt or more is recommended
					                    By default bwa mem requires ~30nt as the minimum mapping length to export the hits
					                    (Default = 50)
		    --seed <value>  		  Seed length to search for flankings. A setting of 7 is recommended for best performance (Default = 7)
		    --untrust1 <value>	  Length threshold of short flankings to enter further check by looking into nearby repeat sequences (Default = 15)
		    --reliab/-r <value>	  Reliability score threshold (0~1) to decide whether to remark short flankings as repeats
					                    Effective for --untrust1 and remark only when acceptable
					                    (Default = 0.8)
		    --untrust2 <value>	  Length threshold of short flankings to enter further check by looking into other supportive pairs
					                    Unsupported short flankings may be remarked as repeats only when acceptable
					                    (Default = 3)
	    	--untrust3 <value>	  Length threshold to forcely remark single-side short flankings as repeats (Default = 0)
		    --clipo <value>	  	  Length threshold to ignore outwards (away from repeats) short clips on the flankings (Default = 5)
		    --clipi <value> 		  Length threshold to ignore inwards (towards the repeats) short clips on the flankings
					                    For bwa mem default settings, 10 is recommended
				                  	  (Default = 10)  
		    --repeatm <value>	    Threshold of repeat ratio to read length (0~1) to initiate extension from remapped hits
					                    A setting = 0.75-2*(minimum mapping length)/readlength is recommended
					                    (Default = 0.35)
		    --repeate <value>	    Threshold of repeat ratio to read length (0~1) to accept the repeats
					                    A setting = 1-2*(minimum mapping length)/readlength is recommended
					                    (Default = 0.6)
		    --back/-b <value>	    Length to step back from remapped ends of reads hitting on repeat references, to search for flankings (Default = 5)
		    --up/-u <value>		    Maximum length of periodic Smith-Waterman extension of repeats towards 5' end of reads
					                    A setting of the minimum mapping length used in previous step plus 20nt is recommended
					                    (Default = 50)
		    --down/-d <value>	    Maximum length of periodic Smith-Waterman extension of repeats towards 3' end of reads
					                    A setting of read length is recommended considering sequencing errors
					                    (Default = 150)
		    --kratio <value>	    Threshold of repeat ratio to read length (0~1) to initiate search for potential repeat units by k-mer method
					                    The k-mer search will cost big running time. Apply it when the STR list is small and a deep realignment is required
                              To apply, a setting = 0.8-2*(minimum mapping length)/readlength is recommended
					                    (Default = 0.4, but will be disabled if neither --kratio nor --kcandi is present)
		    --kcandi <value>	    Besides the best potential repeat unit for each size, consider more repeat units as candidates
					                    For each unit size, candidates with frequencies close to the best within a range determined by this setting will be considered
					                    (Default = 0, but will be disabled if neither --kratio nor --kcandi is present)
		    --multi/-m <value>	  Total number of processes to run parallel (Default = 1)
		    --num/-n <value>	    The order of this run in the parallel processes (Default = 1, but will be disabled if --multi/-m is not present)

(4) Caller
	  Instruction:
		    The Caller module processes the realignment results, and genotypes all possible alleles at each STR.
		    See Step (3) for input requirements, and Step (4) for more details of the output.
	  Usage Sample (without sub output of candidate STRs):
		    perl LUSTR_Caller.pl -i <realign.txt> -o <result.txt> --offtarget
	  Usage Sample (with sub output of candidate STRs):
		    perl LUSTR_Caller.pl -i <realign.txt> -o <result.txt> --offtarget -s 100n --subthres FH
	  Options:
		    --help/-h		          Print usage instructions
		    --input/-i <file>	    Realignment results obtained in the previous step. See details in Step (3)
		    --output/-o <file>	  Output file for STR genotyping results. See format details in Step (4)
		    --untrust1 <value>	  Length threshold of short flankings in direct size calling to enter further check (Default = 3)
		    --release1 <value>	  Reads number threshold to trust short flankings in direct size calling
					                    Effective for --untrust1, and will be automatically modified when repeat size is too long to apply long flankings
					                    (Default = 3)			
		    --untrust2 <value>	  Length threshold of short flankings from reads with single flanking to enter further check (Default = 10)
		    --release2 <value>	  Pair number threshold to trust short flankings from reads with single flanking
					                    If untrusted, the short flankings are considered potential noise from repeat-only pairs
					                    (Default = 1)			
		    --num/-n <value>	    Allowance of reads number variation for fraction calculation and allele possibility check (Default = 3)
		    --ratio/-r <value>	  Ratio variation allowance for allele possibility check (Default = 0.3)
		    --reliabdir1 <value>	Allowance of position variation to determine reliability for repeats shorter than read length (Default = 15)
		    --reliabdir2 <value>	Setting of flanking ratio threshold to determine reliability for repeats shorter than read length (Default = 7)
		    --reliabdir3 <value>	Setting of flanking length threshold to determine reliability for repeats shorter than read length (Default = 10)
		    --reliabest1 <value>	Setting of the shorter flanking length threshold to determine reliability for repeats longer than read length (Default = 20)
		    --reliabest2 <value>	Setting of the longer flanking length threshold to determine reliability for repeats longer than read length (Default = 40)
		    --reliabest3 <value>	Setting of the reads number threshold to determine reliability for alleles longer than read length (Default = 2)
	    	--biasn <value>	    	Reads number threshold for flank bias detection (Default = 20)
		    --biasr <value>		    Calculation ratio threshold (0~1, no bias to complete bias) for flank bias detection (Default = 0.8)
		    --offtarget	        	Switch on single side analysis if flank bias detected (Default = off, but recommended for STRs with high risk of offtargets)
		    --sub/-s <string>	    Switch on a sub output for STR candidates with variations beyond the provided threshold
					                    Use a number or number+"r" to set threshold for repeat number variation, or number+"n" for nucleotide length variation
					                    (Default = off, recommended setting = 100n or 100r if a sub output is wanted)
		    --subthres <string>	  Modification of the sub output effective for --sub/-s
					                    The string is a combination by "Cnumber"+"Apercent"+"F/R"+"H/M/L"+"Y/N"+"E/B" in any order
						                      "Cnumber" - skip candidate STRs not meeting this minimum coverage from non-repeat-only pairs (Default = C10)
						                      "Apercent" - skip candidate STRs when the varied allele has a fraction lower than this minimum percent % (Default = A5)
						                      "F" - only consider varied alleles estimated by non-repeat-only pairs, when "R" consider all estimations (Default = R)
						                      "H"/"M"/"L" - minimum reliability requirement ("High"/"Medium"/"Low") for the alleles with variations (Default = M)
						                      "N" - skip alleles with repeat length close to read length, when "Y" will keep them (Default = N)  
						                      "E" - only consider expanded alleles compared to reference, when "B" consider both expansion/contraction (Default = B)
					                    (Default = C10A5RMNB, but will be disabled if --sub/-s is not present)
