Lara Shonkwiler, Aug 3 2016

This 00README file explains the functioning of:

dnaorg_refseq_assign.pl
dnaorg_evaluate_refseq_assign.pl

----------------------------------------------------------------------------------------------------------------

Program:	dnaorg_refseq_assign.pl


Description:	Given a list of RefSeqs, and a list of accessions, uses the accesions to generate ntlists for 
		each RefSeq.
		
		Specifically:
		Uses E-Utilities to generate FASTA files for all sequences
		Generates an HMM library of RefSeqs
		Runs each accession against the HMM library with nhmmscan
		     - This step is run in parallel on the farm for each sequence
		Parses nhmmscan results, and assigns each accession to one RefSeq
		Generates ntlists for each RefSeq


Usage:          perl dnaorg_refseq_assign.pl <refseq_list> <accn_list>


Input:          <refseq_list> - file containing the accession numbers of all the RefSeqs that should be in the 
                                RefSeq HMM library. Each accession number is on a separate line.

		<accn_list>   - file containing the accession numbers of all the sequences that should be sorted
		                into ntlists.

Output:		<out_folder>  - folder named as "<refseq_list>-ntlists" which contains:

				1) <out_folder>.hmm                - the HMM library of RefSeqs
				2) <out_folder>.log  	           - file containing all output printed to the screen
				   			       	     during the run
				3) <out_folder>.list	           - file containing a list and description of all
				   			             output files
				4) <out_folder>.cmd                - file containg all commands executed by the program
				5) <out_folder>.all.refseqs        - copy of the input list of RefSeqs
				6) <out_folder>.all.seqs           - copy of the input list of sequences to be sorted
				7) <out_folder>.matches.info       - file containing information on all the sequence/RefSeq
				   				     matches				
				8) <out_folder>.non-assigned       - file containing a list of sequences which were not 
				   				     placed in any of the provided RefSeqs' ntlists

				For each accession <accn> in <accn_list>:

				1) <out_folder>.<refseq>.tblout    - summary output of nhmmscan results

				For each RefSeq <refseq> in <refseq_list>:

				1) <out_folder>.<refseq>.ntlist    - properly formatted ntlist for that RefSeq
				   				     (ready for use in dnaorg_annotate.pl)
				

				################## IF --keep OPTION IS ENABLED ###############################
				
				1) <out_folder>.h3m		   - RefSeq HMM library binaries
				   <out_folder>.h3i
				   <out_folder>.h3f
				   <out_folder>.h3p

				For each accession <accn> in <accn_list>, if --keep option is enabled:

				1) <out_folder>.<accn>.qsub.csh  - shell script used to send nhmmscan job to the 
				   				   farm for <accn>
				2) <out_folder>.<accn>.qsub.err  - stderr from farm nhmmscan job for <accn>
				3) <out_folder>.<accn>.qsub.out  - stdout from farm nhmmscan job for <accn>
				4) <out_folder>.<accn>.fasta     - FASTA file for <accn>
				5) <out_folder>.<accn>.nhmmscan.out 
				   				   - complete nhmmscan results for <accn>
								   
				For each RefSeq <refseq> in <refseq_list>, if --keep option is enabled:

				1) <out_folder>.<refseq>.fasta   - FASTA file of <refseq>
				2) <out_folder>.<refseq>.stk     - Stockholm alignment file of <refseq>
				3) <out_folder>.<refseq>.hmm	 - HMM of <refseq>

				###############################################################################

Options:
		-h              : prints info about options
		-c              : genome is closed (a.k.a circular)
  		-f              : force; if dir from --dirout exists, overwrite it
  		-v              : be verbose; output commands to stdout as they're run
  		--dirout <s>    : specify output directory as <s>, not <ref accession>
  		--keep          : do not remove intermediate files, keep them all on disk


Checks:         Can handle blank lines in <accn_list> and in <refseq_list>



Assumptions:    1) The program will sort all the accessions provided. It is assumed that the user has previously
		   filtered out patented sequences, or other undesirable sequences.

-----------------------------------------------------------------------------------------------------------------




-----------------------------------------------------------------------------------------------------------------


Program:	dnaorg_evaluate_refseq_assign.pl


Description:	Given a list of RefSeqs, and a list of accessions, uses the accesions to generate ntlists for 
		each RefSeq.
		
		Specifically:
		Uses E-Utilities to generate FASTA files for all sequences
		Generates an HMM library of RefSeqs
		Runs each accession against the HMM library with nhmmscan
		     - This step is run in parallel on the farm for each sequence
		Parses nhmmscan results, and assigns each accession to one RefSeq
		Generates ntlists for each RefSeq


Usage:          perl dnaorg_evaluate-refseq_assign.pl <results> <standards>


Input:          <results>     - Folder generated by a run of dnaorg_refseq_assign.pl

		<standards>   - Folder containing a list of 'gold standard' ntlists to compare to

Output: 	<out_folder>  - folder named as "<results>-evaluation" which contains:

                                1) <out_folder>.log                - file containing all output printed to the screen
                                                                     during the run
                                2) <out_folder>.list               - file containing a list and description of all
                                                                     output files
                                3) <out_folder>.cmd                - file containg all commands executed by the program

				4) <out_folder>.evaluation         -

Options:
		-h              : prints info about options
  		-f              : force; if dir from --dirout exists, overwrite it
  		-v              : be verbose; output commands to stdout as they're run
  		--dirout <s>    : specify output directory as <s>, not <ref accession>
  		--keep          : do not remove intermediate files, keep them all on disk


Checks:         Can handle seque

Assumptions:    1) There should not be slashes on the end of the folder names of the arguments
		2) ntlists in the <standards> file should be named in the following style: <standards>.<RefSeq Accn#>.ntlist
		 

------------------------------------------------------------------------------------------------------------------





------------------------------------------------------------------------------------------------------------------

Sample run:

Instructions:   dnaorg_refseq_assign should be run first, as dnaorg_evaluate_refseq_assign uses output produced by
		dnaorg_refseq_assign. 

		The following is an example of running these scripts with test_list, a file containing 8 arbitrarily 
		chosen Flaviviral sequences, and refseq_sample_list, a file containing the corresponding RefSeqs of 
		the sequences in test_list. These files should be at the same directory level as the perl scripts.

		[]$ is used to denote command line prompts
		All other lines are onscreen output


###########################################
#
Step 1: Run dnaorg_refseq_assign.pl
#
###########################################


[]$  perl dnaorg_refseq_assign.pl refseq_list_test test_list

# dnaorg_refseq_assign.pl :: Given sequences, decides which RefSeq's ntlist to add them to
# dnaorg 0.1 (August 2016 - ?)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Tue Aug 23 12:59:51 2016
#
# list of refseqs:    refseq_list_test
# list of sequences:  test_list       
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Creating RefSeq HMM Library                                                      ... done. [21.0 seconds]
# Running nhmmscan on sequences                                                    ... Your job 7399076 ("test_list-ntlists.AB078950.qsub.csh") has been submitted
Your job 7399077 ("test_list-ntlists.KJ741267.qsub.csh") has been submitted
Your job 7399078 ("test_list-ntlists.HQ380231.qsub.csh") has been submitted
Your job 7399079 ("test_list-ntlists.FJ502995.qsub.csh") has been submitted
Your job 7399080 ("test_list-ntlists.EU879060.qsub.csh") has been submitted
Your job 7399081 ("test_list-ntlists.JQ308189.qsub.csh") has been submitted
Your job 7399082 ("test_list-ntlists.KP688058.qsub.csh") has been submitted
Your job 7399083 ("test_list-ntlists.JQ418633.qsub.csh") has been submitted
done. [64.0 seconds]
# Creating match info file                                                         ... done. [0.2 seconds]
# Creating ntlists and other output files                                          ... done. [0.0 seconds]
#
# Output printed to screen saved in:                                                                          test_list-ntlists.log
# List of executed commands saved in:                                                                         test_list-ntlists.cmd
# List and description of all output files saved in:                                                          test_list-ntlists.list
# Library of HMMs of RefSeqs. Not press'd saved in:                                                           test_list-ntlists.hmm
# Table summarizing nhmmscan results for AB078950 against test_list-ntlists/test_list-ntlists.hmm saved in:   test_list-ntlists.AB078950.tblout
# Table summarizing nhmmscan results for KJ741267 against test_list-ntlists/test_list-ntlists.hmm saved in:   test_list-ntlists.KJ741267.tblout
# Table summarizing nhmmscan results for HQ380231 against test_list-ntlists/test_list-ntlists.hmm saved in:   test_list-ntlists.HQ380231.tblout
# Table summarizing nhmmscan results for FJ502995 against test_list-ntlists/test_list-ntlists.hmm saved in:   test_list-ntlists.FJ502995.tblout
# Table summarizing nhmmscan results for EU879060 against test_list-ntlists/test_list-ntlists.hmm saved in:   test_list-ntlists.EU879060.tblout
# Table summarizing nhmmscan results for JQ308189 against test_list-ntlists/test_list-ntlists.hmm saved in:   test_list-ntlists.JQ308189.tblout
# Table summarizing nhmmscan results for KP688058 against test_list-ntlists/test_list-ntlists.hmm saved in:   test_list-ntlists.KP688058.tblout
# Table summarizing nhmmscan results for JQ418633 against test_list-ntlists/test_list-ntlists.hmm saved in:   test_list-ntlists.JQ418633.tblout
# Table with statistics for each match saved in:                                                              test_list-ntlists.matches.info
# ntlist for NC_001461 saved in:                                                                              test_list-ntlists.NC_001461.ntlist
# ntlist for NC_001564 saved in:                                                                              test_list-ntlists.NC_001564.ntlist
# ntlist for NC_002657 saved in:                                                                              test_list-ntlists.NC_002657.ntlist
# ntlist for NC_008604 saved in:                                                                              test_list-ntlists.NC_008604.ntlist
# ntlist for NC_027819 saved in:                                                                              test_list-ntlists.NC_027819.ntlist
# List of sequences not assigned to a RefSeq saved in:                                                        test_list-ntlists.non-assigned
# List of RefSeqs in the HMM library saved in:                                                                test_list-ntlists.all.refseqs
# List of sequences that were sorted into ntlists saved in:                                                   test_list-ntlists.all.seqs
#
# All output files created in directory ./test_list-ntlists/
#
# CPU time:  00:01:25.63
#            hh:mm:ss
# 
# DNAORG-SUCCESS





 The output printed to the screen is also saved to the file
 test_list-ntlists.log. It explains the three main
 steps the script performs as they are being performed, and then
 outputs a list of some of the output files created by the script.
 For a complete list see test_list-ntlists.list.
 All output files have been created in the subdirectory 'test_list-ntlists/'.
 A particularly important file is the test_list-ntlists.cmd
 file that includes all of the commands executed by the script during
 its execution.

 The most important file created is the test_list-ntlists.matches.info, which 
 summarizes each assignment of a sequence to a RefSeq. For this run, it should 
 look like this:

 *****************************************************************************

########################################################################################################################################
#
# Query:       Accession number of the sequence
# RefSeq:      The RefSeq that the sequence was assigned to
# Bit score:   Bit score of hit to 'RefSeq'
# E-val:       E-val of hit to 'RefSeq'
# Coverage:    The percentage of the query that the hit to 'RefSeq' covers (Hit length)/(Query length)
# Bias:        TODO
# # Hits:      The number of individual hits to this RefSeq (the program combines stats such as bit score and covg. from separate hits)
# H2: RefSeq:  The RefSeq that had the second strongest hit
# Bit Diff:    The amount by which the bit score of the 'RefSeq' hit is greater than that of the 'H2: RefSeq' hit
# Covg. Diff:  The amount by which the coverage of the 'RefSeq' hit is greater than that of the 'H2: RefSeq' hit
# Num. Correct Hits: The amount of times 'Exp. RefSeq' produced a hit
#
########################################################################################################################################


#H Query                RefSeq                  Bit score               E-val           Coverage                Bias            # Hits          H2: RefSeq              Bit Diff                Covg. Diff
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

AB078950                NC_001461               8983.6          0               0.9993500               302.9           2               NC_002657               0.0013815               1.0000000
KJ741267                NC_001564               11714.7         0.0000000               0.9999064               21.8            1               NC_008604               0.1346190               -1.0000000
HQ380231                NC_002657               12794.1         0.0000000               0.9999187               223.1           1               NC_001461               0.0026021               -1.0000000
FJ502995                NC_008604               11307.2         0.0000000               0.9999077               81.1            1               NC_001564               0.0623789               0.0000000
EU879060                NC_008604               9944.8          0.0000000               0.9999077               81.7            1               NC_001564               0.0642244               0.0000000
JQ308189                NC_008604               11269.9         0.0000000               0.9995388               86.2            1               NC_001564               0.0661378               0.0000000
KP688058                NC_027819               11994.0         0.0000000               0.9999086               32.6            1               NC_008604               0.0279759               0.0000000
JQ418633                NC_001461               9884.9          0               0.9995738               290.7           2               NC_002657               0.0016193               1.00000

 *******************************************************************************

###########################################
#
Step 2: Run dnaorg_evaluate_refseq_assign.pl
#
###########################################

[]$ perl dnaorg_evaluate_refseq_assign.pl test_list-ntlists final_ntlists   
# dnaorg_evaluate_refseq_assign.pl :: Outputs info which helps the user evaluate dnaorg_refseq_assign.pl's performance
# dnaorg 0.11 (August 2016 - ?)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Tue Aug 23 13:34:52 2016
#
# results of dnaorg_refseq_assign:  test_list-ntlists
# 'gold standard' ntlists:          final_ntlists    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Output printed to screen saved in:                                   test_list-ntlists-evaluation.log
# List of executed commands saved in:                                  test_list-ntlists-evaluation.cmd
# List and description of all output files saved in:                   test_list-ntlists-evaluation.list
# Table summarizing the success of dnaorg_refseq_assign.pl saved in:   test_list-ntlists-evaluation.evaluation
#
# All output files created in directory ./test_list-ntlists-evaluation/
#
# CPU time:  00:00:00.02
#            hh:mm:ss
# 
# DNAORG-SUCCESS


The most important file to view from a dnaorg_evaluate_refseq_assign.pl run is the .evaluation file, which
should look like this for this sample run:

 ********************************************************************************

########################################################################################################################################
#
# Query:       Accession number of the sequence
# RefSeq:      The RefSeq that the sequence was assigned to
# Exp RefSeq: The RefSeq that this sequence should have been assigned to
# Pass?: Was the sequence assigned to the right RefSeq? P (pass) or F (fail)
# Bit score:   Bit score of hit to 'RefSeq'
# E-val:       E-val of hit to 'RefSeq'
# Coverage:    The percentage of the query that the hit to 'RefSeq' covers (Hit length)/(Query length)
# Bias:        TODO
# # Hits:      The number of individual hits to this RefSeq (the program combines stats such as bit score and covg. from separate hits)
# H2: RefSeq:  The RefSeq that had the second strongest hit
# Bit Diff:    The amount by which the bit score of the 'RefSeq' hit is greater than that of the 'H2: RefSeq' hit
# Covg. Diff:  The amount by which the coverage of the 'RefSeq' hit is greater than that of the 'H2: RefSeq' hit
# Num. Correct Hits: The amount of times 'Exp. RefSeq' produced a hit
#
########################################################################################################################################


#H Query                RefSeq                  Exp RefSeq              Pass?           Bit score               E-val           Coverage                Bias            # Hits          H2: RefSeq              Bit Diff                Covg. Diff
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------

AB078950                NC_001461               NC_001461               P               8983.6          0               0.9993500               302.9           2               NC_002657               0.0013815               1.0000000
KJ741267                NC_001564               NC_001564               P               11714.7         0.0000000               0.9999064               21.8            1               NC_008604               0.1346190               -1.0000000
HQ380231                NC_002657               NC_002657               P               12794.1         0.0000000               0.9999187               223.1           1               NC_001461               0.0026021               -1.0000000
FJ502995                NC_008604               NC_008604               P               11307.2         0.0000000               0.9999077               81.1            1               NC_001564               0.0623789               0.0000000
EU879060                NC_008604               NC_008604               P               9944.8          0.0000000               0.9999077               81.7            1               NC_001564               0.0642244               0.0000000
JQ308189                NC_008604               NC_008604               P               11269.9         0.0000000               0.9995388               86.2            1               NC_001564               0.0661378               0.0000000
KP688058                NC_027819               NC_027819               P               11994.0         0.0000000               0.9999086               32.6            1               NC_008604               0.0279759               0.0000000
JQ418633                NC_001461               NC_001461               P               9884.9          0               0.9995738               290.7           2               NC_002657               0.0016193               1.0000000

 **********************************************************************************


