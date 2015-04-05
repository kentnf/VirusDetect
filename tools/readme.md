
Virus Classification Pipeline (version 0.1)
--------------------------------------------

####1. download viral sequneces and taxnomy database from GenBank ftp (ftp://ftp.ncbi.nih.gov/genbank/)

	$ perl viral_DB_prepare.pl -t download > download.sh  
	$ bash ./download.sh

\* the 1st command only generate download commands, and the 2nd command will execute the download
commands to download all viral sequences and taxonomy database.

####2. run viral_DB_prepare.pl script to classify virus
   
	$ perl viral_DB_prepare.pl -t category gbvrl*.gz 1>report.txt 2>&1
   
There will be several virus do not have taxon id in download files, but they may have correct taxon id on 
GenBank website. So I guess the website update more frequently than ftp. To generate accurate classification, we need manually update 
these virus taxon id in the script viral_DB_prepare.pl, sorry for that I did not make it as a automatically process.

   * check warning message: if there is any virus do not have taxon id.
                            manually chagne it in [sub correct_org_taxon_division]

######2.1 the file "report.txt" records these virus without taxon id like this:

	$ grep "WARN" report.txt
	[WARN]organism is not virus: KC244112
	[WARN]organism is not virus: KC244113
	[WARN]organism is not virus: AF050065
	[WARN]organism is not virus: AF101979
	... 

######2.2 search the ID in WARN in genbank to find correct taxon id
        
######2.3 open file viral_DB_prepare.pl, add the correct taxon id to subroutine correct_org_taxon_division
        
######2.3 re-run the viral_DB_prepare.pl until there is no WARN message

	$ perl viral_DB_prepare.pl -t category gbvrl*.gz 1>report.txt 2>&1


####3. generate manually classification file

	$ viral_DB_prepare.pl -t manually
	$ viral_DB_prepare.pl -t patch
   
There is no process in this step, a guide about how to manually change the
     	classification is print on screen

After manually change, please patch the changed files to previous file

####4. run viral_DB_prepare.pl script to category virus again using manually changed file

	$ viral_DB_prepare.pl -t category -c 1 gbvr*.gz
	$ viral_DB_prepare.pl -t manually
	$ viral_DB_prepare.pl -t patch

   * manually change the desc file, then patch it

   === notice ===
   step 3, and 4 could be combined into one. Before generate and apply patched genes table,
   the step to generate desc table using word search method will take very long time. After
   apply patched genus table, only ~8000 unclassified virus need to be search using word search
   and blast search method (about half day)

####5. append the patched files

After manually correct all errors in classication, please move three patch files into tools folder of VirusDetect (the step remind me again Virus Classification Pipeline is just a tool of VirusDetect). Then do the same thing as before, but I promise it is the last time for you to run it.

	$ perl  viral_DB_prepare.pl -t category -c 1 gbvr*.gz

####6. unique the classified virus

	$ perl viral_DB_prepare.pl -t unique input_virus.fasta

This is the last step of classication of virus, the unique virus sequences could be used as reference
of VirusDetect program after create index for bwa and blast

####7. manually update the classification

The classification result of this pipline may contain very limited errors. For any specific reason or specific project, user may need to correct or update the classication for only several virus. 
First, just correct the classification information in vrl_genbank.txt. Then run this command:

	$ perl viral_DB_prepare.pl -t refine

The sequence file will be updated according to vrl_genbank.txt, and the updated sequences will be named like 'vrl_Bacteria_all.fasta.new'.

