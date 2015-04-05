
Virus Classification Pipeline (version 0.1)
--------------------------------------------

####1. download viral sequneces and taxnomy database from GenBank ftp (ftp://ftp.ncbi.nih.gov/genbank/)

	$ viral_DB_prepare.pl -t download > download.sh  
	$ bash ./download.sh

\* the 1st command only generate download commands, and the 2nd command will execute the download
commands to download all viral sequences and taxonomy database.

####2. run viral_DB_prepare.pl script to category virus
   
	$ viral_DB_prepare.pl -t category gbvrl*.gz 1>report.txt 2>&1
   
   * There will be several virus do not have taxon id in download files, but they may have correct taxon id on 
     GenBank website (I guess the website update more frequently than ftp). To make 	 we need manually update 
     these virus taxon id in the pipline. 

   * check warning message: if there is any virus do not have taxon id.
                            manually chagne it in [sub correct_org_taxon_division]

   2.1 the file "report.txt" records these virus without taxon id, 
	$ grep "WARN" report.txt

   2.2 search the ID in WARN in genbank to find correct taxon id (img)
        
   2.3 open file viral_DB_prepare.pl, add the correct taxon id to subroutine correct_org_taxon_division
        
   2.3 re-run the viral_DB_prepare.pl until there is no WARN message
	 $ viral_DB_prepare.pl -t category gbvrl*.gz 1>report.txt 2>&1

##3. generate manually classification file
   $ viral_DB_prepare.pl -t manually
   $ viral_DB_prepare.pl -t patch
   * there is no process in this step, a guide about how to manually change the
     classification is print on screen
   * after manually change, please patch the changed files to previous file

##4. run viral_DB_prepare.pl script to category virus again using manually changed file
   $ viral_DB_prepare.pl -t category -c 1 gbvr*.gz
   $ viral_DB_prepare.pl -t manually
   $ viral_DB_prepare.pl -t patch
   * manually change the desc file, then patch it

   === notice ===
   step 3, and 4 could be combined into one. Before generate and apply patched genes table,
   the step to generate desc table using word search method will take very long time. After
   apply patched genus table, only ~8000 unclassified virus need to be search using word search
   and blast search method (about half day)

##5. move the patched file into bin folder, then run classification again
   $ viral_DB_prepare.pl -t category -c 1 gbvr*.gz

6. unique the classified virus
   $ viral_DB_prepare.pl -t unique input_virus.fasta

*7.if the output vrl_genbank.txt has been manually changed
   $ viral_DB_prepare.pl -t refine [no parameters]
   * output is like vrl_Bacteria_all.fasta.new
