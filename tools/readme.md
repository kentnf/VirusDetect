
Virus Classification Pipeline (version 0.1)
--------------------------------------------

####1. download viral sequneces and taxnomy database from GenBank ftp (ftp://ftp.ncbi.nih.gov/genbank/)

	$ perl viral_DB_prepare.pl -t download > download.sh  
	$ bash ./download.sh

The 1st command only generate download commands, and the 2nd command will execute the download
commands to download all viral sequences and taxonomy database.

####2. run viral_DB_prepare.pl script to classify virus
   
	$ perl viral_DB_prepare.pl -t category gbvrl*.gz 1>report.txt 2>&1

####3. manually correction

#####3.1 manually correct the virus taxon id

There will be several virus do not have correct taxon id in download files, but they may have correct taxon id on 
GenBank website. So I guess the website update more frequently than ftp. To generate accurate classification, we need manually update 
these virus taxon id in the script viral_DB_prepare.pl, sorry for that I did not make it as a automatically process.

######3.1.1 the file "report.txt" records these virus without taxon id like this:

	$ grep "WARN" report.txt
	[WARN]organism is not virus: KC244112
	[WARN]organism is not virus: KC244113
	[WARN]organism is not virus: AF050065
	[WARN]organism is not virus: AF101979
	... 

The warnings message means virus sequence, take KC244112 as a example, does not have correct taxon id.

######3.1.2 search the locus ID in genbank to find correct taxon id

In GenBank website, the KC244112 has correct taxon id 1600283.

![img01](http://kentnf.github.io/tools/img/vcp_p1.png)


Obviously, the taxon id 1600283 belong to Influenza B virus

![img02](http://kentnf.github.io/tools/img/vcp_p2.png)


######3.1.3 open file viral_DB_prepare.pl, add the correct taxon id to subroutine correct_org_taxon_division

Add the correct virus locus id and taxon id to subroutine correct_org_taxon_division as below:
       
![img03](http://kentnf.github.io/tools/img/vcp_p3.png)
 
######3.1.4 re-run the viral_DB_prepare.pl until there is no WARN message.

	$ perl viral_DB_prepare.pl -t category gbvrl*.gz 1>report.txt 2>&1

#####3.2 manually corret host name

The pipeline will classify some virus by its host name when it is not classified by means of genus. But the host name
in GenBank may exist in abnormal that could not accepted by taxonomy database. The step 2 will generate a file **manual_hname_table.txt** which records abnormal host name.

The file manual_hname_table.txt lists one abnormal name per line. Please search the abnormal name using google, 
genbank taxon database, and iciba to infer the correct name, this correct name must be include in genbank taxon 
database (in name.dmp.gz). Then add the correct name behind the abnormal name in below format: 
>abnormal name [tab] correct name [return]

>example
>![img03](http://kentnf.github.io/tools/img/vcp_p3.png)
>A. open file manual_hname_table.txt
>B. search the abnormal name using google and genbank
>C. add the correct to behind the abnormal name

#####3.3 manually check the virus genus and classification 

manual_genus_table.txt

The file format is one genus & classification per line. If the genus has one division, that is pretty good result. If genus
has more than one division, make it to keep one or more according frequency of each division, or searching
background of this genus to make correct decision.

If you are not familiar with virus genus classification system, please skip this step. That means you trust
the genus & classification information from GenBank.

#####3.4 manually classification using virus description

manual_desc_table.txt

>The format of manual_desc_table.txt:  
>1 - ID: GFLRNA3  
>2 - Desc: Grapevine fanleaf virus satellite RNA (RNA3), complete cds.  
>3 - Same Desc: Grapevine fanleaf virus  
>4 - Plants  
>5 - No. of same: 1  
>6 - Freq of same: 100.00  
>7 - ID of same: GFLRNA1  
>8 - Div name by blast  
>9 - best hit of blast  
>10- match length of blast  
>11- percentage identity  
>12- match score  

The script seach the same description of unclassified virus against classified virus, then brorrow the info
of classification to the unclassified. Just check the 2nd, 3rd, and 4th column. If it does not make sense,
correct them manually. Please also concern the 5th and 6th column, the error usually happens to the record
without 100% freq.

Beside correction with word search method, we also blast the unclassified virus against the classified virus
(col 8-12). Most of word search result and blast search result are same, less the 100 of them are diff. So
manually corret the different part will saving a lot of time

> Steps
-A. sort to find different in word search and blast search (col 4 and 8).
-B. check the blast match length, identify. Lower than 100 match base, 90% identity show low blast clue for classification
-C. check the word search column 3, 5, 6, and 7.
-D. assign a correct div to column 4 for the virus

After correction, please save the correction in same file with same format.
then patch the changed files to previous file using:
 
	$ perl viral_DB_prepare.pl -t patch


        check the updated file in current folder
        1. update_genus_table_v205.txt
                old :
                        Alphapermutotetravirus  NA      NA
                update: # patch v205
                        Alphapermutotetravirus  G1      Invertebrates
                Then remove the line on old file
        2. update_hname_table_v205.txt
                The manually changed name should be in the end of this file

#####3.5 appned the manually correct file to classification

	$ perl viral_DB_prepare.pl -t patch
   
DESCAfter manually change, please patch the changed files to previous file

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

