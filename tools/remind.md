
Virus Classification Pipeline (version 0.1)
--------------------------------------------

####1. download viral sequneces and taxnomy database from GenBank ftp (ftp://ftp.ncbi.nih.gov/genbank/)

The 1st command only generate download commands, and the 2nd command will execute the download
commands to download all viral sequences and taxonomy database.

```
$ perl viral_DB_prepare.pl -t download > download.sh  
$ bash ./download.sh
```

####2. run viral_DB_prepare.pl script to classify virus

>__Note__
>please add -c 1 if the classification is from previous result, or please remove `-c 1` if a denovo
classification will be performed 

``` 
$ perl viral_DB_prepare.pl -t category -c 1 gbvrl*.gz 1>report.txt 2>&1
```

####3. manually correction

>__Update__
>the correction of virus taxon id will be automatically execute. An new file `add_taxid.txt` will
>record the incorrect taxon id and corresponding correct one. 

#####3.1 manually corret host name

The pipeline will classify some virus by its host name when it is not classified by means of genus. But the host name
in GenBank may exist in abnormal that could not accepted by taxonomy database. So the file __manual_hname_table.txt__
records abnormal host name. It lists one abnormal name per line. Please search the abnormal name using google, 
genbank taxon database, and iciba to infer the correct name, this correct name must be included in genbank taxon 
database (in name.dmp.gz). Then add the correct name behind the abnormal name in below format: 

>abnormal name [tab] correct name [return]

Example 

######3.1.1 open file manual_hname_table.txt, the hosts were named as below:

	grape cultivar 6-23     
	grape cultivar 8612     
	grape cultivar 87-1     
	grape cultivar Atebage  
	grape cultivar Augusta  
	grape cultivar Benifuji 


######3.1.2 the abnormal host name should be 'Vitis vinifera' or 'wine grape' by searching GenBank taxonomy database.

![img04](http://kentnf.github.io/tools/img/vcp_p4.png)


######3.1.3 add the correct host name behind the abnormal name

	grape cultivar 6-23	Vitis vinifera
	grape cultivar 8612	Vitis vinifera
	grape cultivar 87-1	Vitis vinifera
	grape cultivar Atebage	Vitis vinifera
	grape cultivar Augusta	Vitis vinifera
	grape cultivar Benifuji	Vitis vinifera

#####3.2 manually check the virus genus and classification 

After the classification in step2, the pipeline will create some new genus classification information that not presented in web (www.mcb.uct.ac.za/tutorial/ICTV%20Species%20Lists%20by%20host.htm). For example, the becurtovirus GI:169303562 was found in GenBank nt database, and it was categorized into plant virus according to its host is Suger Beet. Then the pipeline will automatically create a rule that becurtovius should be categorized into plant virus (below).

![img05](http://kentnf.github.io/tools/img/vcp_p5.png)


The new genus classification information was stored in file **manual_genus_table.txt**. 
The file format is one genus & classification per line. If the genus has one division, that is pretty good result. If genus
has more than one division, make it to keep one or more according frequency of each division, or searching
background of this genus to make correct decision.

If you are not familiar with virus genus classification system, please skip this step. That means you trust
the genus & classification information from GenBank.

>__Note__  denovo classification for step 3.3 (with -c 1 parameter) will take very long time, due to lot of virus have not been classified after step 3.2. It is better to perform step 3.4 to update two manually correct files, and perform the classification using the updated files (without -c 1 parameter, then run step 3.4). After that, less unclassified virus need to be analysis in step 3.2, and lot time saved. In my first run of vrl_v205, About 8,000 unclassified virus need to be analyzed in step 3.3 after update files.

#####3.3 manually classification using virus description

__World-Search__  
Afer manually correct in step 3.1 and 3.2. Some viral sequences still can not be classified for 
missing host feature and have a unrecognized genus name. But they have almost same description 
with other classified virus. These virus could be classified by comparing their descriptions with 
classified virus. 

For example, the sequence GI:331725 does not have host feature, and the genus name is 
“Small linear single stranded RNA satellites”. It described as “Cucumber mosaic virus satellite RNA”. 
And another sequence GI:1103556 with same description was classified into plant virus. 
It should be plant virus according to the description.

![img06](http://kentnf.github.io/tools/img/vcp_p6.png)

__Blast-Search__
Besides, the un-classified virus will blast against the classified virus. If the blast result shows 
that they share similar sequences, the un-classified could be classifed according to blast result.

>The word-search and blast-search will automatically execute when perform classification. The output
file __manual_desc_table.txt__ will record the word-search and blast-search result.

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

But the file __manual_desc_table.txt__ need to be checked manually in case making incorrect classification.

Just check the 2nd, 3rd, and 4th column. If it does not make sense,
correct them manually. Please also concern the 5th and 6th column, the error usually happens to the record
without 100% freq.

Most of word search result and blast search result are same, less the 100 of them are diff. So
manually corret the different part will saving a lot of time

> **Suggestion method:**

- Add col A behind col 'match score', compare col 4 and 8 in col A using function 'exac', 
then sort col A and mv rows with 'true' value to another table. We think the classification 
result is correct if the word search method and blast method generate same result.
- then sort col 'Div name by blast', remove col with value 'NA-BLAST'. sort the remaing table
with col 'match length of blast' and 'percentage identity'. 


> - A. sort to find different in word search and blast search (col 4 and 8).
> - B. check the blast match length, identify. Lower than 100 match base, 90% identity show low blast clue for classification
> - C. check the word search column 3, 5, 6, and 7.
> - D. assign a correct div to column 4 for the virus


#####3.4 update the manually correct file to classification

After steps 3.1 to 3.3, we need to update the manually correct files. 
First, backup all the txt files in **tools** folder of VirusDetect. 
Then run below command:
 
	$ perl viral_DB_prepare.pl -t patch

Three new update files will be generated in current folder.
We need to check these updated files: 
        
	update_genus_table_v205.txt
	update_hname_table_v205.txt
	update_desc_table_v205.txt

It is better to check these update files again. The updated information is append behind the label line 

	# patch v205**
	*update information will be here

After check, use the three updated files to replace files in **tools** folder of VirusDetect. 
Last, check the  **default_verion** locate in **tools** folder, make sure it same as the update file.
If the default version is **205**, the three update file name should be:

	update_desc_table_v205.txt  
	update_genus_table_v205.txt  
	update_hname_table_v205.txt

After manually correct all errors in classication, please move three patch files into tools folder of VirusDetect (the step remind me again Virus Classification Pipeline is just a tool of VirusDetect). Then do the same thing as before, but I promise it is the last time for you to run it.

####4. run viral_DB_prepare.pl script to category virus again using manually changed file

	$ perl viral_DB_prepare.pl -t category -c 1 gbvr*.gz

####5. unique the classified virus

	$ perl viral_DB_prepare.pl -t unique input_virus.fasta

This is the last step of classication of virus, the unique virus sequences could be used as reference
of VirusDetect program after create index for bwa and blast

####6. extract protein sequences from unique virus

	$ viral_DB_prepare.pl -t extProt -i input_seq gbvrl1.seq.gz gbvrl2.seq.gz ... gbvrlN.seq.gz

Two file will be generate for protein extracting. One is protein sequences, another is id mapping file for nucleotide and protein. 
One virus nucleotide sequence may contain several ORFs. The last step is move the protein sequence and id mapping file into database
folder of virus-detect. Please follow the file name of reference database.

>vrl_plant
>vrl_plant_prot
>vrl_plant_prot_table

####7. manually update the classification (options)

The classification result of this pipline may contain very limited errors. For any specific reason or specific project, user may need to correct or update the classication for only several virus. 
First, just correct the classification information in vrl_genbank.txt. Then run this command:

	$ perl viral_DB_prepare.pl -t refine

The sequence file will be updated according to vrl_genbank.txt, and the updated sequences will be named like 'vrl_Bacteria_all.fasta.new'.

