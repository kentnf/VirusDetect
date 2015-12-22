
Virus Classification Pipeline (version 0.1)
--------------------------------------------

>__Note__ 
This document just describes how __Virus Classification Pipline__ works.

[TOC]

####1. download viral sequneces and taxnomy database from GenBank ftp (ftp://ftp.ncbi.nih.gov/genbank/)

The 1st command only generate download commands, and the 2nd command will execute the download
commands to download all viral sequences and taxonomy database.

```
$ perl viral_DB_prepare.pl -t download > download.sh  
$ bash ./download.sh
```

####2. run viral_DB_prepare.pl script to initial classification

``` 
$ perl viral_DB_prepare.pl -t category gbvrl*.gz 1>report.txt 2>&1
```
The initial classification may contains some errors or unclassified viruses, please correct them according below description


####3. manually correction

#####3.1 manually corret host name

The pipeline will classify some virus by its host name when it is not classified by means of genus. But the host name
in GenBank may present with abnormal description that could not accepted by taxonomy database. So the file **manual_hname_table.txt** contains all abnormal host name that did not accpet by GenBank. It lists one abnormal name per line. Please correct the abnormal name to make them consistent with GenBank.


Example 

######3.1.1 open file manual_hname_table.txt, the hosts were named as below:

>grape cultivar 6-23     
>grape cultivar 8612     
>grape cultivar 87-1     
>grape cultivar Atebage  
>grape cultivar Augusta  
>grape cultivar Benifuji 


######3.1.2 the abnormal host name should be 'Vitis vinifera' or 'wine grape' by searching GenBank taxonomy database.

![img04](http://kentnf.github.io/tools/img/vcp_p4.png)


######3.1.3 add the correct host name behind the abnormal name

>grape cultivar 6-23 [tab] Vitis vinifera
>grape cultivar 8612 [tab] Vitis vinifera
>grape cultivar 87-1 [tab] Vitis vinifera
>grape cultivar Atebage [tab] Vitis vinifera
>grape cultivar Augusta [tab] Vitis vinifera
>grape cultivar Benifuji [tab] Vitis vinifera


#####3.2 manually check the virus genus and classification 

After the classification in step2, the pipeline will create some new genus classification information that not presented in web (www.mcb.uct.ac.za/tutorial/ICTV%20Species%20Lists%20by%20host.htm). For example, the becurtovirus GI:169303562 was found in GenBank nt database, and it was categorized into plant virus according to its host is Suger Beet. Then the pipeline will automatically create a rule that becurtovius should be categorized into plant virus (below).

![img05](http://kentnf.github.io/tools/img/vcp_p5.png)


The new genus classification information was stored in file **manual_genus_table.txt**. 
The file format is one genus & classification per line. If the genus has one division, that is pretty good result. If genus
has more than one division, make it to keep one or more according frequency of each division, or searching
background of this genus to make correct decision.

If you are not familiar with virus genus classification system, please skip this step. That means you trust
the genus & classification information from GenBank.

#####3.3 manually classification using virus description

After correting host name and genus, some viruses still can not be classified for missing host feature and have a unrecognized genus name. But they have almost same description, or show sequence similarity with other classified viruses. These unclassified viruses could be classified by comparing their descriptions and sequences with classified. The below command will do it automatically. The command will compare the description of unclassified virus with classified virus, then borrow the info of classification to the unclassified.

For example, the sequence M18869 does not have host feature, and the genus name is “Small linear single stranded RNA satellites”. But it described as “Cucumber mosaic virus satellite RNA”. Another sequence X86421 with same description was classified into plant virus. Through blast, The M18869 is 100% covered by X86421 with 96% identity. So, it should be classified as plant virus.

![img06](http://kentnf.github.io/tools/img/vcp_p6.png)

```	
$ perl viral_DB_prepare.pl -t category -c 1 gbvrl*.gz 1>report.txt 2>&1
```

The classified result **manual_desc_table.txt** need to be checked manually in case making incorrect classification.

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


> **Suggestion method:**
> - A. sort to find different in word search and blast search (col 4 and 8), only foucus difference in division.
> - B. check the blast match length, identify. Lower than 100 match base, 90% identity should be manually checked.
> - C. check the word search column 3, 5, 6, and 7.
> - D. fill correct division in column 4


#####3.4 update the manually correct file to classification

Next, the manually correct files should be append to previous correction. The parameter **v** indicate the version of GenBank.
```
$ perl viral_DB_prepare.pl -t patch -v 211
```

Three new update files will be generated in current folder.
       
``` 
update_genus_table_v211.txt
update_hname_table_v211.txt
update_desc_table_v211.txt
```

It is better to check these update files again. The updated information is append behind the label line 

```
# patch v211
---- update information will be here ----
```

Copy the three updated files into **tools/genbank_vrl_classification** folder of VirusDetect. 
Check the **current_genbank_verion** locate in **tools/genbank_vrl_classification** folder, make sure it same as the update file.
If the current version is **211**, the three update file name should be:

```
update_desc_table_v211.txt  
update_genus_table_v211.txt  
update_hname_table_v211.txt
```

####4. run viral_DB_prepare.pl script to category virus again

```
$ perl viral_DB_prepare.pl -t category -c 1 gbvr*.gz
```

####5. extract protein sequences

```
$ viral_DB_prepare.pl -t extProt gbvrl1.seq.gz gbvrl2.seq.gz ... gbvrlN.seq.gz
```

Two file will be generate. 

>vrl_genbank_prot: virus protein sequences
>vrl_genbank_tab: virus nucleotide accession and protein accession

####6. remove redundancy (using plant as example)

remove redundancy sequences with 100% sequence similarity
```
$ perl viral_DB_prepare.pl -t unique vrl_Plants_all.fasta -s 100
```

then generate u99 u97 and u95 base on the U100
```
$ perl viral_DB_prepare.pl -t unique -p 64 -s 99 vrl_Plants_u100 
$ perl viral_DB_prepare.pl -t unique -p 64 -s 97 vrl_Plants_u100
$ perl viral_DB_prepare.pl -t unique -p 64 -s 95 vrl_Plants_u100
```

####7. retrieve proteins for each division (using plant u95 as example)

```
$ perl viral_DB_prepare.pl -t genProt vrl_Plants_u95 vrl_genbank_prot vrl_genbank_tab
```
The output protein sequences will be named as: **vrl_Plants_u95_prot**






