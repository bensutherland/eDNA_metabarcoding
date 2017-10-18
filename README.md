# eDNA_metabarcoding
Dependencies:
`python 2.7, gcc, python-dev packages`    
`OBITools` http://metabarcoding.org/obitools/doc/welcome.html       

Add both the obitools binary and the /export/bin folder to the $PATH to make obitools accessible everywhere

Launch Obitools
`obitools`    

Move to 02_raw_data and decompress    
`for i in $(ls *.fastq.gz ) ; do gunzip -c $i > ${i%.gz} ; done`

### Prepare the interpretation files
Note: This must be done for each sequencing lane separately       


This code helps with the multiple indexed libraries in the same file issue (e.g. SOG data)    
`awk -F'\t' '$1=="1" { print "SOG", $2, $5, $6, $7, "F @ NGSlib=1" } ' OFS='\t' ./interp_file_2017-10-17.txt | sed 's/\ //g' > interp_lib1_headless.txt`    
`awk -F'\t' '$1=="2" { print "SOG", $2, $5, $6, $7, "F @ NGSlib=2" } ' OFS='\t' ./interp_file_2017-10-17.txt | sed 's/\ //g' > interp_lib2_headless.txt`    
`awk -F'\t' '$1=="3" { print "SOG", $2, $5, $6, $7, "F @ NGSlib=3" } ' OFS='\t' ./interp_file_2017-10-17.txt | sed 's/\ //g' > interp_lib3_headless.txt`    

Then add header to each
`for i in *headless* ; do cat header.txt $i > ${i%_headless.txt}.txt ; done`

### Merge paired-end reads (PE only)   
Only run this step if reads are paired end. If not, skip ahead to ngsfilter.   

Recover full seq reads from F and R reads
`illuminapairedend --score-min=40 -r 02_raw_data/NGSLib1_S1_L001_R2_001.fastq 02_raw_data/NGSLib1_S1_L001_R1_001.fastq > ./03_merged/NGSLib1.fq`

Only keep aligned sequences    
`obigrep -p 'mode!="joined"' 03_merged/NGSLib1.fq > 03_merged/NGSLib1_ali.fq`    

After running the grep above, you can check how many of your reads were retained using wc   

### Separate Individuals (SE start)   
Use the interpretation files described above to separate your individuals   
`ngsfilter -t ./00_archive/interp_lib1.txt -u 04_samples/unidentified_lib1.fq 03_merged/NGSLib1_ali.fq > ./04_samples/NGSLib1_ali_assi.fq`    

Compare assigned vs unassigned    
`grep -cE '^@' 04_samples/NGSLib1_ali_assi.fq`    


### Retain only unique reads
`obiuniq -m sample ./04_samples/NGSLib1_ali_assi.fq > 04_samples/NGSLib1_ali_assi_uniq.fa`    

Sum up the count value to make sure all reads are accounted for:    
`grep -E '^>' 04_samples/NGSLib1_ali_assi_uniq.fa | awk -F'count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | paste -sd+ - | bc`

Look at the distribution of counts   
`grep -E '^>' 04_samples/NGSLib1_ali_assi_uniq.fa | awk -F'count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | sort -nr | less`

### Denoise (remove artefactual reads)    
Take another look at the distribution    
`obistat -c count 04_samples/NGSLib1_ali_assi_uniq_trim.fa | sort -nk1 | head -20`

Remove any sequences with fewer than 10 reads, and with a 55 < length < 75
`obigrep --lmin 55 --lmax 75 -p 'count>=10' 04_samples/NGSLib1_ali_assi_uniq_trim.fa > 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75.fa`    

Can also test others by streaming into grep: 
`obigrep --lmin 55 --lmax 75 -p 'count>=5' 04_samples/NGSLib1_ali_assi_uniq_trim.fa | grep -cE '^>' - `

### Remove potential PCR/seq errors    
First label tags as either H, I, or S   
`obiclean -s merged_sample -r 0.05 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75.fa > 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_clean.fa`

And filter    
`obigrep -a 'obiclean_status:s|h' 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_clean.fa > 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_cleanHS.fa`

## Assign each sequence to a taxon
Using obitools (not currently working)    
`ecotag -d ~/taxonomy_databases/wolf_tutorial/embl_r117 -R ~/taxonomy_databases/wolf_tutorial/db_v05_r117.fasta ./04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_cleanHS.fa > 05_annotated/NGSLib1_ali_assi_uniq_trim_c10_55-75_cleanHS_annot.fa`    

Using standard BLAST
`blastn -query 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_cleanHS.fa -db ~/blastplus_databases/nt -evalue 1e-3 -outfmt 6 -max_target_seqs 1 > 05_annotated/cleanHS_blastn_results.txt`   


### 


