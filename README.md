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

**Optional**: compare assigned vs unassigned    
`grep -cE '^@' 04_samples/NGSLib1_ali_assi.fq`    


### Retain only unique reads
`obiuniq -m sample ./04_samples/NGSLib1_ali_assi.fq > 04_samples/NGSLib1_ali_assi_uniq.fa`    

Optional: sum up the count value to make sure all reads are accounted for:    
`grep -E '^>' 04_samples/NGSLib1_ali_assi_uniq.fa | awk -F'count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | paste -sd+ - | bc`

Optional: look at the distribution of counts   
`grep -E '^>' 04_samples/NGSLib1_ali_assi_uniq.fa | awk -F'count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | sort -nr | less`

### Denoise (remove artefactual reads)    
Optional: take another look at the distribution    
`obistat -c count 04_samples/NGSLib1_ali_assi_uniq_trim.fa | sort -nk1 | head -20`

Remove any sequences with fewer than 10 reads, and with a 55 < length < 75
`obigrep --lmin 55 --lmax 75 -p 'count>=10' 04_samples/NGSLib1_ali_assi_uniq_trim.fa > 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75.fa`    

For longer reads, use:   
`obigrep --lmin 100 --lmax 300 -p 'count>=10' 04_samples/NGSLib_uniq.fa > 04_samples/NGSLib_uniq_100-300_10.fa`    


Optional: can also test other lengths by streaming into grep: 
`obigrep --lmin 55 --lmax 75 -p 'count>=5' 04_samples/NGSLib1_ali_assi_uniq_trim.fa | grep -cE '^>' - `

### Remove potential PCR/seq errors    
Label tags as either H, I, or S   
`obiclean -s merged_sample -r 0.05 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75.fa > 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_clean.fa`

Filter    
`obigrep -a 'obiclean_status:s|h' 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_clean.fa > 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_cleanHS.fa`

## Export data     
in tabular format per individual    
`obitab --output-seq 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_cleanHS.fa > 04_samples/NGSLib1_cleanHS.txt`   
This txt file will be input into R below.    

## Assign each sequence to a taxon
Obitools approach for this is not currently working    
For MEGAN, output cannot be in outfmt format, but rather just standard output format   

Using standard blastn against the remote nt database   
`blastn -db nt -query 04_samples/NGSLib1_ali_assi_uniq_trim_c10_55-75_cleanHS.fa -out 05_annotated/NGSLib1_cleanHS_hits.txt -remote -num_descriptions 10 -num_alignments 10`    


### Annotate sequences with species   
Launch MEGAN, import blast output and the fasta file used as the blast query  
Apply the following LCA settings:   
`min score 100`    
`max expected 0.00000001`   
`min % ID 97`   
`top % 10`   
`min support % 0 (off)`   
`min support 1`   

Select level of taxonomy to view, possibly use multiple different levels, e.g. species, genus, family       
File -> Export csv      
Choose: `readName_to_taxonName`    

### Connect read counts and annotations    
This will use the R script, suggested currently to run interactively.   
Necessary inputs:   
Output from `MEGAN`, and output from `obitab`   



