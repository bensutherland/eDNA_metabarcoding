# eDNA_metabarcoding
Note: This repo is mainly for the developers purpose, although it is attempted to be robust to different data types, it comes with no guarantees of functionality or usefullness.    


Dependencies:    
`python 2.7`    
`gcc`     
`python-dev packages`        
`OBITools` http://metabarcoding.org/obitools/doc/welcome.html       

To make obitools available everywhere, add the obitools binary and the obitools `/export/bin` folder to your $PATH      


Launch Obitools    
`obitools`    

Decompress fastq.gz files
```
cd 02_raw_data     
for i in $(ls *.fastq.gz ) ; do gunzip -c $i > ${i%.gz} ; done`
```

### Prepare the interpretation files
Note: This must be done for each sequencing lane separately       
Use the file `00_archive/header.txt` as a template and create an interpretation file for your samples. An example interpretation is given in `00_archive/interp_example.txt`       

### Merge paired-end reads (PE only)   
This step is for PE reads only, skip to ngsfilter if data is SE.    

Use automated script to run illuminapairedend      
`01_scripts/01_read_merging.sh` 

Then only retain aligned sequences    
`01_scripts/02_retain_aligned.sh`     

Detect how many reads remain after keeping only merged    
`grep -cE '^\+$' 03_merged/*ali.fq`

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
Save into the folder  `05_annotated`  

### Connect read counts and annotations    
This will use the R script, suggested currently to run interactively.   
Necessary inputs:   
Output from `MEGAN`, and output from `obitab`   



