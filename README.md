# eDNA_metabarcoding
Note: This repo is mainly for the developers purpose, although it is attempted to be robust to different data types, it comes with no guarantees of functionality or usefullness.    


Dependencies:    
`python 2.7`    
`gcc`     
`python-dev packages`        
`OBITools` http://metabarcoding.org/obitools/doc/welcome.html       
`MEGAN 6 (Community Edition)` ttps://ab.inf.uni-tuebingen.de/software/megan6     


To make obitools available everywhere, add the obitools binary and the obitools `/export/bin` folder to your $PATH      


Launch Obitools    
`obitools`    

### Prepare raw data
Copy raw data into `02_raw_data`    

Decompress fastq.gz files
```
cd 02_raw_data     
for i in $(ls *.fastq.gz ) ; do gunzip -c $i > ${i%.gz} ; done`
```

### Prepare the interpretation files
Note: This must be done for each sequencing lane separately       

Importantly, for automated matching of interpretation to fastq files, use the following naming convention:      
Interpretation files should be named using all of the input fastq file name up to the `R1_001.fastq` or `R2_001.fastq` part, which is replaced by `interp.txt`    
e.g. fastq files are: `Lib1_S1_L001_R1_001.fastq` and `Lib1_S1_L001_R2_001.fastq`     
and then your interpretation file will be `Lib1_S1_L001_interp.txt`

Use the file `00_archive/header.txt` as a template and create an interpretation file for your samples. An example interpretation is given in `00_archive/interp_example.txt`       
### Merge paired-end reads (PE only)   
This step is for PE reads only, if data is SE, skip to [Separate Individuals](#separate-individuals-se-start).    

Use automated script to run illuminapairedend      
`01_scripts/01_read_merging.sh` 
Essentially, `illuminapairedend --score-min=40 -r R2.fq R1.fq > output.fq`

Then only retain aligned sequences    
`01_scripts/02_retain_aligned.sh`     
Essentially, `obigrep -p 'mode!="joined"' input.fq > output.fq`   

Detect how many reads remain after keeping only merged    
`grep -cE '^\+$' 03_merged/*ali.fq`

### Separate Individuals (SE start)   
Use `ngsfilter`, and your interpretation files, to identify individuals in your aligned fastq. All unidentified samples will go into an unidentified.fq   
`./01_scripts/03_ngsfilter.sh`    

Essentially runs:  `ngsfilter -t your_interp.txt -u unidenfied.fq input.fq > output_ali_assi.fq`    

Optional: count how many reads were assigned from your samples?   
`for i in $(ls 04_samples/*assi.fq) ; do echo $i ; grep -cE '^\+$' $i ;  done`   


NOTE: Depending on your preferences, this would be a good place to concatenate all of your output files together. To do this, run:    
```
mkdir 04_samples/sep_indiv
mv 04_samples/*.fq 04_samples/sep_indiv
cat 04_samples/sep_indiv/*_ali_assi.fq > 04_samples/all_files_ali_assi.fq
```

Once you've combined all sample files into one, the following one-liner will give you a file with the number reads assigned per sample:   
`for i in $(grep -vE '^#' 00_archive/*_interp.txt | awk '{ print $2 }' - | uniq) ; do echo "sample_$i" ; grep -E "sample=$i;" 04_samples/all_files_ali_assi.fq | wc -l ; done > ./04_samples/assigned_reads_per_sample.txt`    

Then if you want to perform a few calculations (e.g. in excel, min, max, mean, sd): 
`tr '\n' ',' < 04_samples/assigned_reads_per_sample.txt | sed 's/,sample/\nsample/g' - | grep -vE 'NTC|ExtCnt|Mock' - | awk -F"," '{ print $2 }' - > 04_samples/perform_calcs.txt`    


### Retain only unique reads
Use obiuniq to retain unique reads within each sample.   
`./01_scripts/04_retain_unique.sh`   

Essentially runs: `obiuniq -m sample input.fq > output_uniq.fa`    

Optional: sum up the count value to make sure all reads are accounted for:    
`grep -E '^>' 04_samples/NGSLib1_ali_assi_uniq.fa | awk -F'count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | paste -sd+ - | bc`

Optional: look at the distribution of counts   
`grep -E '^>' 04_samples/NGSLib1_ali_assi_uniq.fa | awk -F'count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | sort -nr | less`

Optional: take another look at the distribution    
`obistat -c count 04_samples/NGSLib1_ali_assi_uniq_trim.fa | sort -nk1 | head -20`

### Denoise (remove artefactual reads)    
Use obigrep to only retain reads within a specified size range and minimum count.    
Edit the following script to set the `LMIN`, `LMAX` and `MIN_READS` variables, then run it.  
`./01_scripts/05_denoise.sh`    

Essentially runs: `obigrep --lmin <min> --lmax <max> -p 'count>= <min.reads>' input.fa > output.fa`   

Personal suggested uses:   
Valentini primers=55-75    
Other longer amplicons=100-300      

Optional: can also test other lengths by streaming into grep: 
`obigrep --lmin 55 --lmax 75 -p 'count>=5' 04_samples/yourfile.fa | grep -cE '^>' - `

### Remove potential PCR/seq errors    
Use obiclean to label tags as either Head (H), In-between (I), or Singletons (S). Using a 'r=0.05' value currently.    
Then within the same script, filter the output to only keep the H or S reads.   
`./01_scripts/06_remove_errors.sh`   
Essentially runs: `obiclean -s merged sample -r 0.05 input.fa > output.fa`    
then: `obigrep -a 'obiclean_status:s|h' input.fa > output.fa`    

The output of this script will be used as an input to the BLAST query below.   

## Export data     
Use obitab to output a tab-delimited text file that will be used as an input to the R Script below.   
`./01_scripts/07_obitab_export.sh`    
Essentially runs: `obitab --output-seq input.fa > output.txt`   


## Assign each sequence to a taxon
Use a blastn to align the H and S fasta file against nt (NCBI remote).   
`blastn -db nt -query 04_samples/your_cleaned_HS.fa -out 05_annotated/your_lib_output.txt -remote -num_descriptions 10 -num_alignments 10`    

Note: for MEGAN, the blast output must be standard output format, not outfmt.  

### Annotate sequences with species   
Launch MEGAN, import the blast output and the fasta file used as the blast query  
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


### Extras
`for i in $(grep -vE '^#' 00_archive/*_interp.txt | awk '{ print $2 }' - | uniq) ; do echo "sample_$i" ; grep -E "sample=$i;" 04_samples/NGSLib1_ali_assi.fq | wc -l ; done`
