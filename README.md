# eDNA_metabarcoding
Note: This repo is mainly for the developers purpose, no guarantees of functionality or usefulness.    


Dependencies:    
`python 2.7`    
`gcc`     
`R`     
`python-dev packages`        
`OBITools` http://metabarcoding.org/obitools/doc/welcome.html       
`MEGAN 6 (Community Edition)` https://ab.inf.uni-tuebingen.de/software/megan6     
`blastn` https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download     

To make obitools available everywhere, add the obitools binary and the obitools `/export/bin` folder to your $PATH      


Launch Obitools    
`obitools`    

This pipeline can handle the following:     
* single-end (SE) or paired-end (PE) data      
* demultiplexed or multiplexed     
* multiple amplicons within a single sample file
...depending on the above, there will be different steps to take in this repo.    

See the workflow figure below to show how to start these various analyses (Figure 1).    

![](00_archive/eDNA_metabarcoding_workflow.png)





### Prepare raw data
Copy raw data into `02_raw_data`    

Decompress fastq.gz files
```
cd 02_raw_data     
for i in $(ls *.fastq.gz ) ; do gunzip -c $i > ${i%.gz} ; done
```

Optional: Run FastQC on input files to view quality   
```
mkdir 02_raw_data/fastqc_output    
fastqc -o 02_raw_data/fastqc_output 02_raw_data/*.fastq    
multiqc -o 02_raw_data/fastqc_output/ 02_raw_data/fastqc_output    
```

### 00. Prepare the interpretation files
Note: This must be done for each sequencing lane or chip separately       

*Importantly*, for auto matching of interpretation to fastq files, use the following naming convention:      
Interpretation files, name with the full fastq file name up to `R1_001.fastq` or `R2_001.fastq`, and replace this with `interp.txt`    
e.g.    
fastq files: `Lib1_S1_L001_R1_001.fastq` and `Lib1_S1_L001_R2_001.fastq`     
interp file: `Lib1_S1_L001_interp.txt`

Use the file `00_archive/interp_example.txt` as a template to create an interp file for your samples.
 

### 01. Merge paired-end reads (PE only)   
This step is for PE reads only, if data is SE, skip to [Separate Individuals](#separate-individuals-se-start).    

Launch script to run illuminapairedend to merge read files      
`01_scripts/01a_read_merging.sh`      
Essentially: `illuminapairedend --score-min=40 -r R2.fq R1.fq > output.fq`

Then launch the following to only retain aligned sequences    
`01_scripts/01b_retain_aligned.sh`     
Essentially: `obigrep -p 'mode!="joined"' input.fq > output.fq`   

Detect how many reads remain after keeping only merged    
`grep -cE '^\+$' 03_merged/*ali.fq`

### 02. Separate Individuals (SE start)   
If you are using single-end data, to match the initial PE steps, use the following script (for Mac just use cp, instead of cp -l):   
`cp -l 02_raw_data/your_file_R1_001.fastq 03_merged/your_file_ali.fq`    

Launch script to Use `ngsfilter` with your interp files to separate individuals out of your aligned fq file. All unidentified reads will go into an unidentified.fq   
`./01_scripts/02_ngsfilter.sh`    
Essentially: `ngsfilter -t your_interp.txt -u unidenfied.fq input.fq > output_ali_assi.fq`    

Detect how many reads were assigned to a sample?   
`for i in $(ls 04_samples/*assi.fq) ; do echo $i ; grep -cE '^\+$' $i ;  done`   

Depending on your project, this may be a good place to simplify by concatenating all of your output files together.        
```
mkdir 04_samples/sep_indiv
mv 04_samples/*.fq 04_samples/sep_indiv
cat 04_samples/sep_indiv/*_ali_assi.fq > 04_samples/all_files_ali_assi.fq
```

The following one-liner will provide the number reads assigned per sample using your assigned fastq and interp (here using `all_files_ali_assi.fq` as assigned fastq:   
`for i in $(grep -vE '^#' 00_archive/*_interp.txt | awk '{ print $2 }' - | uniq) ; do echo "sample_$i" ; grep -E "sample=$i;" 04_samples/all_files_ali_assi.fq | wc -l ; done > ./04_samples/assigned_reads_per_sample.txt`    

Clean up this file before bringing to excel to make a few calculations (e.g. in excel, min, max, mean, sd):      
Note, replace the NTC terms with the NTC labels used in your study.    
`tr '\n' ',' < 04_samples/assigned_reads_per_sample.txt | sed 's/,sample/\nsample/g' - | grep -vE 'NTC|ExtCnt|Mock' - | awk -F"," '{ print $2 }' - > 04_samples/perform_calcs.txt`    


### 03. Retain only unique reads
Use obiuniq to keep one record per unique amplicon in the fastq (outputs fasta).   
`./01_scripts/03_retain_unique.sh`   

Essentially: `obiuniq -m sample input.fq > output_uniq.fa`    
One can also add other -m flags, such as `run`, or `pcr_rep`, etc., anything that you may want to summarize over using obitab later.    

Optional: sum up the count value to make sure all reads are accounted for:    
`grep -E '^>' 04_samples/NGSLib1_ali_assi_uniq.fa | awk -F'count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | paste -sd+ - | bc`

Optional: look at the distribution of counts   
`grep -E '^>' 04_samples/NGSLib1_ali_assi_uniq.fa | awk -F'count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | sort -nr | less`

Optional: take another look at the distribution    
`obistat -c count 04_samples/NGSLib1_ali_assi_uniq_trim.fa | sort -nk1 | head -20`

### 04. Denoise (size and count) and remove putative seq/pcr errors
Use obigrep to only retain reads within a specified size range and minimum count. Then use obiclean to only keep the head (H) or singleton (S) amplicons, not the internals (I) (slight deviations from the head). Currently using `r=0.5`    
    
Edit the following script to set the `LMIN`, `LMAX` and `MIN_READS` variables, then run it.  
`./01_scripts/04_denoise_and_remove_err.sh`    

Essentially, altogether this does: 
```
obigrep --lmin <min> --lmax <max> -p 'count>= <min.reads>' input.fa > output.fa     
obiclean -s merged sample -r 0.05 output_obigrep.fa > output_obiclean.fa   
obigrep -a 'obiclean_status:s|h' output_obiclean.fa > output_all.fa
```

Personal suggested uses:   
Valentini primers=55-75    
Other longer amplicons=100-300      

To help determine how many reads make it through each filtering steps, you can use the following commands.    
For the reads in the fasta after size selecting and count filter:   
`grep -E '^>' 04_samples/*_ali_assi_uniq_c10_55-75.fa | awk -F"count=" '{ print $2 }' - | awk -F";" '{ print $1 }' - | paste -sd+ - | bc`     
For the reads in the fasta after only keeping head and singletons:
`grep -E '^>' all_files_ali_assi_uniq_c10_55-75_clean_HS.fa | awk -F"; count=" '{ print $2 }' - | awk -F";" '{ print $1 }' - | paste -sd+ - | bc`    


Optional: can also test other lengths by streaming into grep: 
`obigrep --lmin 55 --lmax 75 -p 'count>=5' 04_samples/yourfile.fa | grep -cE '^>' - `

The final output of this script will be used as an input to the BLAST query below.   
Note: one can try different count parameters, for example, and can run it within the existing directory.    

## 05. Export data     
Use obitab to output a tab-delimited text file that will be used as an input to the R Script below.   
`./01_scripts/05_obitab_export.sh`    
Essentially: `obitab --output-seq input.fa > output.txt`   


## 06. Assign each sequence to a taxon
Use a blastn to align the H and S fasta file against nt (NCBI remote).   
`blastn -db nt -query 04_samples/your_cleaned_HS.fa -out 05_annotated/your_lib_output.txt -remote -num_descriptions 10 -num_alignments 10`    

Note: for MEGAN, the blast output must be standard output format, not outfmt.  

## 07. Annotate sequences with species   
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

## 08. Connect read counts and annotations    
This will use the R script `read_counts_to_annotations.R`, run interactively.   
Necessary inputs:   
Amplicon annotation output from MEGAN, and amplicon read count from `obitab`   

In brief, this will merge these two inputs, attach locations, aggregate different amplicons with same annotation, calculate proportions, save out proportion plots and count/proportion tables.    

Within here, one can apply a low expression filter to remove any counts less than 10.   
