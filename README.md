# eDNA_metabarcoding
Note: This repo is mainly for the developers purpose, no guarantees of functionality or usefulness.    

Dependencies:    
`OBITools` http://metabarcoding.org/obitools/doc/welcome.html       
`MEGAN 6 (Community Edition)` https://ab.inf.uni-tuebingen.de/software/megan6     
`blastn` https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download     
`cutadapt` http://cutadapt.readthedocs.io/en/stable/     
`R`     

To make obitools available everywhere, add the obitools binary and the obitools `/export/bin` folder to your $PATH     

Launch OBITools    
`obitools`    

This pipeline can handle the following, and to find the appropriate pipeline, see Figure 1:     
* single-end (SE) or paired-end (PE) data      
* demultiplexed or multiplexed data     
* multiple amplicons within a single sample file     

![](00_archive/eDNA_metabarcoding_workflow.png)
**Figure 1.** eDNA_metabarcoding workflow, showing the front end different options for either single-end (SE) or paired-end (PE) data, prior to the main analysis section. The grey box pipelines are variants derived from the standard multiplexed workflow (currently more stable).     

### Prepare Raw Data and Interpretation File
Copy raw data into `02_raw_data`, decompress, then run fastqc to view quality.    

```
cd 02_raw_data     
for i in $(ls *.fastq.gz ) ; do gunzip -c $i > ${i%.gz} ; done
mkdir 02_raw_data/fastqc_output    
fastqc -o 02_raw_data/fastqc_output 02_raw_data/*.fastq    
multiqc -o 02_raw_data/fastqc_output/ 02_raw_data/fastqc_output    
```

The interpretation file must be made for each sequencing lane or chip separately.      
Use `00_archive/interp_example.txt` as a template.            
**Importantly**, name interp file with input fastq name, replacing `R[1/2]_001.fastq` with `interp.txt`    
e.g. `Lib1_S1_L001_R1_001.fastq`, `Lib1_S1_L001_R2_001.fastq`, `Lib1_S1_L001_interp.txt`       


## Part 1A. Enter Pipeline - Multiplexed Data
### 1A.1.a. Merge Paired-End Reads 
Paired-end data will undergo read merging first, run the following script:    
`01_scripts/01a_read_merging.sh`     
(in brief: `illuminapairedend --score-min=40 -r 02_raw_data/*R2.fq 02_raw_data/*R1.fq > 03_merged/*merged.fq`)      
Retain only the merged (aligned) reads:     
`01_scripts/01b_retain_aligned.sh`     
(in brief: `obigrep -p 'mode!="joined"' 03_merged/*merged.fq > 03_merged/*ali.fq`)   

Audit: how many reads remain after keeping only merged?     
`grep -cE '^\+$' 03_merged/*ali.fq`


### 1A.1.b. Mimic PE Step For SE Samples
Single-end data, to catch up w/ paired-end, run the following command:   
`cp -l 02_raw_data/your_file_R1_001.fastq 03_merged/your_file_ali.fq`    


### 1A.2. Separate Individuals   
Use ngsfilter with the interp file(s) to demultiplex samples out of the `*.ali.fq` file(s).     
`./01_scripts/02_ngsfilter.sh`    
(in brief: `ngsfilter -t 00_archive/*_interp.txt -u unidenfied.fq 03_merged/*ali.fq > 04_samples/*_ali_assi.fq`)    

Audit: how many reads were assigned to a sample?   
`for i in $(ls 04_samples/*_ali_assi.fq) ; do echo $i ; grep -cE '^\+$' $i ;  done`   

Each output file now should be annotated with a sample ID in the read accession header, and if so, one can concatenate all files together now, as follows:  
```
mkdir 04_samples/sep_indiv
mv 04_samples/*.fq 04_samples/sep_indiv
cat 04_samples/sep_indiv/*_ali_assi.fq > 04_samples/all_files_ali_assi.fq
```

Move on to [Part 2](#part-2-main-analysis).

## Part 1B. Enter Pipeline - De-Multiplexed Data
This section is the preparation of input data section if your data comes de-multiplexed.    
Depending on the data type (see Figure 1), the steps taken here will vary. See Variant A and Variant B.        

### Variant A. De-multiplexed single-amplicon (SE and PE)
### 1B.0. Cutadapt
As the barcodes are not used to de-multiplex in this case, the primer sequence still must be removed:     
`01_scripts/00_cutadapt_primer_seq.sh`    
(in brief: `cutadapt -g adapt1 -G adapt2 -o 02_raw_data/*R1_001_noprime.fastq -p 02_raw_data/*R2_001_noprime.fastq 02_raw_data/*R1_001.fastq.gz 02_raw_data/*R1_001_fastq.gz`)       

### 1B.1.a. Merge Paired-End Reads
See Part 1A, 01a (merge Paired-End Reads for details, but the following is for those without primers:   
`01_scripts/01a_read_merging_noprime.sh`    
`01_scripts/01b_retain_aligned.sh`    

Because we did not use ngsfilter, we need to annotate each read accession with sample IDs, then the read files can be combined into a single file.    
`01_scripts/obiannotate_ident.sh`      
(in brief: `obiannotate -S sample:$i 04_samples/*.fq > 04_samples/*_sannot.fq`)    
(in brief: `cat *datatype_sannot.fq > 04_samples/datatype_merged_data_assi.fq`)    

Move on to [Part 2](#part-2-main-analysis).

### 1B.1.b. Use ngsfilter to Enter Pipeline with Single-End Reads With Unidentified Reads
Single-end data, enter the obitools pipeline via a failed run of ngsfilter and take all of the 'unidentied reads' per sample as your sample's reads:    
`01_scripts/02_ngsfilter_SE_exp_unident.sh`    
(in brief: `ngsfilter -t $interp -u 04_samples/*_unidentified.fq 02_raw_data/*_001.fastq > 04_samples/*_assi.fq`)
As noted, all the unidentified files will have full data in them, and the assigned are empty.      

As above, we need to annotate each read accession with a sample ID (no ngsfilter success), then combine all files into one:    
`01_scripts/obiannotate_unident.sh`
(in brief: `obiannotate -S sample:$i 04_samples/*L001_Rq_unidentified.fq > 04_samples/*_sannot.fq`)
(in brief: `cat 04_samples/*_sannot.fq > 04_samples/merged_data_assi.fq`)

To improve identification of identical amplicons, cut the SE data to a uniform size using cutadapt:     
`cutadapt --length 230 -o 04b_annotated_samples/merged_data_assi_230.fq 04b_annotated_samples/merged_data_assi.fq`      
Move on to [Part 2](#part-2-main-analysis).

### Variant B. De-multiplexed multiple-amplicon (SE option only)
Single-end data, enter the obitools pipeline by using ngsfilter with the primer sequence split into the first six basepairs as a fake 'barcode' and the last sequence of the primer as the primer sequence. This way you can, per sample, de-multiplex your data by amplicon type.    
`01_scripts/02_ngsfilter.sh`    
(as described above)

Then per sample, the data can be split into the two amplicon types (it is currently just named in the accession):   
`01_scripts/00b_split_by_type.sh`    
Note here that this script will need to be edited for your use. Currently uses grep to take the following three lines after the identifier of interest.    

Do the rest separately for your two amplicon types.    
As above, use obiannotate to annotate your sequences.      
`01_scripts/obiannotate_ident.sh`     

Move on to [Part 2](#part-2-main-analysis).

## Part 2. Main Analysis 
![](00_archive/eDNA_metabarcoding_workflow_2.png)

### 2.1. Retain Only Unique Reads
Input is a single fastq file containing all samples for a specific amplicon, annotated with sample name.   

Use obiuniq to keep one record per unique amplicon in the fastq (outputs fasta).   
For paired-end data: `./01_scripts/03_retain_unique_PE.sh`      
For single-end data: `./01_scripts/03_retain_unique_SE.sh`     
(in brief: `obiuniq -m sample 04_samples/*assi.fq > 04_samples/*_uniq.fa`)        
Note: one can also add other -m flags, such as `run`, or `pcr_rep`, etc., anything that you may want to summarize over using obitab later.    

Audit: sum up the count value to make sure all reads are accounted for:    
`grep -E '^>' 04_samples/your_file_ali_assi_uniq.fa | awk -F'\ count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | paste -sd+ - | bc`

Audit: look at the distribution of counts   
`grep -E '^>' 04_samples/your_file_ali_assi_uniq.fa | awk -F'\ count=' '{ print $2 }' - | awk -F';' '{ print $1 }' | sort -nr | less`

### 2.2. Denoise by Filtering by Size, Count, and PCR/Seq Error
Use obigrep to only retain reads within a specified size range and minimum count. Then use obiclean to only keep the head (H) or singleton (S) amplicons, not the internals (I) (slight deviations from the head). Currently using `r=0.5`    

PE data does full filtering as above:    
`./01_scripts/04_denoise_and_remove_err.sh` (edit LMIN, LMAX, MIN_READS)    

SE data filters only on size and count:    
`./01_scripts/04_denoise_and_remove_err_SE.sh` (edit LMIN, LMAX, MIN_READS)    
    
In brief, this does the following: 
```
obigrep --lmin <min> --lmax <max> -p 'count>= <min.reads>' input.fa > output.fa     
obiclean -s merged sample -r 0.05 output_obigrep.fa > output_obiclean.fa   
obigrep -a 'obiclean_status:s|h' output_obiclean.fa > output_all.fa
```

Audit: determine how many reads make it through each filtering steps:    
Reads in the fasta after size selecting/count filter:   
`grep -E '^>' 04_samples/*_ali_assi_uniq_c10_55-75.fa | awk -F"\ count=" '{ print $2 }' - | awk -F";" '{ print $1 }' - | paste -sd+ - | bc`     
Reads in the fasta after only keeping head and singletons:      
`grep -E '^>' all_files_ali_assi_uniq_c10_55-75_clean_HS.fa | awk -F"; count=" '{ print $2 }' - | awk -F";" '{ print $1 }' - | paste -sd+ - | bc`    

Note: can also test other lengths by streaming into grep: 
`obigrep --lmin 55 --lmax 75 -p 'count>=5' 04_samples/yourfile.fa | grep -cE '^>' - `

### 2.3. Export data     
Use obitab to output a tab-delimited text file that will be used as an input to the R Script below.   
`./01_scripts/05_obitab_export.sh`    
(In brief: `obitab --output-seq 04_samples/*clean_HS.fa > 04_samples/*.txt`)   

### 2.4. Assign each sequence to a taxon
Use a blastn to align the H and S fasta file against nt (NCBI remote).   
`blastn -db nt -query 04_samples/your_cleaned_HS.fa -out 05_annotated/your_lib_output.txt -remote -num_descriptions 10 -num_alignments 10`    

Or if running a massive blast, use parallel:   
`SEQUENCE_FILE="04_samples/your_filtered_fasta.fa" ; OUTPUT="05_annotated/your_filtered_fasta_hits.txt" ; cat $SEQUENCE_FILE | parallel -j 12 -k --block 1k --recstart '>' --pipe 'blastn -db /home/ben/blastplus_databases/nt -query - -num_descriptions 10 -num_alignments 10 ' > $OUTPUT`    

Track output:    
`grep -cE 'Query= ' 05_annotated/your_filtered_fasta_hits.txt`

Note: for MEGAN, the blast output must be standard output format, not outfmt.  

### 2.5. Annotate sequences with species   
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

### 2.6. Connect read counts and annotations    
This will use the R script `read_counts_to_annotations.R`, run interactively.   
Necessary inputs:   
Amplicon annotation output from MEGAN, and amplicon read count from `obitab`   

In brief, this will merge these two inputs, attach locations, aggregate different amplicons with same annotation, calculate proportions, save out proportion plots and count/proportion tables.    

Within here, one can apply a low expression filter to remove any counts less than 10.   

(note: currently working on improving this script to be more universal. See a larger version on `read_counts_to_annotations_HABs.R`)    
