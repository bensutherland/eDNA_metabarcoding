# eDNA_metabarcoding
Dependencies:
`python 2.7, gcc, python-dev packages`    
`OBITools` http://metabarcoding.org/obitools/doc/welcome.html       

Add both the obitools binary and the /export/bin folder to the $PATH to make obitools accessible everywhere

Launch Obitools
`obitools`    

Move to 02_raw_data and decompress    
`for i in $(ls *.fastq.gz ) ; do gunzip -c $i > ${i%.gz} ; done`

Prepare the interpretation files
Note: This must be done for each sequencing lane separately       
`awk -F'\t' '$1=="1" { print "SOG", $2, $5, $6, $7, "F @ NGSlib=1" } ' OFS='\t' ./interp_file_2017-10-17.txt | sed 's/\ //g' > interp_lib1_headless.txt`
`awk -F'\t' '$1=="2" { print "SOG", $2, $5, $6, $7, "F @ NGSlib=2" } ' OFS='\t' ./interp_file_2017-10-17.txt | sed 's/\ //g' > interp_lib2_headless.txt`
`awk -F'\t' '$1=="3" { print "SOG", $2, $5, $6, $7, "F @ NGSlib=3" } ' OFS='\t' ./interp_file_2017-10-17.txt | sed 's/\ //g' > interp_lib3_headless.txt`

Then add header to each
`for i in *headless* ; do cat header.txt $i > ${i%_headless.txt}.txt ; done`


Recover full seq reads from F and R reads
``



obitaxonomy 


