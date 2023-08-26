# Reciprocal Blast Full Pipeline
Full pipeline from raw nanopore sequence data through to reciprocal blasting of target genes.

Made to de novo assemble raw nanopore reads with Flye from Dictystoelium discoideum strains, polish with 4x Racon and 1x Medaka and extract tgr genes via BLAST. Checks if the extracted genes are reciprocal best hits with the reference genes.

*Note: this still uses hardcoded filepaths and needs to be cleaned up and finalised.

## General Outline

<img width="492" alt="Screen Shot 2023-08-26 at 4 06 23 PM" src="https://github.com/tehchoobworl/ReciprocalBlastFullPipeline/assets/116311825/3a554d91-9930-4eff-ba4d-ca2a7a9d990c">


## Bash files

### fast5assembler.sh
Input: Directory containing folders with .fast5 (raw nanopore reads) inside. Each folder should be named as with a strain and contain nothing but fast5 files.

What it does: Creates the basecalled, assemblies, and polished directories. Loops through all the raw nanopore read folders and basecalls them with guppy, saving the output to the basecalled directory. This step may take an extremely long time depending on your computer specs. After basecalling, will de novo assemble all basecalled reads with Flye and save the output to the assemblies directory.

Outputs: Basecalled reads and de novo Flye assemblies in their respective folders.

```sh
# Setup file structure as so:
#	|
#	|--> fast5assembler.sh
#	|--> raw_data_directory
#					|
#					|--> EO1331
#					|		|-->FAO3052.fast5
#					|		|-->FAO3053.fast5
#					|		|-->etc...
#					|--> EO1332
#					|		|-->FAO3056.fast5
#					|		|-->FAO3057.fast5
#					|		|-->etc...

#Usage: bash fast5assembler.sh <raw_data_directory>/

#Get input directory from user
echo "Raw data directory: " $1
raw_data_directory=$1

#Make ordered directory
mkdir "basecalled"
mkdir "assemblies"
mkdir "polished"

#Loop through all files in input
for strain in $raw_data_directory*
do
	#Exctract just the title of directory, ie home/sean/data/EO167 becomes EO167
	strain_name=$(basename -- "$strain")
	strain_name="${strain_name%.*}"
	echo $strain_name
	
	#Basecall raw data using guppy
	nice -10 /home/bhargavam/ont-guppy/bin/guppy_basecaller --input_path $strain/ --save_path "basecalled/"$strain_name"_basecalled" --flowcell FLO-MIN106 --kit SQK-LSK109 -x cuda:all:100%

	#Concatenate all resulting fastq files into one
	cat basecalled/"$strain_name"_basecalled/pass/*.fastq > $strain_name.fastq
	mv $strain_name".fastq" "basecalled/"$strain_name"_basecalled/"
	
	#Assemble concatenated fastq's using flye
	nice -10 flye --nano-raw "basecalled/"$strain_name"_basecalled/"$strain_name".fastq" --out-dir "assemblies/"$strain_name"_flye_assembly"

done
```
### medakapolish.sh
Input: The same raw directory for fast5assembler.sh, as well as the de novo Flye assemblies produced by it. Also requires the AX4 reference genome for a final comparison.

What it does: Uses minimap2 to map basecalled reads to de novo assemblies and polish with Racon. Repeats 4x for racon, before mapping the basecalled reads to the new racon assembly and polishing once with Medaka. Finally, uses Quast to generate a comparison folder between the AX4 reference genome, Flye, Racon and Medaka assembly.

Output: All racon polished assemblies in ./polished/racon, all medaka polished assemblies in ./polished and all quast assemblies in ./comparisons

```bash
#Usage: bash medakapolish.sh <raw_data_directory>

echo "Raw data directory: " $1
raw_data_directory=$1

###########################################################
#Racon polishing iterations
maxiter=4

#CPU Threads
threads=24

###########################################################

for strain in $raw_data_directory*
do

	#Extract just the title of directory, ie /data/sean/raw/EO167 becomes EO167
	strain_name=$(basename -- "$strain")
	strain_name="${strain_name%.*}"

	reads="basecalled/"$strain_name"_basecalled/"$strain_name".fastq"
	contigs="assemblies/"$strain_name"_flye_assembly/assembly.fasta"

	echo "############################################"
	echo "Starting Minimap2 and Racon for $strain_name"
	echo "############################################"

	for ((i=1; i<=$maxiter; i++))
	do

		#Map reads to assembly, generate temporary overlap file
		/data/sean/minimap2/minimap2 -x map-ont -t $threads $contigs "$reads" > tmp.ovl.paf

		#Racon polish with overlap file and assembly
		racon -t $threads "$reads" tmp.ovl.paf $contigs > tmp.$i.fasta

		#Continue with new contig file
		contigs=tmp.$i.fasta

		echo "Racon round $i finished."

	done
	mv tmp.$maxiter.fasta polished/racon/$strain_name.racon_out.fasta

	echo "############################################"
	echo "$strain_name Racon finished, starting Medaka"
	echo "############################################"

	echo "Info: Run Minimap for Medaka"
	tmpsam="medakatmp.sam"
	/data/sean/minimap2/minimap2 -a -x map-ont -t $threads polished/racon/$strain_name.racon_out.fasta "$reads" > $tmpsam
	
	echo "Converting SAM to BAM"
	x=$tmpsam
	samtools view -@ 8 -bhS $x > $x.bam
	samtools sort -@ 8 $x.bam -o $x.s.bam
	samtools index $x.s.bam

	# remove intermediate unsorted bam
	rm $x.bam
	
	if [ "$strain_name" != "raw/" ];then
		nice -10 medaka_consensus -i "$reads" -d polished/racon/$strain_name.racon_out.fasta -o /data/sean/polished/$strain_name.polished -t $threads -m r941_min_high_g360 2 > ./logs/$strain_name".flye.medaka.log"
	fi

	echo '#################################################'
	echo "$strain_name polishing finished"
	echo '#################################################'

	echo '#################################################'
	echo "Starting QUAST comparison for $strain_name"
	echo '#################################################'

	quast/quast.py "assemblies/"$strain_name"_flye_assembly/assembly.fasta" polished/racon/$strain_name.racon_out.fasta polished/$strain_name.polished/consensus.fasta -r dictygff.fasta -l "Flye, Racon, Racon+Medaka" -o comparisons/$strain_name.comparison
done

#Clean Minimap and Racon up
rm tmp.ovl.paf
rm $tmpsam
rm tmp.*.fasta

echo '#################################################'
echo "All strains finished!"
echo '#################################################'
```
### getgenes.sh
Input: Directory with reference .gff (In this case, dicty_gff3/ downloaded from (http://dictybase.org/Downloads/)http://dictybase.org/Downloads/)

What it does: Creates a .fasta from each chromosome, greps all annotations that include the term 'tgr' in the title and creates a .fasta for each from the associated coordinates and saves them all

Output: All tgr genes in the gff in fasta format to alltgr.fasta, and all inidividual tgr genes in their own .fasta files in ./alltgr_inidividual

```bash
#!/bin/bash

#Get gff directory
echo "Directory: " $1
directory=$1

for gff in $directory*.gff
do
	strain_name=$(basename -- "$gff")
	strain_name="${strain_name%.*}"

    #Get chromosome fastas from gffs
    sed -e '1,/##FASTA/ d' $gff > $strain_name.fasta

    #Get all annotations with tgr in title
    grep "tgr" $gff > $strain_name.tgrGenes.gff

    #Create fasta files from tgr gene coordinates
    bedtools getfasta -fo $strain_name.tgrGenes.fasta -fi $strain_name.fasta -bed $strain_name.tgrGenes.gff

    #Adding gene name, ie tgr04 to fasta files
    grep ">" $strain_name.tgrGenes.fasta > temp_titles.txt
	awk '{print $9}' $strain_name.tgrGenes.gff > temp_names.txt

    sed -i -e "s/Name=//g" temp_names.txt

    i=1
    while read p;
    do
        name=$(sed -n "${i}p" temp_names.txt)
        search=$p
        replace=${p}":"${name}
        sed -i -e "s/$search/$replace/g" $strain_name.tgrGenes.fasta
        i=$(expr $i + 1)

    done < temp_titles.txt
done

#Tidying up
cat *.tgrGenes.fasta > alltgr.fasta
mkdir generated
mv *.fasta generated
mv *.fai generated
mv *.gff generated
mv generated/alltgr.fasta .

awk -F ';' '/^>/ {F=sprintf("%s.fasta",$2); print > F;next;} {print >> F;}' < alltgr.fasta
mkdir alltgr_individual
mv tgr*.fasta alltgr_individual
```
### reciprocalblast.sh
Input: Directory with polished assemblies, directory with reference genes (tgr genes) and reference genome

What it does: Makes a BLAST database from each polished assembly and loops through all reference genes, blasting each against every database. Hits are saved to ./all_hits and the top hit (listed first in each .tsv) is taken and BLASTed against the reference genome. The reciprocal hits from this search are saved in ./all_recip_hits.

Outputs: All hits from inital blast of each gene to databases are saved to ./all_hits, the top hit from each search is saved to ./best_hits, and the reciprocal hits are saved to ./all_recip_hits.

```bash
echo "Assemblies: " $1
assemblies=$1
echo "References: " $2
references=$2

####### Variables #######
reference_genome="/data/sean/dictygff.fasta"
#########################


#Loop through assemblies and make a database for each
for file in $assemblies*.fasta
do
	filename=$(basename -- "$file")
	filename="${filename%.*}"
	
	/data/sean/ncbi-blast-2.13.0+/bin/makeblastdb -in $file -dbtype nucl
	
	#Loop through references and blast each against database, ouput .results.fasta file containing all hits
	for ref in $references*.fasta
	do
		ref_filename=$(basename -- "$ref")
		ref_filename="${ref_filename%.*}"
		
		/data/sean/ncbi-blast-2.13.0+/bin/blastn -db $file -query $ref -out temp.tsv -outfmt '6 std sseq'
		cut -f 2,13 temp.tsv | awk '{printf(">%s\n%s\n", $1, $2); }' > "$filename.$ref_filename.results.fasta"
		
		#Take only the top result and output best_result.fasta
		head -n 2 "$filename.$ref_filename.results.fasta" > "$filename.$ref_filename.best_result.fasta"
	done
done

#Move all results to own folder
mkdir all_hits
mv *.results.fasta all_hits

#Move top results to own folder
mkdir best_hits
mv *.best_result.fasta best_hits

#Blast the best results back against the reference genome
/data/sean/ncbi-blast-2.13.0+/bin/makeblastdb -in $reference_genome -dbtype nucl
for hit in best_hits/*.fasta
do
	hit_name=$(basename -- "$hit")
	hit_name="${hit_name%.*}"
	
	/data/sean/ncbi-blast-2.13.0+/bin/blastn -db $reference_genome -query $hit -out $hit_name.all_recip.tsv -outfmt '6 std sseq'
done

#Evaluate reciprocals - do the reciprocals return the originial reference?
#Move evaluate.sh here???

#Move recips to own folders
mkdir all_reciprocal_hits
mv *.all_recip.tsv all_reciprocal_hits


#Cleanup temp file
rm temp.tsv
```
### evaluate.sh
Input: Directory containing reference genes in .fasta format and a directory containing reciprocal hits in .tsv format.

What it does: Gets the start/end positions of all reference genes and the start/end positions of all reciprocal hits from reciprocalblast.sh and compares them. If the reciprocal hit is within 50nt of the reference gene then they are considered reciprocal best hits.

Output: A .tsv file for each gene in ./eval_results with the start and end positions of the best reciprocal hit for each strain, with a 'pass' or 'fail' depending on whether they were reciprocal best hits.

```bash
#!/bin/bash
pass_strains=0
total_strains=0
for gene in tgr_genes_final/*.fasta
do
	gene_name=$(basename -- "$gene")
	gene_name="${gene_name%.*}"
    
    #Store start and end positions as variable
    gene_start=$(head -n 1 $gene | cut -d ':' -f 2 | cut -d '-' -f 1)
    gene_end=$(head -n 1 $gene | cut -d ':' -f 2 | cut -d '-' -f 2)

    for recip in all_reciprocal_hits/*.$gene_name.*.*.tsv
    do
	    recip_name=$(basename -- "$recip")
	    recip_name="${recip_name%.*}"

        recip_start=$(awk -F'\t' 'NR==1{print $9}' "$recip")
        recip_end=$(awk -F'\t' 'NR==1{print $10}' "$recip")

        if [ -z "$recip_start" ] || [ -z "$recip_end" ]
        then
            continue
        fi

        # Calculate the maximum and minimum values that recip_start should be within to be considered +-50nt of gene_start
        max_value=$(echo "$gene_start + 50" | bc)
        min_value=$(echo "$gene_start - 50" | bc)

        # Check if recip_start is within 50nt of gene_start
        if [ $(echo "$recip_start >= $min_value && $recip_start <= $max_value" | bc) -eq 1 ]
        then
          start_result=1
        else
          start_result=0
        fi

        # Calculate the maximum and minimum values that recip_end should be within to be considered +-50nt of gene_start
        max_value=$(echo "$gene_end + 50" | bc)
        min_value=$(echo "$gene_end - 50" | bc)

        # Check if recip_end is within 50nt of gene_end
        if [ $(echo "$recip_end >= $min_value && $recip_end <= $max_value" | bc) -eq 1 ]
        then
          end_result=1
        else
          end_result=0
        fi

        if [ $(echo "$start_result == 1 && $end_result == 1" | bc) -eq 1 ]
        then
          result="Pass"
        else
          result="Failed"
        fi

        echo -e "$recip_name\t$recip_start\t$recip_end\t$result" >> $gene_name.recip_eval_results.tsv
    done

    pass_count=$(awk -F '\t' '$4 == "Pass" { count++ } END { print count }' $gene_name.recip_eval_results.tsv)
    total_count=$(wc -l < $gene_name.recip_eval_results.tsv)
    if [ -z "$pass_count" ]
    then
        pass_count=0
    fi
    pass_strains=$(echo "$pass_strains + $pass_count" | bc)
    total_strains=$(echo "$total_strains + $total_count" | bc)
    echo "$gene_name was a reciprocal best hit for $pass_count out of $total_count strains"

done
percent=$(echo "scale=2 ; $pass_strains / $total_strains" | bc)
echo "$pass_strains genes were reciprocal best hits out of $total_strains genes checked ($percent)"

mkdir eval_results
mv *.recip_eval_results.tsv eval_results
```



