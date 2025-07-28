#!/bin/sh

FOLDERS=("0_rawdata" "0_rawdata/raw_qc" "1_trim" "2_mapping" "3_postmapping" "4_peakcall" "5_homer" "6_gem")

# Loop through folder names and create if they don't exist
for folder in "${FOLDERS[@]}"; do
    if [ ! -d "$folder" ]; then
        mkdir -p "$folder"
    fi
done

# Define input directory name
input_dir="0_rawdata"
rawqc_dir="0_rawdata/raw_qc"
trim_dir="1_trim"
map_dir="2_mapping"
postmapping_dir="3_postmapping"
macs2_dir="4_peakcall"
homer_dir="5_homer"
gem_dir="6_gem"

#Define adaptor and ref_genome
mouse_genome="$HOME/refgenome/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa"
human_genome="$HOME/refgenome/homo_genome/Homo_sapiens.GRCh38.dna.primary_assembly_chrM.fa"
adapter_file="$HOME/refgenome/adaptor/Common_Sequence_Adaptors.fa"

# rawdata fastqc
ls "$input_dir"/*.fq.gz | xargs -I {} -P 8 sh -c 'input_file="{}"; name="$(basename "${input_file%.fq.gz}")";
# Define fastqc output file names
fastqc_html="$rawqc_dir/${name}_fastqc.html";
fastqc_zip="$rawqc_dir/${name}_fastqc.zip";
# check if file names exist, skipping
if [ -f "$fastqc_html" ] && [ -f "$fastqc_zip" ]; then
    echo "fqc files already exist for $input_file, skipping this step";
else
    fastqc -t 20 -o "$rawqc_dir" "$input_file";
fi'


# adaptor trimming
for input_file in "$input_dir"/*_1.fq.gz; do

echo "Input file $input_file is paired-end sequencing data"
 # Define fastp output file names
        sample_name="$(basename "${input_file%_1.fq.gz}")"
        input_file_2="$input_dir/${sample_name}_2.fq.gz"
        trim_file_1="$trim_dir/${sample_name}_trim_1.fq.gz"
        trim_file_2="$trim_dir/${sample_name}_trim_2.fq.gz"
        json_file="$trim_dir/${sample_name}.json"
        html_file="$trim_dir/${sample_name}.html"

        if [ -f "$json_file" ] && [ -f "$html_file" ]; then
    echo "${sample_name} already trimmed and filtered, skipping fastp step"
else
# Run fastp with paired-end options
     fastp -i "$input_file" -I "$input_file_2" -o "$trim_file_1" -O "$trim_file_2" --thread 15 --detect_adapter_for_pe --adapter_fasta "$adapter_file" --trim_poly_g --trim_poly_x --json "$json_file" --html "$html_file"  -w 15
     fi

# Define mapping output file names
       map_file="$map_dir/${sample_name}.sam"
       map_log="$map_dir/${sample_name}_map.log"

# Define filter output file names
       bam_file="$map_dir/${sample_name}_first.bam"
       sort_file="$map_dir/${sample_name}_first_sort.bam"
       flag_file="$map_dir/${sample_name}.flagstat"
       dedup_file="$map_dir/${sample_name}_dedup.bam"
       markdup_metrics="$map_dir/${sample_name}_markdup.metrics"
       #rmdup_file="$map_dir/${sample_name}_rmdup.bam"
       rmdup_sort_file="$map_dir/${sample_name}_rmdup_sort.bam"

# Run trim_files with paired-end mapping options
       if [ -f "$markdup_metrics" ] && [ -f "$rmdup_file" ] && [ -f "$rmdup_sort_file" ]; then
    echo "${sample_name} already mapping and sort, skipping mapping and filter step"
else
       bwa mem -t 15 "$human_genome" "$trim_file_1" "$trim_file_2" > "$map_file" 2> "$map_log"

# Run trim_files with paired-end filter options
       samtools view -h -@ 25 -b "$map_file" >  "$bam_file"
       samtools sort -@ 25 "$bam_file"  >   "$sort_file"
       samtools flagstat "$sort_file"  > "$flag_file"
       picard MarkDuplicates INPUT="$sort_file"   OUTPUT="$dedup_file" M="$markdup_metrics" REMOVE_DUPLICATES=false
       samtools view -@ 25 -h -f 2 -F 1804  -b "$dedup_file"  > "$rmdup_sort_file"
       samtools index "$rmdup_sort_file"
       samtools view -@ 25 -b "$rmdup_sort_file" chr{1..22} chrX chrY > tmp.bam
       mv tmp.bam "$rmdup_sort_file"
       samtools index "$rmdup_sort_file"
       #samtools sort -@ 25  "$rmdup_file"  > "$rmdup_sort_file"
fi
done
