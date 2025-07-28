
FOLDERS=("0_rawdata" "0_rawdata/raw_qc" "1_trim" "1_trim/trim_qc" "2_mapping" "3_postmapping" )

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
trim_qc_dir="1_trim/trim_qc"
map_dir="2_mapping"
postmapping_dir="3_postmapping"


#Define adaptor and ref_genome
mouse_genome="$HOME/refgenome/mouse_genome"
human_genome="$HOME/refgenome/homo_genome"
adapter_file="$HOME/refgenome/adaptor/Common_Sequence_Adaptors.fa"

# rawdata fastqc
for input_file in "$input_dir"/*.fq.gz; do
 # Define fastqc output file names
 fastqc_html="$rawqc_dir/$(basename "${input_file%.fq.gz}")_fastqc.html"
 fastqc_zip="$rawqc_dir/$(basename "${input_file%.fq.gz}")_fastqc.zip"
 # check if file names exiss skipping
    if [ -f "$fastqc_html" ] && [ -f "$fastqc_zip" ]; then
    echo "fqc files already exist for $input_file, skipping this step"
    continue
    fi
fastqc -t 20 -o "$rawqc_dir" "$input_file"
 done

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
       map_file="$map_dir/${sample_name}_1_bismark_bt2_pe.bam"
       map_log="$map_dir/${sample_name}_bismark_bam.log"

# Define filter output file names
       dedup_file="$map_dir/${sample_name}.deduplicated.bam"
       dedup_log="$map_dir/${sample_name}_dedup.log"
       sort_file="$map_dir/${sample_name}_dedup_sort.bam"
       fragment_file="$map_dir/${sample_name}_fragment_size.png"

# Run trim_files with paired-end mapping options
       if [ -f "$dedup_file" ] && [ -f "$sort_file" ] && [ -f "$fragment_file" ]; then
    echo "${sample_name} already mapping and sort, skipping mapping and filter step"
else
      bismark  "$human_genome"  -1 "$trim_file_1"  -2  "$trim_file_2" --multicore 2 -p 8 -X 1000 -o "$map_dir"  2> "$map_log"

#remove the duplicate reads after mapping

deduplicate_bismark -p -bam "$map_file" -o "$dedup_file"  2> "$dedup_log"
samtools sort -n "$dedup_file" -@ 25 -o  "$sort_file"
samtools index "$sort_file" -@ 20
bamPEFragmentSize --bamfiles "$sort_file" -hist "$fragment_file" --maxFragmentLength 1000 --samplesLabel "$sample_name"

fi
done
