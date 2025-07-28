FOLDERS=("0_rawdata" "0_rawdata/raw_qc" "1_trim" "1_trim/trim_qc" "2_mapping" "2_mapping/1_log")

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
map_log="2_mapping/1_log"

#Define adaptor and ref_genome
mouse_genome="/share/home/sango/refgenome/mouse_genome/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa"
human_genome="/share/home/sango/refgenome/homo_genome/Homo_sapiens.GRCh38.dna.primary_assembly_chrM.fa"
adapter_file="/share/home/sango/refgenome/adaptor/Common_Sequence_Adaptors.fa"
hg38_index="/share/home/sango/refgenome/homo_genome/star_index/hg38"
gtf_file="/share/home/sango/refgenome/homo_genome/star_index/gencode.v38.annotation.gtf"
hg38_bed="/share/home/sango/refgenome/homo_genome/star_index/hg38_GENCODE.v38.bed"

for i in $input_dir/*_R1.fq.gz; do
    input_file2="${i%_R1.fq.gz}_R2.fq.gz"
    # Define fastp output file names
    prefix="$(basename "${i%_R1.fq.gz}")"
    trim_file1="$trim_dir/${prefix}_trim_R1.fq.gz"
    trim_file2="$trim_dir/${prefix}_trim_R2.fq.gz"
    json_file="$trim_qc_dir/$prefix.json"
    html_file="$trim_qc_dir/$prefix.html"

    if [ -f "$json_file" ] && [ -f "$html_file" ]; then
        echo "${sample_name} already trimmed and filtered, skipping fastp step"
    else
        # Run fastp with paired-end options
        fastp -i "$i" -I "$input_file2" -o "$trim_file1" -O "$trim_file2" --thread 15 --detect_adapter_for_pe --adapter_fasta "$adapter_file" --trim_poly_g --trim_poly_x --json "$json_file" --html "$html_file" -w 15
    fi

    # Run mapping using STAR
    STAR --runThreadN 15 --genomeDir $hg38_index --sjdbGTFfile $gtf_file --limitBAMsortRAM 2609076260 --readFilesCommand zcat --readFilesIn $trim_file1 $trim_file2 --outFileNamePrefix $map_dir/$prefix --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20
done

# Post-mapping QC
for i in $map_dir/*Aligned.sortedByCoord.out.bam; do
    name=$(basename ${i%Aligned.sortedByCoord.out.bam})
    filter_bam=$map_dir/"${name}.filtered.bam"
    gene_tab=$map_dir/${name}ReadsPerGene.out.tab
    raw_count=$map_dir/${name}_Rawcount.txt

 if [ -f "$filter_bam" ] && [ -f "$raw_count" ]; then
        echo "${sample_name} already mapping and genecount, skipping all steps"
    else

    samtools view -@ 15 -F 780 -q 30 -b $i \
    | samtools sort -@ 15 - > $filter_bam

    samtools index $filter_bam

    tail -n +5 $gene_tab | cut -f 1,2 | sort -k1,1 > $raw_count

    fi
done

    mv $map_dir/*Log* $map_dir/*.summary $map_log

    #featureCounts -T 16 -f -t gene -B -p -C -g gene_id -a $gtf_file -o all_sample_geneID.count $map_dir/*.filtered.bam 1> featureCounts.std 2>&1
    #featureCounts -T 16 -f -t gene -B -p -C -g gene_name -a $gtf_file -o all_sample_geneNAME.count $map_dir/*.filtered.bam 1> featureCounts_genecount.std 2>&1
