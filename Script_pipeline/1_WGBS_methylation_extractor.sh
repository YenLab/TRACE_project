FOLDERS=("00_rawdata" "00_rawdata/raw_qc" "01_trim" "01_trim/trim_qc" "02_mapping" "03_postmapping" "04_methylation_extractor" )

# Loop through folder names and create if they don't exist
for folder in "${FOLDERS[@]}"; do
    if [ ! -d "$folder" ]; then
        mkdir -p "$folder"
    fi
done

# Define input directory name
input_dir="00_rawdata"
rawqc_dir="00_rawdata/raw_qc"
trim_dir="01_trim"
trim_qc_dir="01_trim/trim_qc"
map_dir="02_mapping"
postmapping_dir="03_postmapping"
methy_dir="04_methylation_extractor"
export PATH=/share/home/sango/software/cgmaptools:$PATH

#Define adaptor and ref_genome
mouse_genome="$HOME/refgenome/mouse_genome"
human_genome="$HOME/refgenome/homo_genome"
adapter_file="$HOME/refgenome/adaptor/Common_Sequence_Adaptors.fa"
cgmap_bed_bs50="/share/home/sango/refgenome/homo_genome/hg38_CGmap_bed/Homo_sapiens.GRCh38.95_chrM_gtf_bs50_cgmaptools_test.bed"

for input_file in "$map_dir"/*_dedup_sort.bam; do

sample_name="$(basename "${input_file%_dedup_sort.bam}")"
methy_log="$methy_dir/${sample_name}_methy.log"
cpg_file="$methy_dir/${sample_name}.CpG_report.txt.gz"
cgmap_file="$methy_dir/${sample_name}.CGmap"
mfg_file="$methy_dir/${sample_name}.mfg"


# for wgbs data methylation extractor human data
bismark_methylation_extractor "$input_file" -p --gzip --parallel 4 --bedGraph --cytosine_report --genome_folder "$human_genome"  2> "$methy_log"
bismark2report .

cgmaptools convert bismark2cgmap -i "$cpg_file" -o "$cgmap_file"

#gzip "$cgmap_file"

cgmaptools mfg -i "$cgmap_file" -r "$cgmap_bed_bs50" -x CG > "$mfg_file"

(head -1 one.mfg | awk -F "\t" -v OFS="\t" '{$1="Sample"; print $0;}';
for mfg_file in *.mfg; do
    cat ${mfg_file} | awk -F "\t" -v OFS="\t" -v SampleName=$(echo ${mfg_file} | sed s/.mfg//g) '/total_ave_mC/{$1=SampleName;print $0}'
done
) > mfg_merge.xls
