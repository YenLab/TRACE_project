FOLDERS=( "3_postmapping" "3_postmapping/1_bw_cpm" "3_postmapping/2_cor_file" "3_postmapping/3_plot" )

# Loop through folder names and create if they don't exist
for folder in "${FOLDERS[@]}"; do
    if [ ! -d "$folder" ]; then
        mkdir -p "$folder"
    fi
done

# Define input directory name
map_dir="2_mapping"
postmapping_dir="3_postmapping"
cpm_dir="3_postmapping/1_bw_cpm"
cor_dir="3_postmapping/2_cor_file"
plot_dir="3_postmapping/3_plot"

#reference-point
blacklist="/share/home/sango/refgenome/homo_genome/hg38-blacklist.v2.bed"
declare -A ref_points=(
  ["tss"]="/share/home/sango/refgenome/homo_genome/hg38.ENCODE/hg38.TSS.bed"
  # ["tes"]="/share/home/sango/refgenome/homo_genome/hg38.ENCODE/hg38.TTS.bed"
  # ["dELS"]="/share/home/sango/refgenome/homo_genome/hg38.ENCODE/hg38.ENCODE_cCREs_dELS.bed"
  # ["pELS"]="/share/home/sango/refgenome/homo_genome/hg38.ENCODE/hg38.ENCODE_cCREs_pELS.bed"
  # ["PLS"]="/share/home/sango/refgenome/homo_genome/hg38.ENCODE/hg38.ENCODE_cCREs_PLS.bed"
  # ["CTCF"]="/share/home/sango/refgenome/homo_genome/hg38.ENCODE/hg38.ENCODE_cCREs_CTCFonly_CTCFbound.bed"
)

# post_mapping beginning with bamCoverage
    for sort_file in "$map_dir"/*rmdup_sort.bam; do
        bigwig_file="$postmapping_dir/$(basename "${sort_file%_sort.bam}")_bs200.bw"
        bigwig_file_cpm="$cpm_dir/$(basename "${sort_file%_sort.bam}")_bs200_cpm.bw"
        # check if bw files exist, if they do, skip this step
        if [ -f "$bigwig_file" ] ; then
            echo "bw files already exist for postmapping, skipping this step"
            continue
        fi

        bamCoverage --bam "$sort_file" -o "$bigwig_file"  --blackListFileName "$blacklist"  --binSize 200 -p 15
        bamCoverage --bam "$sort_file" -o "$bigwig_file_cpm"  --blackListFileName "$blacklist" --normalizeUsing CPM --binSize 200 -p 15

      done

# post_mapping beginning with correlation
for bs_value in 200 500 1000 2000; do
     npz_file="$cor_dir/bs${bs_value}.npz"
     cpm_npz_file="$cor_dir/bs${bs_value}_cpm.npz"
     pearson_file="$plot_dir/bs${bs_value}_pearson.pdf"
     cpm_pearson_file="$plot_dir/bs${bs_value}_cpm_pearson.pdf"
     spearman_file="$plot_dir/bs${bs_value}_spearman.pdf"
     cpm_spearman_file="$plot_dir/bs${bs_value}_cpm_spearman.pdf"
      # check if npz files exist, if they do, skip this step

    if [ -f "$npz_file" ] && [ -f "$pearson_file" ] && [ -f "$spearman_file" ]; then
        echo "npz files already exist for postmapping, skipping this step"
        continue
    fi
        multiBigwigSummary bins -b $postmapping_dir/*bs200.bw --blackListFileName "$blacklist" -bs "$bs_value" -p 15 -out "$npz_file"
        multiBigwigSummary bins -b $cpm_dir/*bs200_cpm.bw --blackListFileName "$blacklist" -bs "$bs_value" -p 15 -out "$cpm_npz_file"
        plotCorrelation -in "$npz_file" -o "$pearson_file" --removeOutliers  --skipZeros -c pearson --whatToPlot heatmap --plotNumbers
        plotCorrelation -in "$npz_file" -o "$spearman_file" --removeOutliers --skipZeros -c spearman --whatToPlot heatmap --plotNumbers
        plotCorrelation -in "$cpm_npz_file" -o "$cpm_pearson_file" --removeOutliers  --skipZeros -c pearson --whatToPlot heatmap --plotNumbers
        plotCorrelation -in "$cpm_npz_file" -o "$cpm_spearman_file" --removeOutliers --skipZeros -c spearman --whatToPlot heatmap --plotNumbers
    done

wait

# reference-point for mouse tss tes CTCF enhancer and promoter
for ref_point in tss tes dELS pELS PLS CTCF; do
    #matrix_file="$cor_dir/${ref_point}_matrix.mat.gz"
    cpm_matrix_file="$cor_dir/${ref_point}_cpm_matrix.mat.gz"
    genebody_matrix="$cor_dir/genebody_matrix.mat.gz"
    #pdf_file="$plot_dir/${ref_point}.pdf"
    cpm_pdf_file="$plot_dir/${ref_point}_cpm.pdf"
    genebody_file="$plot_dir/genebody.pdf"
    cpm_tab_file="$plot_dir/${ref_point}_cpm.tab"
    genebody_tab_file="$plot_dir/genebody.tab"


    # check if matrix files exist, if they do, skip this step
    if [ -f "$matrix_file" ] && [ -f "$pdf_file" ]; then
        echo "heatmap files already exist for postmapping, skipping this step"
        continue
    fi
    computeMatrix reference-point -S $postmapping_dir/*bs200.bw -R "${ref_points[$ref_point]}" --referencePoint center --missingDataAsZero -a 2000 -b 2000 --skipZeros --blackListFileName "$blacklist" -o "$matrix_file" -p 15 -q
    #computeMatrix reference-point -S $cpm_dir/*.bw -R "${ref_points[$ref_point]}" --referencePoint center --missingDataAsZero -a 2000 -b 2000 --skipZeros --blackListFileName "$blacklist" -o "$cpm_matrix_file" -p 15 -q
    computeMatrix scale-regions -S  $cpm_dir/*bs200_cpm.bw -R /share/home/sango/refgenome/homo_genome/Homo_sapiens.GRCh38.95_chrM.gtf -a 2000 -b 2000 -p 15 --blackListFileName "$blacklist"  --regionBodyLength 5000 --missingDataAsZero --skipZeros -o "$genebody_matrix"

    #plotProfile -m "$matrix_file" -o "$pdf_file"  --yMin 0
    plotProfile -m "$cpm_matrix_file" -o "$cpm_pdf_file"  --yMin 0  --outFileNameData "$cpm_tab_file"
    plotProfile -m "$genebody_matrix" -o "$genebody_file" --outFileNameData "$genebody_tab_file"

done
