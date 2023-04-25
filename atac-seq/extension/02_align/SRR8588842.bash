# SRR20653397   NB4 Cut&Tag PML-RARA antibody replicate 1
# SRR20653396   NB4 Cut&Tag PML-RARA antibody replicate 2
# SRR8588841    anti-PML/RARa (customized by ABclonal Biotechnology)
# SRR8588842    anti-HDAC1 (Millipore, 06-720, Lot 2571576)
# SRR8588843    anti-P300 (Santa Cruz Biotechnology, sc-585x, Lot J2214)
# SRR8588844    anti-H3K27ac (Diagenode, C15410174, Lot A7071-001P)
# SRR8588845    Control one
# SRR8588846    Control two

REQ="SRR8588842"

# # 2. Align to hg 19
BOWTIE="/local/data/mphilcompbio/2022/mw894/gi2_a3_data/software/bowtie2-2.5.1-linux-x86_64/bowtie2"
INDEX="/local/data/mphilcompbio/2022/mw894/gi2_a3_data/software/index3/GRCh38_noalt_as/GRCh38_noalt_as" 

IN="/local/data/mphilcompbio/2022/mw894/gi2_a3_data/extension/01_raw/${REQ}.fastq"

OUTSAM="/local/data/mphilcompbio/2022/mw894/gi2_a3_data/extension/03_align/${REQ}.sam"

$BOWTIE -p 16 -q -x $INDEX -U $IN -S $OUTSAM # removed the k 16

# 3. process the alignement
OUTBAM="/local/data/mphilcompbio/2022/mw894/gi2_a3_data/extension/03_align/${REQ}.bam"
OUTBAMSORT="/local/data/mphilcompbio/2022/mw894/gi2_a3_data/extension/03_align/${REQ}_sorted.bam"

samtools view -h -S -b -@ 16 -o $OUTBAM $OUTSAM
samtools sort $OUTBAM -o $OUTBAMSORT
samtools index -@ 16 $OUTBAMSORT