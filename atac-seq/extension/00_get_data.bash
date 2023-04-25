# SRR20653397   NB4 Cut&Tag PML-RARA antibody replicate 1
# SRR20653396   NB4 Cut&Tag PML-RARA antibody replicate 2
# SRR8588841    anti-PML/RARa (customized by ABclonal Biotechnology)
# SRR8588842    anti-HDAC1 (Millipore, 06-720, Lot 2571576)
# SRR8588843    anti-P300 (Santa Cruz Biotechnology, sc-585x, Lot J2214)
# SRR8588844    anti-H3K27ac (Diagenode, C15410174, Lot A7071-001P)
# SRR8588845    Control one
# SRR8588846    Control two

OUTDIR="/local/data/mphilcompbio/2022/mw894/gi2_a3_data/07_extension/01_raw"
cd $outdir

for REQ in "SRR8588841" "SRR8588842" "SRR8588843" "SRR8588844" "SRR8588845" "SRR8588846" "SRR20653397" "SRR20653396"
do
    fasterq-dump --outdir $OUTDIR --threads 16 --progress $REQ -t $OUTDIR
done