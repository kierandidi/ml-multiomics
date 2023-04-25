# SRR20653397   NB4 Cut&Tag PML-RARA antibody replicate 1
# SRR20653396   NB4 Cut&Tag PML-RARA antibody replicate 2
# SRR8588841    anti-PML/RARa (customized by ABclonal Biotechnology)
# SRR8588842    anti-HDAC1 (Millipore, 06-720, Lot 2571576)
# SRR8588843    anti-P300 (Santa Cruz Biotechnology, sc-585x, Lot J2214)
# SRR8588844    anti-H3K27ac (Diagenode, C15410174, Lot A7071-001P)
# SRR8588845    Control one
# SRR8588846    Control two

# SETTINGS
# q: FDR <= 0.01
# g: hs, for homo sapiens genome size

for REQ in "SRR8588842" "SRR8588843"    #"SRR20653397" "SRR20653396"
do
    IN="/local/data/mphilcompbio/2022/mw894/gi2_a3_data/extension/03_align/${REQ}_sorted.bam"
    OUT="/local/data/mphilcompbio/2022/mw894/gi2_a3_data/extension/04_peaks/"
    macs2 callpeak -p 1e-10 -t $IN -n $REQ  --outdir $OUT
done
