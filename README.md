# ml-multiomics
Project for MPhil Computational Biology at Cambridge University 2023, Genomics II

Reproduction of Villiers *et. al.*: Multi-omics and machine learning reveal context-specific gene
regulatory activities of PML::RARA in acute promyelocytic leukemia, *Nature Communications **14***, 2023

### Resources

- [Paper](https://www.nature.com/articles/s41467-023-36262-0)
- [Code](https://github.com/borimifsud/REBEL)
- [Data Spreadsheet](https://zenodo.org/record/7467566#.ZCLbL-zML0o)

### Tasks

<details>
  <summary>1. Gene expression (RNA-Seq): Martina</summary>
  
1. COMPARE: Source Data Tab 1 U937-PR9 DEGs 
2. 1b: RNA-seq volcano plot 
3. 1c: GO terms 
4. 1d: Pathway enrichment 
5. S1b: Bar plot of numbers of fusion transcripts in patients 
6. S1c: Correlation heatmaps of RNA-seq samples 
7. S1d-S1e: Scatter plot of induced/uninduced replicates 
8. S1f: MDS plot of replicates
</details>

<details>
  <summary>2. PML::RARA binding (Cut&Run): Yingtong</summary>
  
1. [] COMPARE: Source Data Tab 2 U937 PR Binding Sites 
2. 2b: Pie charts of genomic distribution of binding peaks 
3. 2c: Venn diagram of bound/DEG 
4. 2d: Ranked plot: expression of bound genes + bar plot 
5. S2a: Venn diagram of peak overlaps 
6. S2b: Jaccard heatmap of replicates 
7. S2c: Table of peak numbers 
8. S2g-S2h: Venn diagrams of overlapping peaks 
9. S2i: Notched box plot comparing induced/uninduced peak strength with p-value 
10. S2j-l: Venn diagrams of overlapping peaks 
11. S2m-S2n: notched boxplots of strength/width of induced/uninduced peaks 
12. S2o: Boxplots of bound genes in DEG 
13. S2p: Scatter plot of expression in key genes in two sample types 
14. S2s-v: notched boxplots 
</details>

<details>
  <summary>3. Chromatin occupancy (ATAC-seq): Max</summary>
  
1. COMPARE: Source Data 
  - Tab 3 U937 Induce Merg CHiCAGO 
  - Tab 4 U937 Uninduced Merg CHiCAGO 
  - Tab 5 U937 Gain GOTHiC Interact 
  - Tab 6 U937 Lost GOTHiC Interact 
2. 3a: Venn diagram of uninduced/induced interactions 
3. 3b+3c: Bar plots 
4. 3d: Venn diagram of genes associated with consistently lost/gained interactions 
5. 3e+3f: Venn diagram of down/up regulated genes gaining/losing interactions 
6. 3g: Bar plot 
7. 3h: Bar plot of H3k914av levels 
8. S3a-d: Venn diagrams of interactions in reps and conditions 
9. S3e: boxplot of interactions per gene 
10. S3f-i: Venn diagrams of lost/gained interactions overlapped with binding 
  
</details>

<details>
  <summary>4. long-range interactions (Promoter capture Hi-C): Tom</summary>
  
1. COMPARE: Source Data Tab 8 U937 ATAC Merged Peaks 
2. 4a: Volcano plot of differential ATAC 
3. 4c: Venn diagram of binding sites/ATAC overlap 
4. S4a: Venn diagrams of replicates 
5. S4b-S4c: pie charts showing genomic distribution 
6. S4d: Histogram showing positional frequency of binding relative to ATAC 
7. S4e: bar plot showing peak scores overlapping ATAC 
8. S4f-S4g: Venn diagrams overlapping ATA with down/up regulated genes 
9. S4h: boxplots with reads at ATAC peaks 

  
</details>

<details>
  <summary>5. Pattern Recognition (REBEL, XGBoost): Kieran</summary>
  
1. COMPARE: Source Data 
  - Tab 9: Lost UP U937 (Top SHAPs) 
  - Tab 10: Lost DOWN U937 (Top SHAPs) 
  - Tab 11: Lost NC (Top SHAPs) 
  - Tab 12: Gained UP U937 (Top SHAPs) 
  - Tab 13: Gained DOWN U937 (Top SHAPs) 
  - Tab 14: Gained NC U937 (Top SHAPs) 
  - Tab 15: Lost UP U937 (GenesInClus) 
  - Tab 16: Lost DOWN U937 (GenesInClus) 
  - Tab 17: Lost NC U937 (GenesInClus) 
  - Tab 18: Gained UP U937 (GenesInClus) 
  - Tab 19: Gained DOWN U937 (GenesInClus) 
  - Tab 20: Gained NC U937 (GenesInClus) 
2. 5b: tSNE of categories 
3. 5c: Bar plot of AUC for categories 
4. 5d+5e: tSNEs coloured by SHAPELY scores 
5. 5g + 5h: Bar plots and complex Venn diagrams of genes per cluster 
6. S5a: Bar plots showing enrichment of specific motifs in different conditions 
7. S5b: Bar plots showing AUC for each condition 
8. S5c: Pie chart with SHAPELY drivers 
9. S5d: Boxplot with SHAR score 
10. S5e: Motifs for top 5 predictors 
  
</details>

<details>
  <summary>Optional: Figure 6/S6: Expression and long-range interactions in APL patients</summary>
  
1. COMPARE: Source Data 
  - Tab 21: APL Patient 1 CHiCAGO 
  - Tab 22: APL Patient 2 CHiCAGO 
2. 6a: Correlation heatmap of expression in patients 
3. 6b: Bat plot of interaction overlaps 
4. 6c-6d: Venn diagram and scatter plot of patient expression 
5. S6a-S6b: Scatter plots of reads for various sample combinations 
6. S6c-e: Van diagrams overlapping long range interactions 
7. S6g: Notched boxplot Jaccard distance between patients 
8. S6h-S6i: Scatter plot of GEX between U and patients 
  
</details>

<details>
  <summary>Optional: Figure 7: Machine learning model applied to patient</summary>
  
1. COMPARE: Source Data 
  - Tab 23:  Gained UP in CD34 (Top SHAPs) 
  - Tab 24:  Gained DOWN in CD34 (Top SHAPs) 
  - Tab 25:  Gained NC in CD34 (Top SHAPs) 
  - Tab 26:  Gained UP in CD34 (GenesInClus) 
  - Tab 27:  Gained DOWN in CD34 (GenesInClus) 
  - Tab 28:  Gained NC in CD34 (GenesInClus) 
  - Tab 29:  Gained UP in APL (Top SHAPs) 
  - Tab 30:  Gained DOWN in APL (Top SHAPs) 
  - Tab 31:  Gained NC in APL (Top SHAPs) 
  - Tab 32:  Gained UP in APL (GenesInClus) 
  - Tab 33:  Gained DOWN in APL (GenesInClus) 
  - Tab 34:  Gained NC in APL (GenesInClus) 
2. 7a-d: tSNEs of clusters coloured various ways (with AUC bar plot) 
3. 7e-7f: interaction networks and tSNEs 
4. 7g-7h:Venn diagrams with GO enrichment 

  
</details>










