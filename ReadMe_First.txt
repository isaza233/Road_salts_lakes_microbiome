The codes are written to work directly on your computer. To do so, copy the integrality of the Road_salts_lakes_microbiome folder to your computer in any folder you like. Then, open the Road_salts_lakes_microbiom.Rproj file with RStudio to access any other .R file. The filepaths used in the code are relative to the folder were the .Rproj file is, so do not modify the hierarchy of the folders nor change their names. 


***** Data *****

Chemistry
	- chem.csv
	- chem.csv_READ.ME.txt 

Chloroplasts reads OTU table
	- Chloroplasts.csv
	- Chloroplasts.csv_READ.ME.txt

18S rRNA reads OTU table 
	- Euk_samples_PR2.csv
	- Euk_samples_PR2.csv_READ.ME.txt

16S rRNA reads OTU table 
	- Prok_abondance_taxo.csv
	- Prok_abondance_taxo.csv_READ.ME.txt

Microscopy 
	- Sommaire2.csv
	- Sommaire2.csv_READ.ME.txt

Pictures of organisms identified by microscopy 
	- cells_ID_suivi_annuel2.pptx
	- You may use the pictures if you link either my name or the paper, or both if applicable. 


***** R code *****

Figure 1
	- map2.R

Figure 2
	- NMDS_2017.R

Figure 3
	- mix_analysis_color.R

Figure 4 
	- PERMANOVA_TUKEY_2017_v3.R

Figure 5
	- rda2017_seasons3.R

Figures 6 and 7 
	- heatmap_DESeq.R
	- dependant on 4 files: rect_heatmap_16S_cond_L2FC_2017.R, rect_heatmap_18S_cond_L2FC_2017.R, rect_heatmap_16S_seas_L2FC_2017.R, rect_heatmap_18S_seas_L2FC_2017.R