Meta-analysis_NFD
Codes of the manuscript "Niche differences, not fitness differences, explain predicted coexistence across ecological groups" by Lisa Buche, Jürg Spaak, Javier Javillo and Frederik De Laender. 

(1) The data collection and how it was extracted from each paper can be found under "ExtractionData.R". We gathered 29 independent data-sets corresponding to 1018 communities (Appendix, section A.1, Tab S4 in the manuscript). For each article, we extracted all species-specific growth parameters available (competition coefficient, sensitivity, intrinsic growth rate, invasion growth rate, etc.) and the outcome of species interaction, i.e. coexistence, competitive exclusion or priority effects. The ExtractionData.R will create a data.frame "NFD_meta", containing the different columns of 2 differents classes, first the informations concerning the community: 
•	"First.Author": the first author of the paper in which we found the community;
•	"Year": the year od publication of that paper; 
•	"lab.OR.field": the setup of the experiment, either "lab" or "field"; 
•	"Model.use": the population model they used to convert the raw data into competition coefficient;
•	"Initial.Definition": the definition of NFD used in the paper;
•	"Kingdom": the functional group of the community, here either "Bacteria", "Phytoplankton","Annual plant" or "Perennial plant";
•	"treatment" : the treatment apply to the community, by the authors of the paper;
•	"species.i": the focal species of the community;
•	"species.j" : the competitor affecting species.i. 

Second the quantitative parameters computed by the authors: 
•	"LV_a_ii" - "LV_a_ij" - "LV_a_jj" - "LV_a_ji" :  the competitive coefficient obtain by the fitting of the lotka Voltera model; 
•	"AP_a_ii" - "AP_a_ij" - "AP_a_jj" - "AP_a_ji" :  the competitive coefficient obtain by the fitting of the annual plant model;  
•	"AP_lambda_i" - "AP_lambda_j" - "AP_g_i"  - "AP_s_i" - "AP_g_j" -  "AP_s_j" :  the vital parameters based on the annual plant model
•	"AP_a_ii_exp" - "AP_a_ij_exp" - "AP_a_jj_exp" - "AP_a_ji_exp" :  the competitive coefficient obtain by the fitting of the second annual plant model;  
•	"AP_lambda_i_exp" - "AP_lambda_j_exp" - "AP_g_i_exp"  - "AP_s_i_exp" - "AP_g_j_exp" -  "AP_s_j_exp" :  the vital parameters based on the second annual plant model »
•	"g_ji" -  "g_ij" -  "g_ii" -  "g_jj" : the growth rates of the species in either competition or monoculture;
•	"Si" - "Sj" :  Sensitivity of each species ( See Carroll & al, 2011)
•	"ND" - "FD" : the niche and fitness difference computed by the authors;
•	"ND_Zhao" -  "FD_Zhao" - "Zhao_Coex" : the niche and fitness difference computed by Zhao & al, 2016. 

To not share data that do not belong to us, you can find a subset of this dataframe called "Data/NFD_meta_subset.csv", without the quantitative parameters originally computed by the authors. In this dataset you can find the NFD computed with the file “ComputeDefinition_and_approximation.R”. Additionally a summary of all the data collected can be find under "Data/NFD_meta.recap.csv". 

(2) The computation of Spaak & al, 2020 definition and its approximation can be find under “ComputeDefinition_and_approximation.R”. 
It calls 4 others scripts:
•	AllDefinitions_NFD.R
•	AnnualPlant2_NFD.R
•	AnnualPlant_NFD.R
•	LotkaVolterra_NFD.R

The 3 last ones, are based on the python script : "numerical_NFD.py"

(3) The clustering analysis can be find under the "Metaanalysis_clusters_final.R" and "StatisticalTest_meta_clusters_final.R".
![image](https://user-images.githubusercontent.com/60778585/157350475-d44536ae-4fe5-4eff-b8a6-e01ec2ceeebe.png)

Files and what they do:
clustering_number_of_clusters.py
	Computes the clustering and estimate niche and fitness differences for where this is not known

plot_data_comp_outcome.py
	Plot figure 1

plot_clusters.py
	Plots figure 2

plot_correlation_study_size_coexistence.py
	Create figure S2

plot_number_of_clusters.py
	Plots figure S3
