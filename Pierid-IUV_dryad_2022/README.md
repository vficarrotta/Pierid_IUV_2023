# Pierid-IUV_dryad_2022
# Early phylogenetic origin and repeated ecological implementation of iridescent UV signaling in pierid butterflies
# READ ME: directory structure and file descriptions


#     File Structure       

vficarrotta_dryad_IUV_2022
|-Butterfly Photo Credits
| |-butterfly-photo-credits.ods
|-IUV Scoring QC
| |-data-tree_overlap_7-7-2022.csv
| |-genus_uv_7-7-2022.csv
| |-genus_uv_overlap_7-7-2022.csv
| |-uvg_7-7-2022.csv
| |-IUV_scoring_QC.r
|-Phylogeny and Ancestral State Reconstruction
| |-Phylogeny Build
| | |-bad_seqs.py
| | |-bad_taxa.py
| | |-Pieridae_7114.table
| | |-Pieridae_7114_final.outaln
| | |-Pieridae_7114_final_GB.csv
| | |-Pieridae_7114_SHL_names.tre
| |-ASR-IUV_analysis.r
| |-full-tree_species.csv
| |-Pieridae_7114_final_names_trim_dates_2-22-2022.tre
| |-plotting_ASR-IUV_output.r
| |-uv_mat_syn-edit-V2_2-22-22.csv
|-Ridge Density
| |-ridge_dens_master_v2.csv
| |-ridge_density_analysis.r
|-Species Lists
| |-gbif-species-list_synedit.csv
| |-master-sp-list_extant-only_7-7-2022.csv
|-Spectral Analysis
| |-IUV raw text files.gzip
| |-IUV_spec-analysis_5-22-22.r
|-Tables
| |-tables_master_7-7-2022.xlsx

# Description of file contents 


# Butterfly Photo Credits

Credits for butterfly images in Figure 1 ordered by their numbers.


# IUV Scoring QC

# These files produce supplemental figure 2 which demonstrates the reduced dataset
# does not lose analytical power.

data-tree_overlap_7-7-2022: all species in both the scored dataset and the constructed tree.

genus_uv_7-7-2022: all genera/species with IUV in dataset.

genus_uv_overlap_7-7-2022: all genera/species with IUV overlapping the dataset and tree.

uvg_7-7-2022: csv containing the final dataset for plotting the quality check

IUV_scoring_QC: r code used to generate the plot using the above .csv files.


# Phylogeny and Ancestral State Reconstructions
     
#      Phylogeny Build
     
     bad_seqs.py: sequences dropped from the phylogenetic analysis

     bad_taxa.py: taxa dropped from the phylogenetic analysis

     Pieridae_7114.table: Genbank sequence metadata from PyPHLAWD analysis

     Pieridae_7114_final.outaln: Sparse supermatrix from PyPHLAWD

     Pieridae_7114_final_GB.csv: Genbank accessions for sparse supermatrix from PyPHLAWD

     Pieridae_7114_SHL_names.tre: Dated, named tree from congruify analysis

ASR-IUV_analysis.r: R code for generating all ancestral state reconstruction models and analyses

full-tree_species.csv: csv of all species in tree

Pieridae_7114_final_names_trim_dates_2-22-2022.tre: copy and renamed final tree from congruify

plotting_ASR-IUV_output.r.r: R code for plotting Figure 1 tree, Figure

uv_mat_syn-edit-V2_2-22-22.csv: full dataset of screened museum specimens for IUV


# Ridge Density


ridge_dens_master_v2.csv: dataset of measured ridge density for plotting Figure 2B

ridge_density_analysis.r: R code to generate plot in figure 2B


# Species Lists


gbif-species-list_synedit.csv: species list of Pieridae including fossils

master-sp-list_extant-only_7-7-2022.csv: species list of Pieridae without fossils


# Spectral Analysis


IUV raw text files.gzip: raw spectral data in .txt format. must be unzipped

IUV_spec-analysis.r: R code for wrangling and plotting spectral information


# Tables


tables_master_7-7-2022.xlsx: spreadsheet with two panes. first pane shows model scores of strict presence/absence and the estimated
                              transition rates of the highest scoring model. the second pane shows models scores including male and female
                              scores and the estimated transition rates of the highest scoring model.
