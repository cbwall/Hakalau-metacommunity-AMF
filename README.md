# Hakalau-metacommunity-AMF  
CB Wall, CP Egan, SIO Swift, NA Hynson (2020) Three decades post reforestation has not led to the reassembly of arbuscular mycorrhizal fungal communities associated with remnant primary forests. *Molecular Ecology*.  

**Post-disturbance reassembly of arbuscular mycorrhizal fungal communities**  
Background: Here, we leverage a large-scale reforestation project on Hawai‘i Island in the Hakalau Forest National Wildlife Refuge (hereafter, 'Hakalau') underway for over three decades to assess whether arbuscular mycorrhizal (AM) fungal communities have concurrently been restored. The reference ecosystem for this restoration project is a remnant montane native Hawaiian forest that provides critical habitat for endangered birds. We sampled soils from 12 plots within remnant and restored forest patches and characterized AM fungal communities using high throughput Illumina MiSeq DNA sequencing.  
   
## File Directory  
  
**R project** = 'Hakalau metacommunity AMF.Rproj' = load the R project for all code to run from a common directory  
**R.markdown** = 'Hakalau metacommunity.Rmd' = file with all code to produce analyses  
**R.markdown.html** = 'Hakalau-metacommunity.html' = knitr product from Rmd  
**download the html open in Chrome or other html viewer. Not directly viewable in GitHub**

- **data**  
  - subfolder: 'ecology data'  
      - data for producing ecology figures -- heat map, dbh  
      - data include % cover, abbreviation names, and metadata  
  - subfolder: 'soil chem info'  
      - data for soil chemistry parameters  
      - csvs include figures showing where samples were taken in each plot and the samples used in analysis  
- **figures**  
  - figures exported during code execution  
      - subfolder 'execute code'  
  - final figures (i.e., those formatted and used in the manuscript)  
      - subfolder 'manuscript figures'  
- **output**    
  - dataframes exported during code execution  
  - subfolder 'spiec_easi_files'  
      - the network files produced by 'SpiecEasi'. These are the R.data files that can be loaded to produce plots and analyses
- **photos**  
  - contains selected photos from Hakalau forests  
- **scripts**   
  - independent R scripts that were used to compile Rmarkdown  
- **tables**  
  - compiled tables from executed code and represented in manuscript: 'Wall et al_MolEcol_Hak_AMF_Tables.docx'  




