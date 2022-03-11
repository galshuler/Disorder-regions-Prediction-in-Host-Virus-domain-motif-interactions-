## Disorder regions prediction - Host-Virus domain-motif interactions
Run the script "Code_Disorder_regions_prediction_Host-Virus_domain-motif_interactions" and make sure the that all the input files and the Iupred script and data are extract to the same folder.

This project is run in python 3.9 Using the following packages:
 - biopython
 - matplotlib
 - numpy
 - pandas
 - requests
 - scipy
 - seaborn
 - re
 - os
 
The script only work on Python 3 and must have an active internet connection (for download files from the website).

### input data
- elm_classes.tsv
- elm_instances.tsv

### Recommended run order
Run the script in the order of the functions:
- a.) make_dictionary_match_uniprot_regex - A function for analysis of information between two different dataset : elm instances and elm classes by tsv format.
-	b.) download_uniport_seq - A function for download the species protein sequences of the unique uniprot id from uniprot server.
-	c.) apply_iupred_software - A function for run Iupred software for the purpose of predicting regions of protein that are in disorder area.
-	d.) disorder_value_motif_start_end - A function for search the motif in the protein sequnce by the Regex and integration to the result of the Iupred prediction for each motif.
-	e.) make_figures - A Function for creating figures based on the results of the analysis.
-	makdir - A function for creating a directory even if it already exists.

In order that the function "apply_iupred_software" will work you must that python3 app will appear in "App execution aliases" for that you can run python3 from the command promot.
If it is not installed it is possible to download for free from Microsoft Store.  

### Iupred script and data
- The scripts iupred2a iupred2a_lib and data folder use to run iupred software.

### output data
- Fasta_seq_protein_uniprot_virus.fasta.
- Fasta_seq_protein_uniprot_human.fasta.
- folder sequence_uniprot_virus.
- folder sequence_uniprot_human.
- folder result_iupred_virus.
- folder result_iupred_human.
- result_figures - scatter plot , box plot , histogram.

### Summary
- Project overview and conclusion.pdf
- It should take about 30 minutes to finish running the script.

