# phyloGP
phylogenetic pipeline - final assignment for Comparative Genomics classes

## Usage
```
positional arguments:
  i                     File with species names (species_name in each line)
  fasturec              Path to fasturec program

optional arguments:
  -h, --help            show this help message and exit
  -n N_SPECIES, --n_species N_SPECIES
                        Number of species to draw
  --min_proteins MIN_PROTEINS
                        Minimum number of proteins in a proteome
  --max_proteins MAX_PROTEINS
                        Maximum number of proteins in a proteome
  --my_species MY_SPECIES
                        Name for file to save selected species
  --ids_file IDS_FILE   Name for file to save the Uniprot IDs table
  --proteins_file PROTEINS_FILE
                        Name for file to save all protein sequences
  --proteomes_dir PROTEOMES_DIR
                        Name for directory to save the proteomes
  --clustered_dir CLUSTERED_DIR
                        Name for directory to save the clustered sequences
  --trees_file_c TREES_FILE_C
                        Name for file to save protein family trees (bijective)
  --trees_file_s TREES_FILE_S
                        Name for file to save protein family trees (with paralogs)
  -c C_VALUE, --c_value C_VALUE
                        c parameter for MMseqs2
  --cov_mode COV_MODE   cov_mode parameter for MMseqs2
  --cluster_mode CLUSTER_MODE
                        cluster_mode parameter for MMseqs2
  -e E_VALUE, --e_value E_VALUE
                        e parameter for MMseqs2
  --mm_file MM_FILE     Name for file to save MMseqs2 results
  --max_seqs MAX_SEQS   Maximum number of sequences in a cluster
  --min_freq MIN_FREQ   Minimal frequency of a split to add to consensus tree
  --consensus_file CONSENSUS_FILE
                        Name for file to save consensus tree
  --supertree_file1 SUPERTREE_FILE1
                        Name for file to save basic supertree
  --supertree_file2 SUPERTREE_FILE2
                        Name for file to save supertree (from trees with paralogs)
  --boot_file BOOT_FILE
                        Name for file to save trees with bootstrap supports
  --boot_dir BOOT_DIR   Name for directory to compute bootstrap trees
  --bn BN               Number of bootstrap replicates
  --good_bootstrap_file GOOD_BOOTSTRAP_FILE
                        Name for file to keep trees with good bootstrap
                        supports
  --b_consensus_file B_CONSENSUS_FILE
                        Name for file to save bootstrap consensus tree
```
