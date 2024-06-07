# PatDetect
PatDetect is a bioinformatics pipeline that combines various taxonomic classification tools to enable integrated metagenomic analysis of Illumina short paired-end (PE) and Oxford Nanopore Technologies (ONT) long reads.The operating principle is the same as the MetaAll classification methods, but with a larger number of tools and more complex databases, so it is only suitable for operation on HPC.

## Installation & Dependencies
To obtain the scripts, download repository using `git clone` or `wget` and additionally install:
- Snakemake workflow management system (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- Singularity (https://singularity-tutorial.github.io/01-installation/)
- MEGAN (https://software-ab.cs.uni-tuebingen.de/download/megan6/welcome.html)
  
**NOTE:** Build Singularity images from definition files (.def and .sif file must have same name) in `singularity_images` folder.

### Obtain the required databases
Download required databases:
- KrakenUniq MicrobialDB collection (https://benlangmead.github.io/aws-indexes/k2)
- Kaiju nr_euk database (https://bioinformatics-centre.github.io/kaiju/downloads.html)
- MetaPhlan 4 database (http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/)
- virus/chromosome-specific HMMs  (https://figshare.com/ndownloader/files/17904323?private_link=f897d463b31a35ad7bf0)
- geNomad database (https://zenodo.org/records/8339387)
- NCBI nt (https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz)
- NCBI nr (https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz)
- NCBI prot.accession2taxid (https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz)
- NCBI names/nodes (https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip)
- MEGAN (https://software-ab.cs.uni-tuebingen.de/download/megan6/megan-map-Feb2022.db.zip)

To obtain and set up the nt database for the Krakenuniq and Centrifuge tools, follow the instructions for Krakenuniq (https://github.com/fbreitwieser/krakenuniq/blob/master/README.md) and Centrifuge (https://github.com/khyox/recentrifuge/wiki/Centrifuge-nt).

  
**NOTE:** Make sure you have enough disk space. For optimum performance, 1TB RAM and 32 CPUs are required. 

## Example of use
Before every run double check workflow parameters and path to samples and databases.
Once set, simply run `bash run_workflow.sh`
### Short PE reads/contigs classification
In terminal, navigate to the `short_reads_classification/` folder, which contains  `config.yml`,`run_workflow.sh` and `Snakefile`.
Workflow performs quality check, trimming, host removal, assembly, read/contig classification and visualization preparation of results.
Before run, set the parameters in `config.yml` file and `run_workflow.sh` script. Check PE reads name (must end with "_R1.fastq.gz" and "_R2.fastq.gz").
Host reference genome must be indexed (use `bowtie2-build` command). 

**Suggestion:** Before run use "-n" flag in shell scripts, to perform dry-run.

### Long reads/contigs classification
In terminal, navigate to the `long_reads_classification/` folder, which contains  `config.yml`,`run_workflow.sh` and `Snakefile`.
Workflow performs quality check, trimming, host removal, assembly, polishing, read/contig classification and visualization preparation of results.
Before run, set the parameters in `config.yml` file and `run_workflow.sh` script. 

**IMPORTANT:** Check input path (the defined path must end above the folder containing reads).
For example:
if raw reads are located in `../path_to_sequence_run/fastq_pass/barcode01`, the defined path in `config.yml` must be:
```
../path_to_sequence_run/fastq_pass
```
Rename folder if you wish (e.g. rename "barcode01" to "sample01")

**Suggestion:** Before run use "-n" flag in shell scripts, to perform dry-run.
