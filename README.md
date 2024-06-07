# PatDetect
PatDetect is a bioinformatics pipeline that combines various taxonomic classification tools to enable integrated metagenomic analysis of Illumina short paired-end (PE) and Oxford Nanopore Technologies (ONT) long reads.The operating principle is the same as the [MetaAll](https://github.com/NGS-bioinf/MetaAll) classification methods, but with a larger number of tools and more complex databases, so it is only suitable for operation on HPC.

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
- Host reference genome (e.g. hg38)

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

## Output
```
Short PE reads workflow output structure:
           # raw fastqc
           "{results}/preprocess/QC/reports_raw/{sample}_R1.html"
           "{results}/preprocess/QC/reports_raw/{sample}_R1_fastqc.zip"
           "{results}/preprocess/QC/reports_raw/{sample}_R2.html"
           "{results}/preprocess/QC/reports_raw/{sample}_R2_fastqc.zip"
           # raw multiqc
           "{results}/preprocess/QC/combined_raw/multiqc.html"        
           # fastp
           "{results}/preprocess/trimmed/{sample}_trim_R1.fastq.gz"
           "{results}/preprocess/trimmed/{sample}_trim_R2.fastq.gz"
           "{results}/preprocess/trimmed/{sample}_trim_S.fastq.gz"
           "{results}/preprocess/trimmed/tmp/{sample}_trim_S1.fastq.gz"
           "{results}/preprocess/trimmed/tmp/{sample}_trim_S2.fastq.gz"
           "{results}/preprocess/trimmed/reports/{sample}.html"
           "{results}/preprocess/trimmed/reports/{sample}.json"
           # trim fastqc
           "{results}/preprocess/QC/reports_trim/{sample}_R1.html"
           "{results}/preprocess/QC/reports_trim/{sample}_R1_fastqc.zip"
           "{results}/preprocess/QC/reports_trim/{sample}_R2.html"
           "{results}/preprocess/QC/reports_trim/{sample}_R2_fastqc.zip"
           "{results}/preprocess/QC/reports_trim/{sample}_S.html"
           "{results}/preprocess/QC/reports_trim/{sample}_S_fastqc.zip"
           # trim multiqc
           "{results}/preprocess/QC/combined_trim/multiqc.html"
           # bowtie2
           "{results}/preprocess/host_depl/{sample}_clean_R1.fastq.gz"
           "{results}/preprocess/host_depl/{sample}_clean_R2.fastq.gz"
           "{results}/preprocess/host_depl/{sample}_clean_S.fastq.gz"
           "{results}/preprocess/host_depl/{sample}.bam"
           "{results}/preprocess/host_depl/tmp/{sample}_clean_R%.fastq.gz"
           # bbmap
           "{results}/preprocess/host_depl/{sample}_clean_interleaved.fastq.gz"
           # clean fastqc
           "{results}/preprocess/QC/reports_clean/{sample}_R1.html"
           "{results}/preprocess/QC/reports_clean/{sample}_R1_fastqc.zip"
           "{results}/preprocess/QC/reports_clean/{sample}_R2.html"
           "{results}/preprocess/QC/reports_clean/{sample}_R2_fastqc.zip"
           "{results}/preprocess/QC/reports_clean/{sample}_S.html"
           "{results}/preprocess/QC/reports_clean/{sample}_S_fastqc.zip"
           # clean multiqc
           "{results}/preprocess/QC/combined_clean/multiqc.html"
           # krakenuniq
           "{results}/read_classification_results/krakenuniq_results/krakenuniq_taxonomic/{sample}_krakenuniq.krk"
           "{results}/read_classification_results/krakenuniq_results/pavian_reports/{sample}_krakenuniq.report"
           "{results}/read_classification_results/krakenuniq_results/unclassified_reads/{sample}_unclassified.fastq.gz"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/krakenuniq_taxonomic/{sample}_krakenuniq_unclass.krk"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/pavian_reports/{sample}_krakenuniq_unclass.report"
           # bracken
           "{results}/read_classification_results/krakenuniq_results/bracken_estimate_abundance/taxonomic_level_family/{sample}_family_abundance.bracken"
           "{results}/read_classification_results/krakenuniq_results/bracken_estimate_abundance/taxonomic_level_family/{sample}_bracken_family_abundance.report"
           "{results}/read_classification_results/krakenuniq_results/bracken_estimate_abundance/taxonomic_level_genus/{sample}_genus_abundance.bracken"
           "{results}/read_classification_results/krakenuniq_results/bracken_estimate_abundance/taxonomic_level_genus/{sample}_bracken_genus_abundance.report"
           "{results}/read_classification_results/krakenuniq_results/bracken_estimate_abundance/taxonomic_level_species/{sample}_species_abundance.bracken"
           "{results}/read_classification_results/krakenuniq_results/bracken_estimate_abundance/taxonomic_level_species/{sample}_bracken_species_abundance.report"
           "{results}/read_classification_results/krakenuniq_results/bracken_estimate_abundance/bracken_metaphlan_profile/merged_samples/merged_abundance_table_family.txt"
           "{results}/read_classification_results/krakenuniq_results/bracken_estimate_abundance/bracken_metaphlan_profile/merged_samples/merged_abundance_table_genus.txt"
           "{results}/read_classification_results/krakenuniq_results/bracken_estimate_abundance/bracken_metaphlan_profile/merged_samples/merged_abundance_table_species.txt"
           # centrifuge
           "{results}/read_classification_results/centrifuge_results/centrifuge_taxonomic/{sample}_centrifuge.tsv"
           "{results}/read_classification_results/centrifuge_results/centrifuge_taxonomic/{sample}_centrifuge_report.tsv"
           "{results}/read_classification_results/centrifuge_results/pavian_reports/{sample}_centrifuge_kreport.tsv"
           # metaphlan
           "{results}/read_classification_results/metaphlan_results/pavian_reports/{sample}_metaphlan_profile.txt"
           "{results}/read_classification_results/metaphlan_results/taxonomic_level_family/{sample}_metaphlan_profile_family.txt"
           "{results}/read_classification_results/metaphlan_results/taxonomic_level_genus/{sample}_metaphlan_profile_genus.txt"
           "{results}/read_classification_results/metaphlan_results/taxonomic_level_species/{sample}_metaphlan_profile_species.txt"
           "{results}/read_classification_results/metaphlan_results/taxonomic_level_sgb/{sample}_metaphlan_profile_sgb.txt"
           "{results}/read_classification_results/metaphlan_results/metaphlan_taxonomic/{sample}_bowtie.bt2"
           "{results}/read_classification_results/metaphlan_results/metaphlan_taxonomic/{sample}_bowtie_family.bt2"
           "{results}/read_classification_results/metaphlan_results/metaphlan_taxonomic/{sample}_bowtie_genus.bt2"
           "{results}/read_classification_results/metaphlan_results/metaphlan_taxonomic/{sample}_bowtie_species.bt2"
           "{results}/read_classification_results/metaphlan_results/metaphlan_taxonomic/{sample}_bowtie_sgb.bt2"
           "{results}/read_classification_results/metaphlan_results/taxonomic_level_gtdb/{sample}_metaphlan_profile_gtdb.txt"
           # kaiju
           "{results}/read_classification_results/kaiju_results/kaiju_taxonomic/{sample}_kaiju.tsv"
           # krona reads
           "{results}/read_classification_results/krakenuniq_results/krona_visualization/{sample}_krakenuniq.krona"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/krona_visualization/{sample}_krakenuniq_unclass.krona"
           "{results}/read_classification_results/centrifuge_results/krona_visualization/{sample}_centrifuge.krona"
           "{results}/read_classification_results/krakenuniq_results/krona_visualization/{sample}_krona_reads_krakenuniq.html"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/krona_visualization/{sample}_krona_reads_krakenuniq_unclass.html"
           "{results}/read_classification_results/centrifuge_results/krona_visualization/{sample}_krona_reads_centrifuge.html"
           "{results}/read_classification_results/metaphlan_results/krona_visualization/{sample}_metaphlan.krona"
           "{results}/read_classification_results/metaphlan_results/krona_visualization/{sample}_krona_reads_metaphlan.html"   
           "{results}/read_classification_results/kaiju_results/krona_visualization/{sample}_kaiju.krona"
           "{results}/read_classification_results/kaiju_results/krona_visualization/{sample}_krona_reads_kaiju.html"
           # recentrifuge
           "{results}/read_classification_results/krakenuniq_results/recentrifuge_visualization/{sample}_rcf_krakenuniq.html"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/recentrifuge_visualization/{sample}_rcf_krakenuniq_unclass.html"
           "{results}/read_classification_results/centrifuge_results/recentrifuge_visualization/{sample}_rcf_centrifuge.html"
           # krakentools
           "{results}/read_classification_results/krakenuniq_results/microbial_reads/{sample}_clean_R1.fastq"
           "{results}/read_classification_results/krakenuniq_results/microbial_reads/{sample}_clean_R2.fastq"
           # spades
           "{results}/contig_classification_results/assembly/spades/{sample}"
           "{results}/contig_classification_results/assembly/tmp/draft/{sample}/{sample}_spades_contigs.fasta"
           "{results}/contig_classification_results/assembly/tmp/draft/{sample}/{sample}_spades_scaffolds.fasta"
           "{results}/contig_classification_results/assembly/tmp/{sample}/{sample}_spades_contigs.fasta"
           "{results}/contig_classification_results/assembly/tmp/{sample}/{sample}_spades_scaffolds.fasta"       
           # metaspades
           "{results}/contig_classification_results/assembly/metaspades/{sample}"
           "{results}/contig_classification_results/assembly/tmp/draft/{sample}/{sample}_metaspades_contigs.fasta"
           "{results}/contig_classification_results/assembly/tmp/draft/{sample}/{sample}_metaspades_scaffolds.fasta"
           "{results}/contig_classification_results/assembly/tmp/{sample}/{sample}_metaspades_contigs.fasta"
           "{results}/contig_classification_results/assembly/tmp/{sample}/{sample}_metaspades_scaffolds.fasta"
           # metaviralspades
           "{results}/contig_classification_results/assembly/metaviralspades/{sample}"
           "{results}/contig_classification_results/assembly/tmp/draft/{sample}/{sample}_metaviralspades_contigs.fasta"
           "{results}/contig_classification_results/assembly/tmp/draft/{sample}/{sample}_metaviralspades_scaffolds.fasta"
           "{results}/contig_classification_results/assembly/tmp/{sample}/{sample}_metaviralspades_contigs.fasta"
           "{results}/contig_classification_results/assembly/tmp/{sample}/{sample}_metaviralspades_scaffolds.fasta"
           # megahit
           "{results}/contig_classification_results/assembly/megahit/{sample}"
           "{results}/contig_classification_results/assembly/tmp/draft/{sample}/{sample}_megahit_contigs.fasta"
           "{results}/contig_classification_results/assembly/tmp/draft/{sample}/{sample}_megahit_contigs_edit.fasta"
           "{results}/contig_classification_results/assembly/tmp/{sample}/{sample}_megahit_contigs.fasta"
           # cat assembly
           "{results}/contig_classification_results/assembly/tmp/combined_assembly/{sample}.fasta"
           # seqtk
           "{results}/contig_classification_results/assembly/short_filtered_assembly_threshold_{threshold1}/{sample}_threshold_{threshold1}.fasta"
           "{results}/contig_classification_results/assembly/long_filtered_assembly_threshold_{threshold2}/{sample}_threshold_{threshold2}.fasta"
           # viralverify
           "{results}/contig_classification_results/viralverify_classification_7_methods_threshold_{threshold1}/{sample}"
           # genomad
           "{results}/contig_classification_results/genomad_classification_7_methods_threshold_{threshold1}/{sample}"
           # diamond blastx
           "{results}/contig_classification_results/short_contigs_diamond_blastx_7_methods_threshold_{threshold1}/{sample}_diamond_blastx_contigs_7_methods_threshold_{threshold1}.daa"
           "{results}/contig_classification_results/long_contigs_diamond_blastx_7_methods_threshold_{threshold2}/{sample}_diamond_blastx_contigs_7_methods_threshold_{threshold2}.daa"
           # daa-meganizer
           "{results}/contig_classification_results/tmp/short_contigs_{sample}_diamond_blastx_7_methods_threshold_{threshold1}.log"
           "{results}/contig_classification_results/tmp/long_contigs_{sample}_diamond_blastx_7_methods_threshold_{threshold2}.log"
           # diamond view
           "{results}/contig_classification_results/short_contigs_diamond_view_7_methods_threshold_{threshold1}/{sample}_diamond_blastx_contigs_7_methods_threshold_{threshold1}.tab"
           "{results}/contig_classification_results/long_contigs_diamond_view_7_methods_threshold_{threshold2}/{sample}_diamond_blastx_contigs_7_methods_threshold_{threshold2}.tab"
           # krona contigs
           "{results}/contig_classification_results/short_contigs_krona_visualization_7_methods_threshold_{threshold1}/{sample}_krona_contigs_7_methods_threshold_{threshold1}.html"
           "{results}/contig_classification_results/long_contigs_krona_visualization_7_methods_threshold_{threshold2}/{sample}_krona_contigs_7_methods_threshold_{threshold2}.html"

Long reads workflow output structure:
           # cat raw fastq 
           "{results}/preprocess/catfastq/{sample}.fastq.gz"
           # nanoplot raw
           "{results}/preprocess/QC/reports_raw/{sample}"
           # nanocomp raw
           "{results}/preprocess/QC/combined_raw"
           # porechop + fastp
           "{results}/preprocess/trimmed/tmp/{sample}.fastq.gz"
           "{results}/preprocess/trimmed/{sample}_trim.fastq.gz"
           # nanoplot trim
           "{results}/preprocess/QC/reports_trim/{sample}"
           # nanocomp trim
           "{results}/preprocess/QC/combined_trim"
           # minimap2 + samtools
           "{results}/preprocess/host_depl/bam_files/{sample}_unal_sort.bam"
           "{results}/preprocess/host_depl/{sample}_clean.fastq.gz"
           # nanoplot clean
           "{results}/preprocess/QC/reports_clean/{sample}"
           # nanocamp clean
           "{results}/preprocess/QC/combined_clean"
           # krakenuniq
           "{results}/read_classification_results/krakenuniq_results/krakenuniq_taxonomic/{sample}_krakenuniq.krk"
           "{results}/read_classification_results/krakenuniq_results/pavian_reports/{sample}_krakenuniq.report"
           "{results}/read_classification_results/krakenuniq_results/unclassified_reads/{sample}_unclassified.fastq.gz"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/krakenuniq_taxonomic/{sample}_krakenuniq_unclass.krk"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/pavian_reports/{sample}_krakenuniq_unclass.report"
           # centrifuge
           "{results}/read_classification_results/centrifuge_results/centrifuge_taxonomic/{sample}_centrifuge.tsv"
           "{results}/read_classification_results/centrifuge_results/centrifuge_taxonomic/{sample}_centrifuge_report.tsv"
           "{results}/read_classification_results/centrifuge_results/pavian_reports/{sample}_centrifuge_kreport.tsv"
           # kaiju
           "{results}/read_classification_results/kaiju_results/kaiju_taxonomic/{sample}_kaiju.tsv"
           # krona reads
           "{results}/read_classification_results/krakenuniq_results/krona_visualization/{sample}_krakenuniq.krona"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/krona_visualization/{sample}_krakenuniq_unclass.krona"
           "{results}/read_classification_results/centrifuge_results/krona_visualization/{sample}_centrifuge.krona"
           "{results}/read_classification_results/krakenuniq_results/krona_visualization/{sample}_krona_reads_krakenuniq.html"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/krona_visualization/{sample}_krona_reads_krakenuniq_unclass.html"
           "{results}/read_classification_results/centrifuge_results/krona_visualization/{sample}_krona_reads_centrifuge.html"
           "{results}/read_classification_results/kaiju_results/krona_visualization/{sample}_kaiju.krona"
           "{results}/read_classification_results/kaiju_results/krona_visualization/{sample}_krona_reads_kaiju.html"
           # recentrifuge
           "{results}/read_classification_results/krakenuniq_results/recentrifuge_visualization/{sample}_rcf_krakenuniq.html"
           "{results}/read_classification_results/krakenuniq_results/check_unclassified/recentrifuge_visualization/{sample}_rcf_krakenuniq_unclass.html"
           "{results}/read_classification_results/centrifuge_results/recentrifuge_visualization/{sample}_rcf_centrifuge.html"
           # krakentools
           "{results}/read_classification_results/krakenuniq_results/microbial_reads/{sample}_clean.fastq"
           # metaflye
           "{results}/contig_classification_results/assembly/tmp/denovo_assembly/{sample}"
           "{results}/contig_classification_results/assembly/tmp/draft_assembly/{sample}.fasta"
           # medaka 
           "{results}/contig_classification_results/assembly/tmp/polishing/{sample}"
           "{results}/contig_classification_results/assembly/tmp/polished_assembly/{sample}.fasta"
           "{results}/contig_classification_results/assembly/final_assembly/{sample}.fasta"
           # viralverify
           "{results}/contig_classification_results/viralverify_classification/{sample}"
           # genomad
           "{results}/contig_classification_results/genomad_classification/{sample}"
           # diamond blastx
           "{results}/contig_classification_results/diamond_blast/{sample}_diamond_blast_contigs.daa"
           # daa-meganizer
           "{results}/contig_classification_results/meganizer_logs/{sample}.log"
           # diamond view
           "{results}/contig_classification_results/diamond_view/{sample}_diamond_blast_contigs.tab"
           # krona contigs
           "{results}/contig_classification_results/krona_visualization/{sample}_krona_contigs.html"
```
## List of tools used
[FastQC](https://github.com/s-andrews/FastQC)
[MultiQC](https://github.com/ewels/MultiQC)
[NanoPack](https://github.com/wdecoster/nanopack)
[Porechop](https://github.com/rrwick/Porechop)
[BBMap](https://github.com/BioInfoTools/BBMap)
[bowtie2](https://github.com/BenLangmead/bowtie2)
[minimap2](https://github.com/lh3/minimap2)
[samtools](https://github.com/samtools/samtools)
[SPAdes](https://github.com/ablab/spades)
[Flye](https://github.com/fenderglass/Flye)
[medaka](https://github.com/nanoporetech/medaka)
[seqtk](https://github.com/lh3/seqtk)
[KrakenUniq](https://github.com/fbreitwieser/krakenuniq)
[Krona](https://github.com/marbl/Krona)
[Pavian](https://github.com/fbreitwieser/pavian)
[viralVerify](https://github.com/ablab/viralVerify)
[DIAMOND](https://github.com/bbuchfink/diamond)
[MEGAN](https://github.com/husonlab/megan-ce)
