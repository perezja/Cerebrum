# Project Cerebrum

- Perform eQTL analysis and aggregate genotypic effects on expression using a meta-analysis
- This pipeline is adapted from GTEx (2017) 

# Docker on MGI

The pipeline environment is configured in the Docker image ```apollodorus/brain-eqtl:eqtl-analysis```. Pulling down the image and working interactively in a container can be done with the following job submission request on MGI:

```
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -Is -q research-hpc -J work-space -M 12000000 'select [mem>12000 && gtmp > 4] rusage[mem=12000, gtmp=4]' -a 'docker(apollodorus/brain-eqtl:eqtl-image)' /bin/bash
```

- the memory string using the ```-M``` option is in KB while the memory string in the ```rusage``` statement is in MB.

# Analysis

# I. eQTL-pipeline 

## 1. Prepare Expression

Gene counts (e.g.,, ```*.ReadsPerGene.out.tab``` from STAR) are merged into Gene Count Tables (gct) and Gene TPM Tables (gtt) for each cortex sub-region in a study  

- Environment: FENIX
- Inputs:

  - ```$COUNT_JSON```: The input JSON holds paths to the head directory where individual sample-level STAR or Salmon outputs are stored. An example JSON is provided in ```/gscmnt/gc2645/wgs/qtl/BrainQTL/Cerebrum/json/count_tpm_files.json``` 
  - ```$GTF```: GTF annotation file corresponding to genome version used in alignment (e.g., ```/40/pipelines/RNAseq/STAR/Hg19_gencodev19_spikein/gencode.v19.annotation.spike-in.gtf```)

- For combining gene counts

```
python3 combine_expression.py $COUNT_JSON $GTF --mode cts -o $GCT_DIR 
```

- And combining transcript level counts (e.g., from Salmon) 

```
python3 combine_expression.py $COUNT_JSON $GTF --mode tpm -o $GTT_DIR 
```

- Outputs:

  - Region-level (gene x sample) expression matrices with matrix dimension as header 
  - see ```data/individual_level_stats/gct``` and  ```data/individual_level_stats/gtt```

- Estimated compile time
  - < 1 min

## 2. Expression Normalization

Gene expression results are subset with count and tpm thresholds and 'TMM' normalization is applied to get normalized expression matrices 

- Environment: FENIX

- Inputs:

  - GCT and GTT matrices present in ```GCT_DIR``` and ```GTT_DIR``` directories, generating ```.bed``` and ```.txt``` files 
  - ```$SAMPLE2PARTICIPANT```:  

```
python3 normalize_expression.py $GCT_DIR $GTT_DIR $SAMPLE2PARTICIPANT --annotation_gtf $COLLAPSED_GTF 
```

```$COLLAPSED_GTF``` must be collapsed with [this script](https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py) by Francois Auguet

## 3. Prepare Covariates

PEER factors can be generated using ```run_PEER.R```

- Environment: Docker (apollodorus/brain-eqtl:eqtl-analysis) 

- Inputs:

  - $EXPR: Either the ```.bed``` or ```.gz``` expression matrix from step 2 
  - $PREFIX: Outfile prefix
  - $N: Number of hidden covariates to estimate 
  - $OUTPUT_DIR: Output directory

```
Rscript run_PEER.R $EXPR $PREFIX $N -o $OUTPUT_DUR 
```

- Outputs:

  - $PREFIX.PEER_alpha.txt
  - $PREFIX.PEER_covariates.txt
  - $PREFIX.PEER_residuals.txt

- Estimated compile time: < 5m

Individual phenotypes are combined with PEER expression covariates to create covariate tables for each sub-region in a study. 

- Environment: Docker (apollodorus/brain-eqtl:eqtl-analysis) 

- Inputs:

    - $PHENOTYPES: A table with sample names in column 1 followed by their corresponding k phenotypes in columns (2,...,k+1) 

participant_id  PC1     PC2     PC3     Age     Gender
DIAN_1DKYRE^8012957511^201711_CoreExome -0.00596815     -0.00293676     -0.00766941     57.0    0
DIAN_2968OM^8012957462^2017NeuroX2      -0.00665722     0.0008246369999999999   -0.00591766     47.0    0
DIAN_62CYZP^8012957085^201711_CoreExome -0.00551273     -0.00310693     0.00270309      63.0    1
DIAN_9TPSKM^8012957707^201711_CoreExome -0.00550265     -0.00264916     0.0110613       40.0    0

    - $SAMPLE2REGION: TSV lookup table mapping participants (sample names) to a region ID (same as used for naming expression matrices and peer covariate files) 

MayoADGS_786^unk^2017Omni25Exome        TCTX
MayoADGS_787^unk^2017Omni25Exome        TCTX
GTEX_ZUA1^unk^2017_GTEXv7_HiSeqWGS      PFCTX
MAP_63377^09AD19850^2013OmniEx_NACCR3   PAR

    - $EXPR_DIR: Directory with expression matrices

```
python3 combine_covariates.py $PHENOTYPES $SAMPLE2REGION $EXPR_DIR 
```

- Outputs:

  - $REGION.combined_covariates.txt: A covariate matrix (covariates x sample) 

- Estimated compile time < 1min

It is CRUCIAL that the order of sample IDs columns in ```$REGION.combined_covariates.txt``` file match order of samples in expression matrices. The script ```order_participants.R``` in ```src/prepare_individual_data/prepare_covariates``` reorders the former according to latter. 

## 4. Prepare genotypes

Assuming genotypes for participants are already imputed and passed QC, some data manipulation is necessary in Plink to:
  - subset a merged genotype bfile for region/study specific participants
    - ```extract_ids.sh``` is an sample script of how to do this simply with bash. 
  - recode variants into minor allele dosage counts
  - transpose the matrix into (SNP x participant_id) orientation
  - drop some non-id columns

```
./extract_ids.sh region_participant_dir
```
- writes region/study id lists ```$REGION_id_list.txt``` to directory ```region_participant_dir```

- Environment: FENIX

- Inputs:

  - $PARTICIPANT_LIST_DIR: Directory with participant ID lists for every region/study 
  - $OUTPUT_DIR: Output directory 

```
prepare_genotypes.py $BFILE $PARTICIPANT_LIST_DIR --output_dir $OUTPUT_DIR 
```

- Outputs:

  - genotype file ```<region_id>.snps.txt``` is formatted for qtl analysis

## 4. Region-level regressions of expression on genotype 

MatrixQTL is used to calculate cis and trans eQTLs in cortex regions per study

- Environment: MGI on virtual-workstation. Uses Docker image (apollodorus/brain-eqtl:eqtl-analysis) 

- Inputs:

  - $JSON: JSON file with paths to expression, genotype and covariate data. Look at ```sample_data/mqtl_inputs.json``` 
  - $CHUNKSIZE: Number of chunks (e.g., 500) to split expression matrices into for parallel processing
  - $OUTPUT_DIR: Output directory

```
python3 qtl_analysis.py paths.json $CHUNKSIZE --lsf -o $OUTPUT_DIR
```

The ```--lsf``` flag wraps MatrixQTL call in a LSF job submission command to run on MGI

- Outputs: 
 
  - $REGION_chunk$N.assoc.txt.gz 

SNP     gene    beta    t-stat  p-value FDR
12:31226835:A:T ENSG00000013573 1.2170867316562 22.5695254857688        2.66679285555292e-33    6.91028431826582e-24
12:31244846:C:G ENSG00000013573 1.17756467048307        19.9822594463213        3.60187054210892e-30    2.40499185917088e-22
12:31251544:C:T ENSG00000013573 1.17756467048307        19.9822594463213        3.60187054210892e-30    2.40499185917088e-22
12:31251803:A:G ENSG00000013573 1.17756467048307        19.9822594463213        3.60187054210892e-30    2.40499185917088e-22

- Estimated compile time:
    - 4-6 hrs

## 5. Meta-analysis

Chunk-level Estimates of QTL associations are aggregated gene-wise to prepare input for Metasoft. 

- Environment: MGI in a virtual-workstation 

I. Fragment chunks into gene-level association files ```fragment_to_gene_chunks.py```

- Inputs:

  - $CHUNK_DIR: Directory with MatrixQTL association chunk data 
  - $PREFIX: Prefix (region or study) to use as head directory when writing gene-level chunk files 

- Outputs:

  - MatrixQTL gene-wise chunkized files prefixed with gene name under ```$OUTPUT_DIR/$PREFIX``` directory 

- Estimated compile time:

  - 2-N hrs Depends on size of input. Periodically check job status on LSF

II. Merge gene-level chunks across studies into METASOFT input format:

```
rsAAAAAA study1beta study1stderr study2beta study2stderr study3beta study3stderr
rsBBBBBB study1beta study1stderr study2beta study2stderr study3beta study3stderr
rsCCCCCC study1beta study1stderr study2beta study2stderr study3beta study3stderr
```

- Inputs:

  - $GENE_CHUNKDIR: Head directory (```$OUTPUT_DIR``` in step 5I) 
  - $ID2SYMBOL: TSV file with ENSMEBL IDs in column1 and gene names in column 2. ONLY genes listed in this file will be merged. 
  - $OUTOUT_DIR: Output directory 

```
python3 build_metasoft.py $GENE_CHUNKDIR $ID2SYMBOL -o $OUTPUT_DIR 
```

- Outputs:

  - METASOFT input chunks ```ms_chunk_$CHUNKID.txt.gz``` 

- Estimated compile time:

  - 2-N hrs Depends on size of input. Periodically check job status on LSF

III. Run METASOFT 

- Inputs: 

  - $MS_INPUT_DIR: METASOFT chunk input directory (```$OUTPUT_DIR``` from step 5II)
  - $OUTPUT_DIR: Output directory

```
python3 parallel_run_metasoft.py $MS_INPUT_DIR -o $OUTPUT_DIR
```

- Outputs:

  - METASOFT chunk outputs ```$OUTPUT_DIR/$CHUNKID.metasoft.txt```

- Estimated compile time:

  - 2-N hrs Depends on size of input. Periodically check job status on LSF. You will likely run this over and over in waves until all genes are successfully processed

Downstream post-processing (concatenation of METASOFT chunks) is left to the user
