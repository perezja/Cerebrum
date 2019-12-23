# Colocalization of GWAS loci

## Step 1. The indexSNP file is a tab-delimiated file with five columns:

CHR: chromosome
SNP: SNP ID of index SNP
BP: bp position of index SNP
STARTBP: start of interval around index SNP (bp)
ENDBP: end of interval around index SNP (bp)

- The inputs are the following:

GWAS file:

- Summary statistics of leading GWAS SNPs

SNP	RSID	Chr	Bp	Pvalue	Trait
21:44333234	rs75087725	21	44333234	3e-10	Amyotrophic lateral sclerosis
3:39492990	rs616147	3	39492990	4e-10	Amyotrophic lateral sclerosis
14:30678292	rs10139154	14	30678292	3e-08	Amyotrophic lateral sclerosis
12:64488187	rs74654358	12	64488187	7e-08	Amyotrophic lateral sclerosis
9:27543384	rs3849943	9	27543384	4e-19	Amyotrophic lateral sclerosis

BED file:

- Coordinates of gene TSS in BED format (strand is not used) 

1	11869	11870	ENSG00000223972	1
1	24886	24887	ENSG00000227232	-1
1	29554	29555	ENSG00000243485	1
1	36073	36074	ENSG00000237613	-1
1	52473	52474	ENSG00000268020	1

eQTL summary statistics:

- only cis eQTLs are used

13:90103810:A:G	ENSG00000184371	CSF1	TRANS	A	G	0.3272	0.095911	0.017577799999999998	4.8593699999999995e-08
10:44354068:C:A	ENSG00000184371	CSF1	TRANS	A	C	0.05446	0.196798	0.0358699	4.1007699999999996e-08
10:44347103:A:G	ENSG00000184371	CSF1	TRANS	A	G	0.05446	0.196798	0.0358699	4.1007699999999996e-08
10:44360748:C:T	ENSG00000184371	CSF1	TRANS	C	T	0.05446	0.196798	0.0358699	4.1007699999999996e-08

```
python3 coloc_index_snps.py $GWAS $BED $EQTL -o $OUTDIR 
```

will generate the index file ```$(PREFIX).index```

# Step 2. Intersect eQTL associations with loci's

- Intersect cis-eQTLs falling within GWAS loci and generate candidate eQTL locus files

$CIS: METASOFT file holding all cis associations
$INDEX: GWAS locus intervals created above
$OUT: Output directory for eQTL locus files
$PREFIX: Outfile prefices

```
CIS=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/eqtl_discovery/metasoft/cortex_region.cis.metasoft.txt
INDEX=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/coloc/neurodeg_loci.index
OUT=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/coloc/locus_eqtl_files
PREFIX=cortex_region
python3 eqtl_loci_segmentation.py $CIS $INDEX $PREFIX -o $OUT
```

# Step 3. Run colocalisation routine

- $OUTCOME_MAP: TSV file mapping an outcome ID in first column to file path in second

```
/gscmnt/gc2645/wgs/qtl/BrainQTL/src/new_src/coloc/run_coloc.R /gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/coloc/locus_eqtl_files/cortex_region.17.59817366.60017366 /gscmnt/gc2645/wgs/qtl/BrainQTL/data/genome/plink/maf/cortex_regions.maf.frq /gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/outcome_data/AD_kunkle.outcome.txt AD -o /gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/coloc/results/cortex_region.17.59817366.60017366_AD.coloc
```

```
LOCUS_DIR=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/coloc/locus_eqtl_files
FREQ=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/genome/plink/maf/cortex_regions.maf.frq
OUTDIR=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/coloc/results/
OUTCOME_MAP=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/outcome_data/outcome_data.txt
python3 coloc_run_parallel.py $LOCUS_DIR $FREQ $OUTCOME_MAP -o $OUTDIR
```

# Step 4. Aggregate results

```
OUTFILE=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/coloc/neurodeg_coloc.txt
RESDIR=/gscmnt/gc2645/wgs/qtl/BrainQTL/data/gwas_colocalization/coloc/results/*.coloc
echo -e 'locus_name\tsentinel_eqtl\tsentinel_eQTL_gene\ttrait\tnsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf' > $OUTFILE 
find $RESDIR -maxdepth 1 -type f | sort -V | xargs -I{} tail -n+2 {} >> $OUTFILE
```
