#######################################
# CONFIG
#######################################
# Input data files we need - downloaded in rule download_data
INPUT_METHYLATION_BED = "data/Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC.bed"
MATERNAL_GENE_BED     = "data/maternal.bed"
PATERNAL_GENE_BED     = "data/paternal.bed"
FILTER_CPG = 10   # minimum number of CpGs per gene to keep

PARENTS = ["maternal", "paternal"]
TYPES = ["methylation"]

#######################################
# RULE ALL
#######################################

rule all:
    input:
        INPUT_METHYLATION_BED,
        MATERNAL_GENE_BED,
        PATERNAL_GENE_BED,
        #expand("results/intersect/{parent}_overlap.bed", parent=PARENTS),
        expand("results/{type}/{parent}_gene_{type}.bed", parent=PARENTS, type=TYPES),
        expand("results/stats/{type}_stats.tsv", type=TYPES),
        expand("results/summary/{type}_summary.tsv", type=TYPES),
        expand("results/stats/{type}_stats_top100.tsv", type=TYPES)

rule download_methylation_data:
    output:
        methylation_bed=INPUT_METHYLATION_BED
    shell:
        """
        wget -O {output.methylation_bed} https://public.gi.ucsc.edu/~mcechova/HG002/Q100_ONT_5mC_HG002v1.1_winnowmap_q10_10kb_modkit5mC.bed
        """
rule download_annotation_data:
    output:
        maternal_bed=MATERNAL_GENE_BED,
        paternal_bed=PATERNAL_GENE_BED
    shell:
        """
        wget -O {output.maternal_bed} https://is.muni.cz/www/bulantova.l/pv269_project/maternal.bed?lang=en;stahnout=1;dk=rRwhcRTd
        wget -O {output.paternal_bed} https://is.muni.cz/www/bulantova.l/pv269_project/paternal.bed?lang=en;stahnout=1;dk=T2qlS5XN
        """
    
#######################################
# SPLIT HAPLOTYPES
#######################################

rule split_haplotypes_and_sort:
    input:
        bed=INPUT_METHYLATION_BED
    output:
        maternal=temp("results/split/maternal_methylation.bed"),
        maternal_sorted="results/split/maternal_methylation_sorted.bed",
        paternal=temp("results/split/paternal_methylation.bed"),
        paternal_sorted="results/split/paternal_methylation_sorted.bed"
    conda:
        "envs/bedtools.yml"
    threads: 8
    shell:
        """
        # SPLIT without pipes
        grep "MATERNAL" {input.bed} > {output.maternal}
        grep "PATERNAL" {input.bed} > {output.paternal}

        # SORT
        bedtools sort -i {output.maternal} > {output.maternal_sorted}
        bedtools sort -i {output.paternal} > {output.paternal_sorted}
        """

# EXTRACT GENES COORDINATES

rule extract_genes_coordinates:
    input:
        "data/{parent}.bed"
    output:
        "data/{parent}_genes.bed"
    conda:
        "envs/bedtools.yml"
    shell:
        """
        cut -f1-4 {input} \
        | bedtools sort -i - \
        > {output}
        """

# MEAN METHYLATION PER GENE

rule map_methylation:
    input:
        genes="data/{parent}_genes.bed",
        methylation="results/split/{parent}_methylation_sorted.bed"
    output:
        "results/methylation/{parent}_gene_methylation.bed"
    conda:
        "envs/bedtools.yml"
    shell:
        """
        mkdir -p results/mapped 
        bedtools map \
            -a {input.genes} \
            -b {input.methylation} \
            -c 11 -o mean,count \
            > {output}
        """


rule filter_low_cpg:
    input:
        "results/methylation/{parent}_gene_methylation.bed"
    output:
        "results/methylation/{parent}_gene_methylation_filtered.bed"
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"}} $6 > 10 {{print $0}}' {input} > {output}
        """

# MERGE TABLES both methylation and normalized methylation
rule merge_methylation_tables:
    input:
        maternal="results/{type}/maternal_gene_{type}_filtered.bed",
        paternal="results/{type}/paternal_gene_{type}_filtered.bed"
    output:
        "results/summary/{type}_summary.tsv"
    shell:
        r"""
        awk 'BEGIN {{ OFS="\t" }}
            # Load maternal: save methylation + count
            NR==FNR {{
                mat_meth[$4] = $5;
                mat_count[$4] = $6;
                next
            }}
            # Process paternal
            {{
                pat_meth = $5;
                pat_count = $6;
                gene = $4;

                if (gene in mat_meth)
                    print gene, mat_meth[gene], pat_meth, mat_count[gene], pat_count;
            }}' {input.maternal} {input.paternal} \
        | sed '1iGene_ID\tMaternal_Methylation\tPaternal_Methylation\tMaternal_Count\tPaternal_Count' \
        > {output}
        """

rule stats_and_plots:
    input:
        summary="results/summary/{type}_summary.tsv"
    output:
        results_table="results/stats/{type}_stats.tsv",
        results_xlsx = "results/stats/{type}_stats.xlsx",
        results_table_top100="results/stats/{type}_stats_top100.tsv",
        results_xlsx_top100 = "results/stats/{type}_stats_top100.xlsx",
        histogram="results/plots/{type}_histogram.png",
        scatter="results/plots/{type}_scatter.png",
        boxplot="results/plots/{type}_boxplot.png",
        pdf="results/plots/{type}_combined_plots.pdf",
        piechart="results/plots/{type}_piechart.png"

    conda:
        "envs/stats.yml"       # must contain pandas, scipy, statsmodels, seaborn, matplotlib
    script:
        "scripts/stats_and_plots.py"

