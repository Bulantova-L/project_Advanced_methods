version 1.0

workflow haplotype_methylation_workflow {
  input {
    File gene_bed
    File maternal_methylation_bed
    File paternal_methylation_bed
    File summary_script
  }

  # Intersect for maternal
  call bedtools_intersect as intersect_maternal {
    input:
      gene_bed = gene_bed,
      methylation_bed = maternal_methylation_bed,
      output_name = "genes_maternal_overlap.txt"
  }

  # Intersect for paternal
  call bedtools_intersect as intersect_paternal {
    input:
      gene_bed = gene_bed,
      methylation_bed = paternal_methylation_bed,
      output_name = "genes_paternal_overlap.txt"
  }

  # Run summarization Python script
  call summarize_methylation {
    input:
      maternal_overlap = intersect_maternal.overlap_file,
      paternal_overlap = intersect_paternal.overlap_file,
      script = summary_script
  }

  output {
    File final_summary = summarize_methylation.summary
  }
}


task bedtools_intersect {
  input {
    File gene_bed          # gene coordinates (e.g., GENCODE in BED)
    File methylation_bed   # haplotype-specific methylation BED
    String output_name     # output file name
  }

  command <<<
    bedtools intersect -a ${gene_bed} -b ${methylation_bed} -wa -wb > ${output_name}
  >>>

  output {
    File overlap_file = "${output_name}"
  }

  runtime {
    docker: "quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_3"
  }
}
task summarize_methylation {
  input {
    File maternal_file
    File paternal_file
    String output_name = "haplotype_methylation_summary.csv"
    File script
  }

  command <<<
    python3 ${script} \
      --maternal ${maternal_file} \
      --paternal ${paternal_file} \
      --output ${output_name}
  >>>

  output {
    File summary = "${output_name}"
  }

  runtime {
    docker: "python:3.10"
  }
}
