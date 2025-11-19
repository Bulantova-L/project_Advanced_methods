version 1.0

workflow haplotype_methylation_workflow {
  input {
    File gene_bed_maternal
    File gene_bed_paternal
    File full_methylation_bed
    File summary_script
  }

  # Split maternal / paternal methylation
  call split_haplotypes {
    input:
      full_methylation_bed = full_methylation_bed,
      maternal_name = "maternal_methylation",
      paternal_name = "paternal_methylation"
  }

  # Intersect for maternal
  call bedtools_intersect as intersect_maternal {
    input:
      gene_bed = gene_bed_maternal,
      methylation_bed = split_haplotypes.maternal_bed,
      output_name = "genes_maternal_overlap"
  }

  # Intersect for paternal
  call bedtools_intersect as intersect_paternal {
    input:
      gene_bed = gene_bed_paternal,
      methylation_bed = split_haplotypes.paternal_bed,
      output_name = "genes_paternal_overlap"
  }

  # Run summarization Python script
  call summarize_methylation {
    input:
      maternal_overlap = intersect_maternal.overlap_file,
      paternal_overlap = intersect_paternal.overlap_file,
      script = summary_script,
      output_name = "gene_methylation_summary"
  }

  output {
    File final_summary = summarize_methylation.summary
  }
}

task split_haplotypes {
  input {
    File full_methylation_bed
    String maternal_name 
    String paternal_name
  }

  command <<<
    set -euo pipefail
    grep "MATERNAL" ~{full_methylation_bed} > ~{maternal_name}.bed
    grep "PATERNAL" ~{full_methylation_bed} > ~{paternal_name}.bed
  >>>

  output {
    File maternal_bed = "~{maternal_name}.bed"
    File paternal_bed = "~{paternal_name}.bed"
  }

  runtime {
    docker: "quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_3"   # or any base image with grep
  }
}


task bedtools_intersect {
  input {
    File gene_bed          # gene coordinates (e.g., GENCODE in BED)
    File methylation_bed   # haplotype-specific methylation BED
    String output_name     # output file name
  }

  command <<<
    bedtools intersect -a ~{gene_bed} -b ~{methylation_bed} -wa -wb > ~{output_name}.bed
  >>>

  output {
    File overlap_file = "~{output_name}.bed"
  }

  runtime {
    docker: "quay.io/biocontainers/bedtools:2.31.0--hf5e1c6e_3"
  }
}

task summarize_methylation {
  input {
    File maternal_overlap
    File paternal_overlap
    String output_name
    File script
  }

  command <<<
    python3 ${script} \
      --maternal ${maternal_overlap} \
      --paternal ${paternal_overlap} \
      --output ${output_name}.csv
  >>>

  output {
    File summary = "~{output_name}.csv"
  }

  runtime {
    docker: "python:3.10"
  }
}
