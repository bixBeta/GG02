# GG02
development envir for gg02 biohpc


# beta6.sh <- RNA-seq workflow
dependencies:
  - trim_galore: version 0.6 or above (should be in path on [biohpc](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=663#c))
  - multiqc: install multiqc using the following command:
      `pip install --user multiqc`  
  - STAR: version 2.7.0e or above (add to path `export PATH=/programs/STAR:$PATH`)
  - RSeQC: version 2.61 (add to path using this guideline [biohpc](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=135#c))
