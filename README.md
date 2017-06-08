# simba-publication-pipeline

Pipeline used to generate SimBa publication results.

To use the pipeline you will need to modify the following variables :
- **GENOME_DIR**
- **REFERENCE_BASENAME**
- **ANNOTATIONS**
- **POLYMORPHISMS** 
- **FLUX_ERROR_MODEL**
- **TMP_DIR**

You will also need to uncomment the simct rule that will generate simulated dataset. 
The rule simct causes problem with the rest of the pipeline due to complex use of Snakefile wildcards.
