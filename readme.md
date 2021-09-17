# GATK4 haplotypecaller comparitive pipeline

* written in WDL and Nextflow (DSL2)
* utilizes the identical docker containers
* utilizes the identical bash code

WDL pipeline is `WDL/haplotypecaller.wdl`.\
Nextflow (DSL2) pipeline is `main.nf` and `nextflow.config`.

*The current version does not test scatter functionality*

Current comparison:
|Step|wdl|nf_dsl2|
|:---|:---:|:---:|
|common|<1m|-|
|mapping|2m|2m|
|mardup|2m|1m|
|fastaIndex|4m|4m|
|recal|8m|2m|
|bamIndex|1m|<1m|
|haplotypecaller|43m|43m|
|genotype|7m|8m|
|output|<1m|-|
|Process Time|1h 9m|1h 1m|
|Wall Time|1h 10m|1h 2m|
|Estimated Cost|$0.81|$0.86| 