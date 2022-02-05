# [Nextflow Tutorial - Variant Calling Edition](https://sateeshperi.github.io/nextflow_varcal/)

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/sateeshperi/nextflow_tutorial.git)

Once the pod launches, it will present a VS-Code interface

* To create `varcal` environment based on yml file
```bash
bash envconfig.sh
```
>If conda is not available in PATH:
>
>```bash
>conda init bash
>```
>
>```bash
>source ~/.bashrc
>```
>**Close terminal and open a new terminal. You should be able to see `(base)` at the beginning of terminal indicating active conda base environment.**


* To download reference genome and raw reads
```bash
bash data/fetch_raw_data.sh
```

* To trim the reads using trimmomatic
```bash
conda activate varcal
```
```bash
bash data/trim.sh
```

---

>When executing nextflow if you see `Unable to initialize nextflow environment` error:
>
>```bash
>unset JAVA_TOOL_OPTIONS
>export NFX_OPTS=$JAVA_TOOL_OPTIONS
>```

---
