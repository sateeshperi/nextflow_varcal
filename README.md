# [Nextflow Tutorial - Variant Calling Edition](https://sateeshperi.github.io/nextflow_varcal/)

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/sateeshperi/nextflow_tutorial.git)

* **The GitPod comes with Nextflow, Conda and Docker pre-installed**

Once the pod launches, it will present a VS-Code interface

* To create `varcal` environment based on yml file
```bash
mamba env update -n base -f environment.yml
```

* To download reference genome and raw reads
```bash
bash data/fetch_raw_data.sh
```

* To trim the reads using trimmomatic
```bash
bash data/trim.sh
```

---
