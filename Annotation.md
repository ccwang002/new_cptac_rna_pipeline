## Build RSEM index

```
docker pull quay.io/biocontainers/rsem:1.3.1--pl526r341h4f16992_0

docker run -u (id -u):(id -g) \
    -v /diskmnt/Datasets/Reference:/genome_ref \
    -t \
    quay.io/biocontainers/rsem:1.3.1--pl526r341h4f16992_0 \
    rsem-prepare-reference \
        --gtf /genome_ref/GDC/gencode.v22.annotation.gtf \
        --num-threads 8 \
        /genome_ref/GRCh38.d1.vd1/GRCh38.d1.vd1.fa \
        /genome_ref/GDC/rsem/rsem_genome_d1_vd1_gtfv22 \
    2> rsem/rsem_v1.3.1_prepare_ref.log 1>&2
```