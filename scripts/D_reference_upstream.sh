#!/usr/bin/env bash
# Part D: clean and quantify the reference libraries (Pratx et al. 2025,
# Golden Promise WT, control N and heat shock T) exactly like ours
# (README Steps 2-3), with kallisto in unstranded mode.
# Run from the repository root:  bash scripts/D_reference_upstream.sh
# Needs fastp 1.3.6 and kallisto 0.52.0 on the PATH (or set FASTP / KALLISTO).
set -euo pipefail

FASTP=${FASTP:-fastp}
KALLISTO=${KALLISTO:-kallisto}
raw=data/external/pratx2025
clean=data/processed/pratx2025
out=results/D_validation/upstream
mkdir -p $clean $out/fastp/log $out/kallisto

for id in SRR35220185 SRR35220184 SRR35220183 SRR35220179 SRR35220178 SRR35220177; do
    echo "== ${id}"
    $FASTP \
        --in1 $raw/${id}_1.fastq.gz \
        --in2 $raw/${id}_2.fastq.gz \
        --out1 $clean/${id}_R1.fastp.fastq.gz \
        --out2 $clean/${id}_R2.fastp.fastq.gz \
        --json $out/fastp/${id}.fastp.json \
        --html $out/fastp/${id}.fastp.html \
        --thread 6 --detect_adapter_for_pe \
        -q 30 -u 15 -y -g -Y 10 -l 100 -b 150 -B 150 \
        2> $out/fastp/log/${id}.fastp.log
    mkdir -p $out/kallisto/${id}
    $KALLISTO quant --threads 12 --index references/MorexV3.idx \
        --gtf references/MorexV3.gtf \
        -o $out/kallisto/${id} \
        $clean/${id}_R1.fastp.fastq.gz $clean/${id}_R2.fastp.fastq.gz \
        2> $out/kallisto/${id}/kallisto_quant.log
done
