rule align_model:
    input:
        pos = "refs/positives/{model}.faa"
    output:
        aln = "refs/positives/{model}.aln"
    threads: 8
    conda: "../../envs/base.yml"
    shell:
        r"""
        mafft --maxiterate 1000 --localpair --thread {threads} {input.pos} > {output.aln} 2> logs/mafft.{wildcards.model}.log
        """

rule trim_model:
    input:
        aln = rules.align_model.output.aln
    output:
        trimaln = "refs/positives/{model}.trim.aln"
    conda: "../../envs/base.yml"
    shell:
        r"""
        trimal -in {input.aln} -out {output.trimaln} -automated1 2> logs/trimal.{wildcards.model}.log
        """

rule hmmbuild_model:
    input:
        trimaln = rules.trim_model.output.trimaln
    output:
        hmm = "refs/hmm/{model}.hmm"
    conda: "../../envs/base.yml"
    shell:
        r"""
        hmmbuild {output.hmm} {input.trimaln} > logs/hmmbuild.{wildcards.model}.log 2>&1
        """

