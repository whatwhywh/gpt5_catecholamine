configfile: "../config/config.yaml"

rule build_backbone_tree:
    input:
        align = "refs/positives/backbone.aln"
    output:
        tree = "results/phylogeny/aaad_backbone_tree.nwk"
    threads: 8
    conda: "../envs/base.yml"
    shell:
        """
        if command -v iqtree2 >/dev/null 2>&1; then \
            iqtree2 -s {input.align} -fast -T AUTO -nt {threads} -m MFP -pre results/phylogeny/backbone; \
            cp results/phylogeny/backbone.treefile {output.tree}; \
        else \
            FastTree -gamma -wag < {input.align} > {output.tree}; \
        fi
        """

