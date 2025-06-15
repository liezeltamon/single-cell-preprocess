#rule qc_sc:
#    input:
#        parent_counts_dir = lambda wc: os.path.dirname(DATA_DIR),
#        group_id_dir = lambda wc: os.path.join(RESULTS_DIR, "dehash"),
#        whitelist_dir = lambda wc: os.path.join(RESULTS_DIR, "filter")
#    output:
#       agg_heatmap = os.path.join(RESULTS_DIR, "qc_sc", "aggregate", "heatmaps.pdf"),
#       agg_metrics = os.path.join(RESULTS_DIR, "qc_sc", "aggregate", "metrics.csv"),
#        metrics = os.path.join(RESULTS_DIR, "qc_sc", "{sample}", "metrics.csv"),
#        session = os.path.join(RESULTS_DIR, "qc_sc", "{sample}", "plots.pdf")
#    params:
#        main = lambda wc: params_to_cli_args(config["qc_sc"]),
#        parent_out_dir = lambda wc: os.path.join(RESULTS_DIR, "qc_sc"),
#        out_dir = lambda wc: os.path.join(RESULTS_DIR, "qc_sc", wc.sample)
#    log:
#        os.path.join(RESULTS_DIR, "qc_sc", "{sample}", "log.out")
#    shell:
#        """
#        mkdir -p {params.parent_out_dir} && \
#        Rscript scripts/{rule}.R \
#            --parent_counts_dir {input.parent_counts_dir} \
#            --group_id_dir {input.group_id_dir} \
#            --parent_out_dir {params.parent_out_dir} \
#            --whitelist_dir {input.whitelist_dir} \
#            {params.main} &> {log}
#        """

rule qc_sc_sample:
    input:
        parent_counts_dir = lambda wc: os.path.dirname(DATA_DIR),
        group_id_dir = lambda wc: os.path.join(RESULTS_DIR, "dehash"),
        whitelist_dir = lambda wc: os.path.join(RESULTS_DIR, "filter")
    output:
        metrics = os.path.join(RESULTS_DIR, "qc_sc", "{sample}", "metrics.csv"),
        session = os.path.join(RESULTS_DIR, "qc_sc", "{sample}", "plots.pdf")
    params:
        main = lambda wc: params_to_cli_args(config["qc_sc"]),
        parent_out_dir = lambda wc: os.path.join(RESULTS_DIR, "qc_sc"),
        out_dir = lambda wc: os.path.join(RESULTS_DIR, "qc_sc", wc.sample)
    log:
        os.path.join(RESULTS_DIR, "qc_sc", "{sample}", "log.out")
    shell:
        """
        mkdir -p {params.out_dir} && \
        Rscript scripts/qc_sc.R \
            --parent_counts_dir {input.parent_counts_dir} \
            --group_id_dir {input.group_id_dir} \
            --parent_out_dir {params.parent_out_dir} \
            --whitelist_dir {input.whitelist_dir} \
            {params.main} &> {log}
        """

rule qc_sc_aggregate:
    input:
        metrics = expand(os.path.join(RESULTS_DIR, "qc_sc", "{sample}", "metrics.csv"), sample=SAMPLES)
    output:
        agg_heatmap = os.path.join(RESULTS_DIR, "qc_sc", "aggregate", "heatmaps.pdf"),
        agg_metrics = os.path.join(RESULTS_DIR, "qc_sc", "aggregate", "metrics.csv")
    params:
        main = lambda wc: params_to_cli_args(config["qc_sc"]),
        parent_out_dir = lambda wc: os.path.join(RESULTS_DIR, "qc_sc"),
    log:
        os.path.join(RESULTS_DIR, "qc_sc", "aggregate", "log.out")
    shell:
        """
        mkdir -p {params.parent_out_dir}/aggregate && \
        Rscript scripts/qc_sc.R \
            --parent_counts_dir {params.parent_out_dir} \
            --group_id_dir {params.parent_out_dir}/aggregate \
            --parent_out_dir {params.parent_out_dir} \
            --aggregate_mode TRUE \
            {params.main} &> {log}
        """
