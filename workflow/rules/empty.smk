import os

# Parameters

PROJECT_DIR = os.path.abspath(os.path.join(os.getcwd(), "../../"))
DATA_DIR = os.path.join(PROJECT_DIR, "data", "cellranger")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "process_droplets")

configfile: os.path.join("envs", "config.yml")

SAMPLES = [d for d in os.listdir(DATA_DIR) if os.path.isdir(os.path.join(DATA_DIR, d))]

# Functions

def params_to_cli_args(param_dict):
    args = []
    for k, v in param_dict.items():
        flag = f"--{k.replace('_', '-')}"
        if isinstance(v, bool):
            if v:
                args.append(flag)
        else:
            args.append(f"{flag} {v}")
    return " ".join(args)

def params_to_string(params_dict, parent=None):
    parts = []
    for k in sorted(params_dict):
        v = params_dict[k]
        if isinstance(v, float):
            v_str = f"{v:.4f}".replace('.', '')
        else:
            v_str = str(v)
        parts.append(f"{k}{v_str}")
    
    param_str = "_".join(parts)
    
    if parent:
        return f"{parent}/{param_str}"
    else:
        return param_str

# Rules

rule all:
  input:
      expand(
          os.path.join(RESULTS_DIR, "{param_set}", "{sample}", "{filename}"),
          sample=SAMPLES,
          filename=[
              "blacklist.txt",
              "whitelist.txt",
              "output.qs",
              "plots.pdf",
              "session_info.txt"
          ]
      )
rule empty_emptydrops:
    input:
        counts_dir = lambda wc: os.path.join(DATA_DIR, wc.sample),
    output:
        blacklist = lambda wc: os.path.join(RESULTS_DIR, wc.param_set, wc.sample, "blacklist.txt"),
        whitelist = lambda wc: os.path.join(RESULTS_DIR, wc.param_set, wc.sample, "whitelist.txt"),
        output = lambda wc: os.path.join(RESULTS_DIR, wc.param_set, wc.sample, "output.qs"),
        plot = lambda wc: os.path.join(RESULTS_DIR, wc.param_set, wc.sample, "plots.pdf"),
        session = lambda wc: os.path.join(RESULTS_DIR, wc.param_set, wc.sample, "session_info.txt"),
    params:
        param_set = lambda wildcards, snakemake: params_to_string(config[snakemake.rule.name]),
        main = lambda wildcards, snakemake: params_to_cli_args(config[snakemake.rule.name]),
        output_dir = lambda wc: os.path.join(RESULTS_DIR, wc.param_set, wc.sample),
    log:
        lambda wc: os.path.join(RESULTS_DIR, wc.param_set, wc.sample, f"{rule.name}.log")
    shell:
        """
        mkdir -p {params.output_dir} && \
        Rscript {rule.name}.R \
            --counts_dir {input.counts_dir} \
            --out_dir {params.output_dir} \
            {params.main} &> {log}
        """

