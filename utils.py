from pathlib import Path
import os
import re
import shlex
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython.display import FileLink, display
from scipy.optimize import linear_sum_assignment


REPO = Path.cwd()
IGSR_PATH = REPO / "1000Genomes/igsr_samples.tsv"
K_PLOT_LIST = [4, 5, 6, 7]
PRUNED_VCF_TEMPLATE = "1000Genomes/1000G_chr{chr}_pruned.vcf.gz"
RAW_VCF_TEMPLATE = "1000Genomes/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
SPECIAL_CHR = 20
VERBOSE = False  # Set True to print pipeline progress (Running..., log links, etc.)


def plink_file(prefix_path: Path | str, extension: str) -> Path:
    return Path(f"{prefix_path}{extension}")


def experiment_span_label(chr_start: int, chr_end: int) -> str:
    if chr_start == chr_end:
        return f"chr{chr_start}"
    return "ALL"


def experiment_dataset_prefix(
    out_dir: str,
    prefix: str,
    *,
    chr_start: int,
    chr_end: int,
    output_pruned: bool,
) -> Path:
    suffix = ".pruned" if output_pruned else ""
    span_label = experiment_span_label(chr_start, chr_end)
    return Path(out_dir) / f"{prefix}_{span_label}{suffix}"


def show_log_link(message: str, log_path: str | Path, verbose: bool | None = None) -> None:
    if verbose is None:
        verbose = VERBOSE
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.touch(exist_ok=True)
    if verbose:
        print(f"{message}: {log_path}")
        display(FileLink(str(log_path)))


def run_preprocess_if_needed(
    out_dir: str,
    prefix: str,
    sample_info: str = "1000Genomes/igsr_samples.tsv",
    chr_start: int = 1,
    chr_end: int = 22,
    maf: float | None = 0.01,
    ld_window: int = 50,
    ld_step: int = 10,
    ld_r2: float = 0.1,
    ld_pruned: bool = True,
    output_pruned: bool | None = None,
    vcf_template: str = RAW_VCF_TEMPLATE,
    skip_normalize: bool = False,
    force: bool = False,
    log_path: str | Path | None = None,
    verbose: bool | None = None,
) -> Path:
    if output_pruned is None:
        output_pruned = ld_pruned
    v = verbose if verbose is not None else VERBOSE

    dataset_prefix = experiment_dataset_prefix(
        out_dir,
        prefix,
        chr_start=chr_start,
        chr_end=chr_end,
        output_pruned=output_pruned,
    )
    required_paths = [
        plink_file(dataset_prefix, ".bed"),
        plink_file(dataset_prefix, ".bim"),
        plink_file(dataset_prefix, ".fam"),
    ]
    if all(path.exists() for path in required_paths) and not force:
        if v:
            print(f"Found preprocessed data at {dataset_prefix}, skip preprocess.")
        return dataset_prefix

    chromosome_suffix = ".pruned" if output_pruned else ""
    chromosome_ready = all(
        all(
            (Path(out_dir) / f"{prefix}_chr{chr_number}{chromosome_suffix}{extension}").exists()
            for extension in [".bed", ".bim", ".fam"]
        )
        for chr_number in range(chr_start, chr_end + 1)
    )

    cmd = [
        "bash",
        "preprocess.sh",
        "--out-dir",
        out_dir,
        "--prefix",
        prefix,
        "--sample-info",
        sample_info,
        "--chr-start",
        str(chr_start),
        "--chr-end",
        str(chr_end),
        "--ld-window",
        str(ld_window),
        "--ld-step",
        str(ld_step),
        "--ld-r2",
        str(ld_r2),
        "--vcf-template",
        vcf_template,
    ]
    if maf is None:
        cmd.append("--skip-maf")
    else:
        cmd.extend(["--maf", str(maf)])
    cmd.append("--output-pruned" if output_pruned else "--output-unpruned")
    if skip_normalize:
        cmd.append("--skip-normalize")
    if not ld_pruned:
        cmd.append("--skip-ld-prune")

    append_log = False
    if chromosome_ready and not force:
        cmd.append("--skip-chr-process")
        append_log = True
        if v:
            print("Found chromosome-level preprocess outputs, reusing them and rebuilding the dataset only.")
    elif force and v:
            print("Force enabled: rerunning chromosome preprocessing from scratch.")

    log_path = Path(log_path or Path(out_dir) / "preprocess.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_mode = "a" if append_log else "w"
    if v:
        print("Running:", shlex.join(cmd))
        show_log_link("Preprocess log", log_path, verbose=True)
        if append_log:
            print("Appending to the existing preprocess log because chromosome outputs are being reused.")
        print("If preprocessing fails, inspect the log above for details.")
    else:
        show_log_link("Preprocess log", log_path, verbose=False)
    with log_path.open(log_mode) as log_file:
        if append_log:
            log_file.write("\n===== dataset rebuild =====\n")
        subprocess.run(cmd, check=True, stdout=log_file, stderr=subprocess.STDOUT)
    if v:
        print(f"Preprocess completed: {dataset_prefix}")
    return dataset_prefix


def build_structure_ind_file(fam_path: str | Path, ind_file_path: str | Path) -> Path:
    igsr = pd.read_csv(IGSR_PATH, sep="\t")
    phase3 = igsr[
        igsr["Data collections"].str.contains(
            "1000 Genomes phase 3 release",
            case=False,
            na=False,
        )
    ]
    sample_to_pop = dict(zip(phase3["Sample name"], phase3["Population code"]))

    fam = pd.read_csv(
        fam_path,
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PID", "MID", "SEX", "PHENO"],
        dtype=str,
    )
    assert (fam["FID"] == fam["IID"]).all(), "FID != IID: --double-id may not have been used during preprocessing"

    fam["pop"] = fam["IID"].map(sample_to_pop)
    pop_order = sorted(fam["pop"].dropna().unique())
    pop_to_order = {population: i + 1 for i, population in enumerate(pop_order)}
    fam["order"] = fam["pop"].map(pop_to_order)

    ind_file_path = Path(ind_file_path)
    ind_file_path.parent.mkdir(parents=True, exist_ok=True)
    fam[["IID", "pop", "order"]].to_csv(
        ind_file_path,
        sep="\t",
        header=False,
        index=False,
    )
    if VERBOSE:
        print(
            f"Ind file written: {len(fam)} individuals, {len(pop_order)} populations, "
            f"missing_pop={fam['pop'].isna().sum()} -> {ind_file_path}"
        )
    return ind_file_path


def run_method_script(
    script_name: str,
    *,
    out_dir: str,
    input_bed: str | Path,
    prefix: str,
    log_path: str | Path,
    ind_file: str | Path | None = None,
    verbose: bool | None = None,
) -> None:
    if verbose is None:
        verbose = VERBOSE
    env = os.environ.copy()
    env.update(
        {
            "OUT_DIR": str(out_dir),
            "INPUT_BED": str(input_bed),
            "PREFIX": str(prefix),
        }
    )
    if ind_file is not None:
        env["IND_FILE"] = str(ind_file)

    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["bash", script_name]
    if verbose:
        print("Running:", shlex.join(cmd))
        print(f"  OUT_DIR={out_dir}")
        print(f"  INPUT_BED={input_bed}")
        print(f"  PREFIX={prefix}")
        if ind_file is not None:
            print(f"  IND_FILE={ind_file}")
        show_log_link(f"{script_name} log", log_path, verbose=True)
        print("If the command fails, inspect the log above for details.")
    else:
        show_log_link(f"{script_name} log", log_path, verbose=False)
    with log_path.open("w") as log_file:
        subprocess.run(cmd, check=True, env=env, stdout=log_file, stderr=subprocess.STDOUT)
    if verbose:
        print(f"{script_name} finished successfully.")


def prepare_experiment_inputs(
    *,
    experiment_name: str,
    preprocess_out_dir: str,
    preprocess_prefix: str,
    admixture_out_dir: str,
    admixture_prefix: str,
    structure_out_dir: str,
    structure_prefix: str,
    chr_start: int,
    chr_end: int,
    maf: float | None,
    ld_pruned: bool,
    output_pruned: bool,
    vcf_template: str,
    skip_normalize: bool,
    force_preprocess: bool = False,
) -> dict[str, Path | str | int | float | bool | None]:
    dataset_prefix = run_preprocess_if_needed(
        out_dir=preprocess_out_dir,
        prefix=preprocess_prefix,
        chr_start=chr_start,
        chr_end=chr_end,
        maf=maf,
        ld_pruned=ld_pruned,
        output_pruned=output_pruned,
        vcf_template=vcf_template,
        skip_normalize=skip_normalize,
        force=force_preprocess,
    )
    return {
        "experiment_name": experiment_name,
        "dataset_prefix": dataset_prefix,
        "chr_start": chr_start,
        "chr_end": chr_end,
        "maf": maf,
        "ld_pruned": ld_pruned,
        "output_pruned": output_pruned,
        "vcf_template": vcf_template,
        "admixture_out_dir": admixture_out_dir,
        "admixture_prefix": admixture_prefix,
        "structure_out_dir": structure_out_dir,
        "structure_prefix": structure_prefix,
    }


def run_admixture_experiment(result: dict[str, Path | str | int | float | bool | None]) -> dict[str, Path | str | int | float | bool | None]:
    dataset_prefix = Path(result["dataset_prefix"])
    run_method_script(
        "admixture.sh",
        out_dir=result["admixture_out_dir"],
        input_bed=plink_file(dataset_prefix, ".bed"),
        prefix=result["admixture_prefix"],
        log_path=Path(result["admixture_out_dir"]) / "admixture_pipeline.log",
    )
    return result


def run_structure_experiment(result: dict[str, Path | str | int | float | bool | None]) -> dict[str, Path | str | int | float | bool | None]:
    dataset_prefix = Path(result["dataset_prefix"])
    ind_file = build_structure_ind_file(
        plink_file(dataset_prefix, ".fam"),
        Path(result["structure_out_dir"]) / f"{result['structure_prefix']}_ind_file.tsv",
    )
    run_method_script(
        "structure.sh",
        out_dir=result["structure_out_dir"],
        input_bed=plink_file(dataset_prefix, ".bed"),
        prefix=result["structure_prefix"],
        ind_file=ind_file,
        log_path=Path(result["structure_out_dir"]) / "structure_pipeline.log",
    )
    return {
        **result,
        "ind_file": ind_file,
    }


def show_experiment_sanity(
    dataset_prefix: str | Path,
    admixture_dir: str | Path,
    structure_dir: str | Path,
    structure_prefix: str,
    admixture_prefix: str | None = None,
) -> None:
    dataset_prefix = Path(dataset_prefix)
    fam_path = plink_file(dataset_prefix, ".fam")
    bim_path = plink_file(dataset_prefix, ".bim")

    fam_cols = ["FID", "IID", "PID", "MID", "SEX", "PHENO"]
    fam_df = pd.read_csv(fam_path, sep=r"\s+", header=None, names=fam_cols)
    print(f"[{fam_path.name}]")
    print("shape:", fam_df.shape)
    display(fam_df.head(3))

    bim_cols = ["CHR", "SNP", "CM", "BP", "A1", "A2"]
    bim_df = pd.read_csv(bim_path, sep=r"\s+", header=None, names=bim_cols)
    print(f"\n[{bim_path.name}]")
    print("shape:", bim_df.shape)
    display(bim_df.head(3))

    if admixture_prefix is None:
        admixture_prefix = "admixture"

    q_files = sorted(
        Path(admixture_dir).glob(f"{admixture_prefix}_ALL.pruned.*.Q")
    )
    q_summary: list[dict[str, object]] = []
    for qf in q_files:
        match = re.search(r"\.([0-9]+)\.Q$", qf.name)
        if match is None:
            continue
        q_df = pd.read_csv(qf, sep=r"\s+", header=None)
        q_summary.append(
            {
                "file": qf.name,
                "K": int(match.group(1)),
                "n_samples": q_df.shape[0],
                "n_components": q_df.shape[1],
                "matches_fam_rows": q_df.shape[0] == fam_df.shape[0],
            }
        )
    print("\n[ADMIXTURE Q files]")
    q_df_summary = pd.DataFrame(q_summary)
    if q_df_summary.empty:
        raise FileNotFoundError(
            f"No ADMIXTURE Q files found in {admixture_dir} "
            f"matching pattern '{admixture_prefix}_ALL.pruned.*.Q'"
        )
    display(q_df_summary.sort_values("K"))

    meanq_dir = Path(structure_dir) / f"{structure_prefix}_ALL"
    meanq_files = sorted(meanq_dir.glob("fS_run_K.*.meanQ"))
    meanq_summary: list[dict[str, object]] = []
    for meanq_file in meanq_files:
        match = re.search(r"K\.([0-9]+)\.meanQ$", meanq_file.name)
        if match is None:
            continue
        meanq_df = pd.read_csv(meanq_file, sep=r"\s+", header=None)
        meanq_summary.append(
            {
                "file": meanq_file.name,
                "K": int(match.group(1)),
                "n_samples": meanq_df.shape[0],
                "n_components": meanq_df.shape[1],
                "matches_fam_rows": meanq_df.shape[0] == fam_df.shape[0],
            }
        )
    print("\n[fastSTRUCTURE meanQ files]")
    meanq_df_summary = pd.DataFrame(meanq_summary)
    if meanq_df_summary.empty:
        raise FileNotFoundError(
            f"No fastSTRUCTURE meanQ files found in {meanq_dir} "
            "matching pattern 'fS_run_K.*.meanQ'"
        )
    display(meanq_df_summary.sort_values("K"))


def build_population_order(fam_path: str | Path, igsr_path: str | Path = IGSR_PATH) -> pd.DataFrame:
    sampleinfo = pd.read_csv(igsr_path, sep="\t")
    sample_to_pop = dict(zip(sampleinfo["Sample name"], sampleinfo["Population code"]))

    samples = pd.read_csv(fam_path, sep=r"\s+", header=None, usecols=[0]).iloc[:, 0].tolist()
    populations = [sample_to_pop.get(sample, "NA") for sample in samples]

    order_df = pd.DataFrame({"sample": samples, "pop": populations})
    order_df = order_df.sort_values(["pop", "sample"], kind="stable").reset_index(drop=True)
    return order_df



def load_q_matrix_raw(path: str | Path, eps: float = 1e-12) -> np.ndarray:
    """Load Q matrix as raw numpy array (single path only). Clips and normalizes rows."""
    q = pd.read_csv(path, sep=r"\s+", header=None).to_numpy(dtype=float)
    q = np.clip(q, eps, None)
    q /= q.sum(axis=1, keepdims=True)
    return q


def load_iids_from_fam(path: str | Path) -> list[str]:
    """Load IID column (second column) from a PLINK .fam file as a list of strings."""
    fam = pd.read_csv(path, sep=r"\s+", header=None, usecols=[1], names=["IID"], dtype=str)
    return fam["IID"].tolist()


def load_q_matrix(
    q_path: str | Path,
    order_df: pd.DataFrame,
    fam_path: str | Path,
) -> pd.DataFrame:
    q = pd.read_csv(q_path, sep=r"\s+", header=None)
    fam_samples = pd.read_csv(fam_path, sep=r"\s+", header=None, usecols=[0]).iloc[:, 0].tolist()

    q["sample"] = fam_samples
    q = order_df.merge(q, on="sample", how="left")

    k_cols = [col for col in q.columns if isinstance(col, int)]
    return q[k_cols]



import numpy as np
from scipy.optimize import linear_sum_assignment

def align_q_columns(reference_q, target_q, eps=1e-12):
    ref = reference_q.to_numpy(dtype=float)
    tgt = target_q.to_numpy(dtype=float)

    k = ref.shape[1]
    sim = np.zeros((k, k), dtype=float)

    for i in range(k):
        x = ref[:, i]
        sx = np.std(x)

        for j in range(k):
            y = tgt[:, j]
            sy = np.std(y)

            if sx < eps and sy < eps:
                sim[i, j] = 1.0 if np.allclose(x, y, atol=eps, rtol=0.0) else 0.0
            elif sx < eps or sy < eps:
                sim[i, j] = 0.0
            else:
                sim[i, j] = abs(np.corrcoef(x, y)[0, 1])

    _, col_ind = linear_sum_assignment(-sim)

    aligned = target_q.iloc[:, col_ind].copy()
    aligned.columns = reference_q.columns
    return aligned



def get_population_ticks(order_df: pd.DataFrame):
    population_blocks = order_df.groupby("pop", sort=False).size()
    population_centers = population_blocks.cumsum() - (population_blocks.to_numpy() + 1) / 2
    return population_blocks, population_centers



def style_membership_axis(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    ax.tick_params(axis="y", labelsize=11)
    ax.set_xlabel("")
    ax.set_ylabel("")

def get_root_plots_dir() -> Path:
    """Return the root-level plots directory, creating it if needed."""
    plots_dir = Path("plots")
    plots_dir.mkdir(parents=True, exist_ok=True)
    return plots_dir

def build_membership_plot_path(data_path_template: str, title_prefix: str) -> Path:
    parent = Path(data_path_template).parent          # e.g. .../special_chr20_maf001_pruned/admixture
    method_name = parent.name                         # admixture or structure
    experiment_name = parent.parent.name              # e.g. special_chr20_maf001_pruned

    safe_title = re.sub(r"[^A-Za-z0-9]+", "_", title_prefix).strip("_").lower()
    prefix = f"{experiment_name}_{method_name}"

    return get_root_plots_dir() / f"{prefix}_{safe_title}_membership_panels.png"



def plot_membership_panels(
    *,
    fam_path: str | Path,
    data_path_template: str,
    title_prefix: str,
    k_list: list[int] | None = None,
    save_path: str | Path | None = None,
    show_plot: bool = True,
    reference_data_path_template: str | None = None,
) -> Path:
    if k_list is None:
        k_list = K_PLOT_LIST

    order_df = build_population_order(fam_path)
    population_blocks, population_centers = get_population_ticks(order_df)

    fig = plt.figure(figsize=(16, 3.5 * len(k_list)))
    for plot_index, k_value in enumerate(k_list, start=1):
        ax = fig.add_subplot(len(k_list), 1, plot_index)
        data = load_q_matrix(
            data_path_template.format(K=k_value),
            order_df,
            fam_path,
        )
        if reference_data_path_template is not None:
            reference_q = load_q_matrix(
                reference_data_path_template.format(K=k_value),
                order_df,
                fam_path,
            )
            data = align_q_columns(reference_q, data)

        data.plot.bar(stacked=True, ax=ax, width=1)
        style_membership_axis(ax)
        leg = ax.get_legend()
        if leg is not None:
            leg.set_loc("center")

        if plot_index < len(k_list):
            ax.set_xticks([])
        else:
            ax.set_xticks(population_centers.to_list())
            ax.set_xticklabels(population_blocks.index.to_list(), rotation=90, fontsize=12)
        ax.set_title(f"{title_prefix}, K={k_value}")

    save_path = Path(save_path) if save_path is not None else build_membership_plot_path(data_path_template, title_prefix)

    fig.tight_layout()
    fig.savefig(save_path, dpi=330, bbox_inches="tight")
    pdf_path = save_path.with_suffix(".pdf")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Membership plot saved to: {save_path}")
    print(f"Membership plot saved to: {pdf_path}")

    if show_plot:
        plt.show()
    plt.close(fig)
    return save_path


def plot_combined_membership_panels(
    *,
    fam_path: str | Path,
    panels_config: list[dict],
    k_list: list[int] | None = None,
    save_path: str | Path | None = None,
    show_plot: bool = True,
) -> Path:
    """Plot multiple methods in one figure: each row = one K, each column = one method."""
    if k_list is None:
        k_list = K_PLOT_LIST
    n_rows = len(k_list)
    n_cols = len(panels_config)
    order_df = build_population_order(fam_path)
    population_blocks, population_centers = get_population_ticks(order_df)
    fig = plt.figure(figsize=(16 * n_cols, 3.5 * n_rows))
    for row_i, k_value in enumerate(k_list):
        for col_j, cfg in enumerate(panels_config):
            ax = fig.add_subplot(n_rows, n_cols, row_i * n_cols + col_j + 1)
            data_path_template = cfg["data_path_template"]
            title_prefix = cfg["title_prefix"]
            ref_tpl = cfg.get("reference_data_path_template")
            data = load_q_matrix(
                data_path_template.format(K=k_value),
                order_df,
                fam_path,
            )
            if ref_tpl is not None:
                reference_q = load_q_matrix(ref_tpl.format(K=k_value), order_df, fam_path)
                data = align_q_columns(reference_q, data)
            data.plot.bar(stacked=True, ax=ax, width=1)
            style_membership_axis(ax)
            leg = ax.get_legend()
            if leg is not None:
                leg.set_loc("center")
            if row_i < n_rows - 1:
                ax.set_xticks([])
            else:
                ax.set_xticks(population_centers.to_list())
                ax.set_xticklabels(population_blocks.index.to_list(), rotation=90, fontsize=12)
            ax.set_title(f"{title_prefix}, K={k_value}")
    template0 = panels_config[0]["data_path_template"]
    parent0 = Path(template0).parent
    method0 = parent0.name
    experiment0 = parent0.parent.name
    prefix = f"{experiment0}_{method0}"

    out_path = Path(save_path) if save_path is not None else get_root_plots_dir() / f"{prefix}_combined_membership_panels.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=330, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Combined membership plot saved to: {out_path}")
    if show_plot:
        plt.show()
    plt.close(fig)
    return out_path