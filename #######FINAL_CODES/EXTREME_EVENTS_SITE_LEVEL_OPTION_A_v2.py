"""
Extreme events analysis (Option A) - v2.

Enhancements over v1:
1) Keeps profile-level analysis (SITE x HAB) as the primary unit.
2) Computes and visualizes peak, duration, and frequency separately for
   Inner North (ITA*), Inner Central (TIM*), and Inner South (UCR).
3) Replaces Fig4 example layout with five two-panel figures where each panel
   overlays DLI, SST, and Chl-a anomalies using three independent y-axes.
4) Preserves CV-window HAB/sector figure and generates an automatic
   manuscript-style report in English.
"""

import glob
import os
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

import EXTREME_EVENTS_SITE_LEVEL_OPTION_A as core


# --------------------------------------------
# v2 configuration
# --------------------------------------------
OUTPUT_DIR_V2 = os.path.join(
    core.BASE_DIR, "#######FINAL_RESULTS", "Extreme_Events_Site_Level_v2"
)
os.makedirs(OUTPUT_DIR_V2, exist_ok=True)

INNER_SECTOR_ORDER = ["inner_north", "inner_central", "inner_south"]
INNER_SECTOR_LABELS = {
    "inner_north": "Inner North (ITA)",
    "inner_central": "Inner Central (TIM)",
    "inner_south": "Inner South (UCR)",
}
INNER_SECTOR_COLORS = {
    "inner_north": "#E69F00",
    "inner_central": "#009E73",
    "inner_south": "#0072B2",
}

SECTOR_ORDER = ["inner_north", "inner_central", "inner_south", "outer"]
SECTOR_LABELS = {
    "inner_north": "Inner North",
    "inner_central": "Inner Central",
    "inner_south": "Inner South",
    "outer": "Outer",
}
SECTOR_LABELS_SHORT = {
    "inner_north": "Inner N",
    "inner_central": "Inner C",
    "inner_south": "Inner S",
    "outer": "Outer",
}
HAB_ORDER = ["TP", "PA", "RR"]


def compute_inner_sector_stats(metrics_df):
    """Compute inner-sector statistics for each variable and metric."""
    kw_rows = []
    pair_rows = []

    for var in ["SST", "DLI", "Chl-a"]:
        base = metrics_df[
            (metrics_df["VARIABLE"] == var)
            & (metrics_df["ARC"] == "inner")
            & (metrics_df["INNER_SECTOR"].isin(INNER_SECTOR_ORDER))
        ].copy()

        for metric in [
            "n_events_per_year",
            "mean_duration_d",
            "mean_peak_z",
            "pct_extreme_days",
        ]:
            sub = base.dropna(subset=[metric]).copy()
            if metric in ["mean_duration_d", "mean_peak_z"]:
                sub = sub[sub["n_events"] > 0].copy()

            sec_vals = {}
            for sec in INNER_SECTOR_ORDER:
                vals = (
                    sub[sub["INNER_SECTOR"] == sec][metric]
                    .dropna()
                    .values.astype(float)
                )
                sec_vals[sec] = vals

            kr_groups = [
                sec_vals[s] for s in INNER_SECTOR_ORDER if len(sec_vals[s]) >= 2
            ]
            if len(kr_groups) >= 2:
                try:
                    h_stat, p_kw = stats.kruskal(*kr_groups)
                except Exception:
                    h_stat, p_kw = np.nan, np.nan
            else:
                h_stat, p_kw = np.nan, np.nan

            kw_rows.append(
                {
                    "Variable": var,
                    "Metric": metric,
                    "Mean_inner_north": np.nanmean(sec_vals["inner_north"])
                    if len(sec_vals["inner_north"])
                    else np.nan,
                    "Mean_inner_central": np.nanmean(sec_vals["inner_central"])
                    if len(sec_vals["inner_central"])
                    else np.nan,
                    "Mean_inner_south": np.nanmean(sec_vals["inner_south"])
                    if len(sec_vals["inner_south"])
                    else np.nan,
                    "N_inner_north": len(sec_vals["inner_north"]),
                    "N_inner_central": len(sec_vals["inner_central"]),
                    "N_inner_south": len(sec_vals["inner_south"]),
                    "Kruskal_H": h_stat,
                    "Kruskal_p": p_kw,
                }
            )

            for sec_a, sec_b in combinations(INNER_SECTOR_ORDER, 2):
                vals_a = sec_vals[sec_a]
                vals_b = sec_vals[sec_b]
                perm = core.permutation_test(vals_a, vals_b)
                p_mw = core.mann_whitney_p(vals_a, vals_b)
                pair_rows.append(
                    {
                        "Variable": var,
                        "Metric": metric,
                        "Sector_A": sec_a,
                        "Sector_B": sec_b,
                        "Mean_A": np.nanmean(vals_a) if len(vals_a) else np.nan,
                        "Mean_B": np.nanmean(vals_b) if len(vals_b) else np.nan,
                        "Obs_Difference_A_minus_B": perm["obs_diff"],
                        "P_Value_Perm": perm["p_value"],
                        "P_Value_MWU": p_mw,
                        "Effect_Size_CohenD": perm["effect_size"],
                        "N_A": len(vals_a),
                        "N_B": len(vals_b),
                    }
                )

    return pd.DataFrame(kw_rows), pd.DataFrame(pair_rows)


def summarize_event_metrics_by_sector_hab(metrics_df):
    """Long-format summary for event metrics by variable, sector, and habitat."""
    rows = []
    metrics = [
        "n_events_per_year",
        "mean_duration_d",
        "mean_peak_z",
        "pct_extreme_days",
    ]

    for var in ["SST", "DLI", "Chl-a"]:
        sub_var = metrics_df[
            (metrics_df["VARIABLE"] == var)
            & (metrics_df["INNER_SECTOR"].isin(SECTOR_ORDER))
        ].copy()
        for metric in metrics:
            sub = sub_var.copy()
            if metric in ["mean_duration_d", "mean_peak_z"]:
                sub = sub[sub["n_events"] > 0]

            grouped = (
                sub.groupby(["ARC", "INNER_SECTOR", "HAB"], dropna=False)[metric]
                .agg(["count", "mean", "median", "std", "min", "max"])
                .reset_index()
            )
            grouped["VARIABLE"] = var
            grouped["METRIC"] = metric
            grouped = grouped.rename(
                columns={
                    "count": "N",
                    "mean": "MEAN",
                    "median": "MEDIAN",
                    "std": "SD",
                    "min": "MIN",
                    "max": "MAX",
                }
            )
            rows.append(grouped)

    if not rows:
        return pd.DataFrame(
            columns=[
                "VARIABLE",
                "METRIC",
                "ARC",
                "INNER_SECTOR",
                "HAB",
                "N",
                "MEAN",
                "MEDIAN",
                "SD",
                "MIN",
                "MAX",
            ]
        )

    out = pd.concat(rows, ignore_index=True)
    cols = [
        "VARIABLE",
        "METRIC",
        "ARC",
        "INNER_SECTOR",
        "HAB",
        "N",
        "MEAN",
        "MEDIAN",
        "SD",
        "MIN",
        "MAX",
    ]
    return out[cols]


def create_inner_sector_metric_figure(
    metrics_df,
    value_col,
    ylabel,
    title,
    output_path,
    events_only=False,
):
    """Create 3-panel figure (SST, DLI, Chl-a) split by inner sectors."""
    var_order = ["SST", "DLI", "Chl-a"]
    var_colors = {"SST": core.COLOR_SST, "DLI": core.COLOR_DLI, "Chl-a": core.COLOR_CHL}

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    for ax, var in zip(axes, var_order):
        sub = metrics_df[
            (metrics_df["VARIABLE"] == var)
            & (metrics_df["ARC"] == "inner")
            & (metrics_df["INNER_SECTOR"].isin(INNER_SECTOR_ORDER))
        ].copy()
        if events_only:
            sub = sub[sub["n_events"] > 0].copy()

        sub = sub.dropna(subset=[value_col])

        plot_data = []
        positions = []
        colors = []
        for i, sec in enumerate(INNER_SECTOR_ORDER, start=1):
            vals = sub[sub["INNER_SECTOR"] == sec][value_col].dropna().values
            if len(vals) == 0:
                continue
            plot_data.append(vals)
            positions.append(i)
            colors.append(INNER_SECTOR_COLORS[sec])

        if len(plot_data) == 0:
            ax.text(
                0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center"
            )
            ax.set_title(var)
            continue

        bp = ax.boxplot(
            plot_data,
            positions=positions,
            widths=0.6,
            patch_artist=True,
            showfliers=False,
        )
        for b, c in zip(bp["boxes"], colors):
            b.set_facecolor(c)
            b.set_alpha(0.75)

        for pos, vals, c in zip(positions, plot_data, colors):
            xj = np.random.normal(pos, 0.05, len(vals))
            ax.scatter(
                xj,
                vals,
                c=c,
                s=24,
                alpha=0.55,
                edgecolors="white",
                linewidth=0.35,
                zorder=3,
            )

        # Kruskal-Wallis annotation
        kr_groups = [arr for arr in plot_data if len(arr) >= 2]
        if len(kr_groups) >= 2:
            try:
                _, p_kw = stats.kruskal(*kr_groups)
            except Exception:
                p_kw = np.nan
        else:
            p_kw = np.nan

        ax.text(
            0.97,
            0.95,
            f"Kruskal p={core.format_p_value(p_kw)}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.55),
        )

        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(
            [INNER_SECTOR_LABELS[s] for s in INNER_SECTOR_ORDER], rotation=15
        )
        ax.set_title(var, color=var_colors[var], fontweight="bold")
        ax.set_ylabel(ylabel)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def create_sector_hab_metric_figure(
    metrics_df,
    value_col,
    ylabel,
    title,
    output_path,
    events_only=False,
):
    """
    Create 3-panel figure (SST, DLI, Chl-a) with sector + HAB structure.

    X-axis groups: Inner North, Inner Central, Inner South, Outer.
    Within each group: HAB-specific boxplots (TP, PA, RR) where data are available.
    """
    var_order = ["SST", "DLI", "Chl-a"]
    var_colors = {"SST": core.COLOR_SST, "DLI": core.COLOR_DLI, "Chl-a": core.COLOR_CHL}
    hab_offsets = {"TP": -0.25, "PA": 0.00, "RR": 0.25}

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), sharey=False)
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    for ax, var in zip(axes, var_order):
        sub = metrics_df[
            (metrics_df["VARIABLE"] == var)
            & (metrics_df["INNER_SECTOR"].isin(SECTOR_ORDER))
        ].copy()
        if events_only:
            sub = sub[sub["n_events"] > 0].copy()
        sub = sub.dropna(subset=[value_col])

        if len(sub) == 0:
            ax.text(
                0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center"
            )
            ax.set_title(var)
            continue

        box_data = []
        box_pos = []
        box_hab = []

        for i, sec in enumerate(SECTOR_ORDER):
            sec_sub = sub[sub["INNER_SECTOR"] == sec]
            for hab in HAB_ORDER:
                vals = sec_sub[sec_sub["HAB"] == hab][value_col].dropna().values
                if len(vals) == 0:
                    continue
                box_data.append(vals)
                box_pos.append(i + hab_offsets[hab])
                box_hab.append(hab)

        if len(box_data) == 0:
            ax.text(
                0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center"
            )
            ax.set_title(var)
            continue

        bp = ax.boxplot(
            box_data,
            positions=box_pos,
            widths=0.18,
            patch_artist=True,
            showfliers=False,
        )
        for b, hab in zip(bp["boxes"], box_hab):
            b.set_facecolor(core.HAB_COLORS.get(hab, "gray"))
            b.set_alpha(0.75)

        for pos, vals, hab in zip(box_pos, box_data, box_hab):
            xj = np.random.normal(pos, 0.025, len(vals))
            ax.scatter(
                xj,
                vals,
                c=core.HAB_COLORS.get(hab, "gray"),
                s=18,
                alpha=0.5,
                edgecolors="white",
                linewidth=0.35,
                zorder=3,
            )

        # Kruskal across all sector x HAB cells with n>=2
        grouped = []
        for sec in SECTOR_ORDER:
            for hab in HAB_ORDER:
                vals = (
                    sub[(sub["INNER_SECTOR"] == sec) & (sub["HAB"] == hab)][value_col]
                    .dropna()
                    .values
                )
                if len(vals) >= 2:
                    grouped.append(vals)
        if len(grouped) >= 2:
            try:
                _, p_kw = stats.kruskal(*grouped)
            except Exception:
                p_kw = np.nan
        else:
            p_kw = np.nan

        ax.text(
            0.98,
            0.95,
            f"Kruskal p={core.format_p_value(p_kw)}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.55),
        )

        ax.set_xticks(np.arange(len(SECTOR_ORDER)))
        ax.set_xticklabels([SECTOR_LABELS[s] for s in SECTOR_ORDER], rotation=15)
        ax.set_title(var, color=var_colors[var], fontweight="bold")
        ax.set_ylabel(ylabel)

    legend = [
        Patch(facecolor=core.HAB_COLORS[h], edgecolor="black", alpha=0.75, label=h)
        for h in HAB_ORDER
    ]
    fig.legend(
        handles=legend, title="HAB", loc="upper right", bbox_to_anchor=(0.99, 0.99)
    )
    plt.tight_layout(rect=(0, 0, 0.97, 0.96))
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def _plot_sector_hab_boxpanel(
    ax,
    sub,
    variable,
    value_col,
    ylabel,
    title,
    add_ylabel=False,
    add_xlabel=False,
):
    """Plot a single sector x HAB boxpanel with robust annotation placement."""
    fs_no_data = 14
    fs_anno = 13
    fs_xtick = 14
    fs_ytick = 14
    fs_title = 18
    fs_axis_label = 16

    box_data = []
    box_pos = []
    box_hab = []
    hab_offsets = {"TP": -0.22, "PA": 0.0, "RR": 0.22}

    if variable == "DLI":
        for i, sec in enumerate(SECTOR_ORDER):
            sec_sub = sub[sub["INNER_SECTOR"] == sec]
            for hab in HAB_ORDER:
                vals = sec_sub[sec_sub["HAB"] == hab][value_col].dropna().values
                if len(vals) == 0:
                    continue
                box_data.append(vals)
                box_pos.append(i + hab_offsets[hab])
                box_hab.append(hab)
        box_width = 0.24
    else:
        # SST and Chl-a are pooled across HAB within SITE.
        for i, sec in enumerate(SECTOR_ORDER):
            sec_sub = sub[sub["INNER_SECTOR"] == sec]
            if len(sec_sub) == 0:
                continue
            vals = (
                sec_sub.groupby("SITE_ID", dropna=False)[value_col]
                .mean()
                .dropna()
                .values
            )
            if len(vals) == 0:
                continue
            box_data.append(vals)
            box_pos.append(i)
            box_hab.append("POOLED")
        box_width = 0.48

    if len(box_data) == 0:
        ax.text(
            0.5,
            0.5,
            "No data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=fs_no_data,
        )
        ax.set_title(title)
        return

    bp = ax.boxplot(
        box_data,
        positions=box_pos,
        widths=box_width,
        patch_artist=True,
        showfliers=False,
    )
    for b, hab in zip(bp["boxes"], box_hab):
        if hab == "POOLED":
            b.set_facecolor("#A9A9A9")
        else:
            b.set_facecolor(core.HAB_COLORS.get(hab, "gray"))
        b.set_alpha(0.78)

    for pos, vals, hab in zip(box_pos, box_data, box_hab):
        xj = np.random.normal(pos, 0.025, len(vals))
        c = "#808080" if hab == "POOLED" else core.HAB_COLORS.get(hab, "gray")
        ax.scatter(
            xj,
            vals,
            c=c,
            s=15,
            alpha=0.45,
            edgecolors="white",
            linewidth=0.25,
            zorder=3,
        )

    # Add vertical separators between sectors
    for sep in [0.5, 1.5, 2.5]:
        ax.axvline(sep, color="#DDDDDD", linewidth=0.8, zorder=0)

    # Kruskal across groups with n>=2
    grouped = []
    if variable == "DLI":
        for sec in SECTOR_ORDER:
            for hab in HAB_ORDER:
                vals = (
                    sub[(sub["INNER_SECTOR"] == sec) & (sub["HAB"] == hab)][value_col]
                    .dropna()
                    .values
                )
                if len(vals) >= 2:
                    grouped.append(vals)
    else:
        for sec in SECTOR_ORDER:
            vals = (
                sub[sub["INNER_SECTOR"] == sec]
                .groupby("SITE_ID", dropna=False)[value_col]
                .mean()
                .dropna()
                .values
            )
            if len(vals) >= 2:
                grouped.append(vals)
    if len(grouped) >= 2:
        try:
            _, p_kw = stats.kruskal(*grouped)
        except Exception:
            p_kw = np.nan
    else:
        p_kw = np.nan

    # Dynamic y-limits with extra headroom to avoid text overlap
    all_vals = np.concatenate([arr for arr in box_data if len(arr) > 0])
    ymin = np.nanmin(all_vals)
    ymax = np.nanmax(all_vals)
    yr = ymax - ymin
    if yr <= 0:
        yr = max(1.0, abs(ymax) * 0.2 + 0.2)
    pad_bottom = 0.10 * yr
    pad_top = 0.30 * yr
    if ymin >= 0:
        ax.set_ylim(max(0, ymin - pad_bottom), ymax + pad_top)
    else:
        ax.set_ylim(ymin - pad_bottom, ymax + pad_top)

    # Put significance at top-left, clear of legend and boxes
    ax.text(
        0.02,
        0.98,
        f"Kruskal p={core.format_p_value(p_kw)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs_anno,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="#BBBBBB"),
    )

    ax.set_xticks(np.arange(len(SECTOR_ORDER)))
    ax.set_xticklabels(
        [SECTOR_LABELS_SHORT[s] for s in SECTOR_ORDER], rotation=12, fontsize=fs_xtick
    )
    ax.tick_params(axis="y", labelsize=fs_ytick)
    ax.set_title(title, fontweight="bold", fontsize=fs_title)
    if add_ylabel:
        ax.set_ylabel(ylabel, fontsize=fs_axis_label, labelpad=16)
    if add_xlabel:
        ax.set_xlabel("Sector", fontsize=fs_axis_label)


def create_composite_event_metrics_figure(metrics_df, output_path):
    """
    Build a 3x3 composite:
    - Rows: Frequency, Duration, Peak intensity
    - Columns: SST, DLI, Chl-a
    """
    row_specs = [
        ("n_events_per_year", "Event frequency (events/year)", False),
        ("mean_duration_d", "Event duration (days)", True),
        ("mean_peak_z", "Peak intensity (z)", True),
    ]
    var_order = ["SST", "DLI", "Chl-a"]
    var_colors = {"SST": core.COLOR_SST, "DLI": core.COLOR_DLI, "Chl-a": core.COLOR_CHL}

    fig, axes = plt.subplots(3, 3, figsize=(19, 15.8), sharex=False)
    fig.suptitle(
        "Extreme Event Metrics by Sector, Habitat, and Arc (SITE x HAB)",
        fontsize=24,
        fontweight="bold",
        y=0.995,
    )

    for r, (metric, ylabel, events_only) in enumerate(row_specs):
        for c, var in enumerate(var_order):
            ax = axes[r, c]
            sub = metrics_df[
                (metrics_df["VARIABLE"] == var)
                & (metrics_df["INNER_SECTOR"].isin(SECTOR_ORDER))
            ].copy()
            if events_only:
                sub = sub[sub["n_events"] > 0]
            sub = sub.dropna(subset=[metric])

            _plot_sector_hab_boxpanel(
                ax=ax,
                sub=sub,
                variable=var,
                value_col=metric,
                ylabel=ylabel,
                title=var,
                add_ylabel=(c == 0),
                add_xlabel=(r == 2),
            )
            # Color column titles by variable
            ax.title.set_color(var_colors[var])

    hab_legend = [
        Patch(
            facecolor=core.HAB_COLORS[h],
            edgecolor="black",
            alpha=0.78,
            label=f"HAB {h}",
        )
        for h in HAB_ORDER
    ]
    hab_legend.append(
        Patch(facecolor="#A9A9A9", edgecolor="black", alpha=0.78, label="HAB pooled")
    )
    fig.legend(
        handles=hab_legend,
        title="Habitats",
        loc="upper right",
        bbox_to_anchor=(0.995, 0.995),
        frameon=True,
        fontsize=14,
        title_fontsize=15,
    )

    plt.tight_layout(rect=(0.06, 0, 0.96, 0.975))
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def _plot_three_variable_overlay(ax, profile_id, anomaly_dict, event_dict):
    """Plot DLI, SST and Chl-a anomalies with three independent y-axes."""
    dli = anomaly_dict["DLI"].get(profile_id)
    sst = anomaly_dict["SST"].get(profile_id)
    chl = anomaly_dict["Chl-a"].get(profile_id)

    if dli is None and sst is None and chl is None:
        ax.text(
            0.5,
            0.5,
            f"No data for {profile_id}",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )
        return

    ax2 = ax.twinx()
    ax3 = ax.twinx()
    ax3.spines["right"].set_position(("axes", 1.12))
    ax3.spines["right"].set_visible(True)

    if dli is not None:
        ax.plot(dli.index, dli.values, color=core.COLOR_DLI, linewidth=0.8, alpha=0.9)
        ax.axhline(
            core.THRESHOLD_Z,
            color=core.COLOR_DLI,
            linestyle="--",
            linewidth=0.9,
            alpha=0.5,
        )
        ev = event_dict["DLI"].get(profile_id, pd.DataFrame())
        if ev is not None and len(ev) > 0:
            for _, r in ev.iterrows():
                ax.axvspan(
                    r["start_date"], r["end_date"], color=core.COLOR_DLI, alpha=0.06
                )

    if sst is not None:
        ax2.plot(sst.index, sst.values, color=core.COLOR_SST, linewidth=0.8, alpha=0.85)
        ax2.axhline(
            core.THRESHOLD_Z,
            color=core.COLOR_SST,
            linestyle="--",
            linewidth=0.9,
            alpha=0.5,
        )
        ev = event_dict["SST"].get(profile_id, pd.DataFrame())
        if ev is not None and len(ev) > 0:
            for _, r in ev.iterrows():
                ax.axvspan(
                    r["start_date"], r["end_date"], color=core.COLOR_SST, alpha=0.04
                )

    if chl is not None:
        ax3.plot(chl.index, chl.values, color=core.COLOR_CHL, linewidth=0.8, alpha=0.85)
        ax3.axhline(
            core.THRESHOLD_Z,
            color=core.COLOR_CHL,
            linestyle="--",
            linewidth=0.9,
            alpha=0.5,
        )
        ev = event_dict["Chl-a"].get(profile_id, pd.DataFrame())
        if ev is not None and len(ev) > 0:
            for _, r in ev.iterrows():
                ax.axvspan(
                    r["start_date"], r["end_date"], color=core.COLOR_CHL, alpha=0.04
                )

    ax.set_ylabel("DLI anomaly (z)", color=core.COLOR_DLI)
    ax2.set_ylabel("SST anomaly (z)", color=core.COLOR_SST)
    ax3.set_ylabel("Chl-a anomaly (z)", color=core.COLOR_CHL)
    ax.tick_params(axis="y", colors=core.COLOR_DLI)
    ax2.tick_params(axis="y", colors=core.COLOR_SST)
    ax3.tick_params(axis="y", colors=core.COLOR_CHL)
    ax.axhline(0.0, color="gray", linewidth=0.8, alpha=0.45)
    ax.set_title(profile_id, loc="left", fontweight="bold")

    legend_items = [
        Line2D([0], [0], color=core.COLOR_DLI, lw=1.2, label="DLI"),
        Line2D([0], [0], color=core.COLOR_SST, lw=1.2, label="SST"),
        Line2D([0], [0], color=core.COLOR_CHL, lw=1.2, label="Chl-a"),
        Line2D(
            [0],
            [0],
            color="black",
            lw=1.0,
            linestyle="--",
            label=f"Threshold z={core.THRESHOLD_Z}",
        ),
    ]
    ax.legend(handles=legend_items, loc="upper right", fontsize=8, frameon=True)


def create_fig4_overlay_pair(
    anomaly_dict, event_dict, profile_top, profile_bottom, title, output_path
):
    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    _plot_three_variable_overlay(axes[0], profile_top, anomaly_dict, event_dict)
    _plot_three_variable_overlay(axes[1], profile_bottom, anomaly_dict, event_dict)
    axes[1].set_xlabel("Date")
    fig.suptitle(title, fontsize=13, fontweight="bold")
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def load_cv_long_table_multivar():
    """Load CV_2 / CV_30 / CV_ALL tables for DLI, SST and Chl-a in long format."""
    specs = [
        (
            core.CV2_PATH,
            "CV_2",
            {"DLI": "dli_cv_2", "SST": "sst_cv_2", "Chl-a": "chl_cv_2"},
        ),
        (
            core.CV30_PATH,
            "CV_30",
            {"DLI": "dli_cv_30", "SST": "sst_cv_30", "Chl-a": "chl_cv_30"},
        ),
        (
            core.CVALL_PATH,
            "CV_ALL",
            {"DLI": "cv_DLI_local", "SST": "sst_cv_all", "Chl-a": "chl_cv_all"},
        ),
    ]

    out = []
    for path, window, var_map in specs:
        df = pd.read_csv(
            path, sep=";", decimal=",", encoding="utf-8-sig", low_memory=False
        )
        df.columns = [c.strip() for c in df.columns]

        base_cols = ["Reef_name", "Site_name", "HAB", "Arch"]
        missing_base = [c for c in base_cols if c not in df.columns]
        if missing_base:
            raise ValueError(
                f"Missing required columns in {os.path.basename(path)}: {missing_base}"
            )

        for var, col in var_map.items():
            if col not in df.columns:
                raise ValueError(f"Missing column {col} in {os.path.basename(path)}")

            d = df[base_cols + [col]].copy()
            d = d.rename(
                columns={
                    "Reef_name": "REEF_NAME",
                    "Site_name": "SITE_ID",
                    "Arch": "ARC",
                    col: "CV_VALUE",
                }
            )
            d["ARC"] = d["ARC"].astype(str).str.strip().str.lower()
            d = d[d["ARC"].isin(["inner", "outer"])].copy()
            d["CV_VALUE"] = pd.to_numeric(d["CV_VALUE"], errors="coerce")
            d = d.dropna(subset=["CV_VALUE"])
            d["VARIABLE"] = var
            d["WINDOW"] = window
            d["INNER_SECTOR"] = d.apply(
                lambda r: core.classify_inner_sector(
                    r["SITE_ID"], r["ARC"], r["REEF_NAME"]
                ),
                axis=1,
            )
            out.append(d)

    return pd.concat(out, ignore_index=True)


def create_cv_hab_sector_figure_multivar(cv_long, output_path):
    """Create Fig5 as a 3x3 composite (rows=CV windows, cols=variables)."""
    print("\n  Creating multi-variable CV figure (Fig5 3x3)...")

    # Font sizes tuned for publication readability (applies to this figure only).
    fs_suptitle = 24
    fs_col_title = 19
    fs_row_ylabel = 16
    fs_xlabel = 16
    fs_xtick = 14
    fs_ytick = 14
    fs_anno = 13
    fs_legend = 14
    fs_legend_title = 15

    var_order = ["SST", "DLI", "Chl-a"]
    window_order = ["CV_2", "CV_30", "CV_ALL"]
    var_colors = {"SST": core.COLOR_SST, "DLI": core.COLOR_DLI, "Chl-a": core.COLOR_CHL}
    hab_offsets = {"TP": -0.25, "PA": 0.0, "RR": 0.25}

    fig, axes = plt.subplots(3, 3, figsize=(19, 15), sharex=False, sharey=False)
    fig.suptitle(
        "Coefficient of Variation by Window, Sector, Habitat, and Variable",
        fontsize=fs_suptitle,
        fontweight="bold",
        y=0.995,
    )

    for r, window in enumerate(window_order):
        for c, var in enumerate(var_order):
            ax = axes[r, c]
            sub = cv_long[
                (cv_long["WINDOW"] == window)
                & (cv_long["VARIABLE"] == var)
                & (cv_long["INNER_SECTOR"].isin(SECTOR_ORDER))
            ].copy()

            if len(sub) == 0:
                ax.text(
                    0.5,
                    0.5,
                    "No data",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=fs_xtick,
                )
                ax.set_title(f"{var} | {window}")
                continue

            box_data = []
            box_pos = []
            box_hab = []
            if var == "DLI":
                for i, sec in enumerate(SECTOR_ORDER):
                    sec_sub = sub[sub["INNER_SECTOR"] == sec]
                    for hab in HAB_ORDER:
                        vals = (
                            sec_sub[sec_sub["HAB"] == hab]["CV_VALUE"].dropna().values
                        )
                        if len(vals) == 0:
                            continue
                        box_data.append(vals)
                        box_pos.append(i + hab_offsets[hab])
                        box_hab.append(hab)
                box_width = 0.24
            else:
                # HAB pooled (single gray column) for SST and Chl-a.
                for i, sec in enumerate(SECTOR_ORDER):
                    sec_sub = sub[sub["INNER_SECTOR"] == sec]
                    vals = (
                        sec_sub.groupby("SITE_ID", dropna=False)["CV_VALUE"]
                        .mean()
                        .dropna()
                        .values
                    )
                    if len(vals) == 0:
                        continue
                    box_data.append(vals)
                    box_pos.append(i)
                    box_hab.append("POOLED")
                box_width = 0.56

            if len(box_data) == 0:
                ax.text(
                    0.5,
                    0.5,
                    "No data",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=fs_xtick,
                )
                ax.set_title(f"{var} | {window}")
                continue

            bp = ax.boxplot(
                box_data,
                positions=box_pos,
                widths=box_width,
                patch_artist=True,
                showfliers=False,
            )
            for b, hab in zip(bp["boxes"], box_hab):
                if hab == "POOLED":
                    b.set_facecolor("#A9A9A9")
                else:
                    b.set_facecolor(core.HAB_COLORS.get(hab, "gray"))
                b.set_alpha(0.78)

            for pos, vals, hab in zip(box_pos, box_data, box_hab):
                xj = np.random.normal(pos, 0.025, len(vals))
                c = "#808080" if hab == "POOLED" else core.HAB_COLORS.get(hab, "gray")
                ax.scatter(
                    xj,
                    vals,
                    c=c,
                    s=10,
                    alpha=0.45,
                    edgecolors="none",
                    zorder=3,
                )

            # vertical separators between sectors
            for sep in [0.5, 1.5, 2.5]:
                ax.axvline(sep, color="#DDDDDD", linewidth=0.8, zorder=0)

            grouped = []
            if var == "DLI":
                for sec in SECTOR_ORDER:
                    for hab in HAB_ORDER:
                        vals = (
                            sub[(sub["INNER_SECTOR"] == sec) & (sub["HAB"] == hab)][
                                "CV_VALUE"
                            ]
                            .dropna()
                            .values
                        )
                        if len(vals) >= 2:
                            grouped.append(vals)
            else:
                for sec in SECTOR_ORDER:
                    vals = (
                        sub[sub["INNER_SECTOR"] == sec]
                        .groupby("SITE_ID", dropna=False)["CV_VALUE"]
                        .mean()
                        .dropna()
                        .values
                    )
                    if len(vals) >= 2:
                        grouped.append(vals)
            if len(grouped) >= 2:
                try:
                    _, p_kw = stats.kruskal(*grouped)
                except Exception:
                    p_kw = np.nan
            else:
                p_kw = np.nan

            # y headroom and annotation position
            all_vals = np.concatenate([arr for arr in box_data if len(arr) > 0])
            ymin = np.nanmin(all_vals)
            ymax = np.nanmax(all_vals)
            yr = ymax - ymin
            if yr <= 0:
                yr = max(1.0, abs(ymax) * 0.2 + 0.2)
            ax.set_ylim(max(0, ymin - 0.10 * yr), ymax + 0.30 * yr)

            ax.text(
                0.02,
                0.98,
                f"Kruskal p={core.format_p_value(p_kw)}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=fs_anno,
                bbox=dict(
                    boxstyle="round", facecolor="white", alpha=0.8, edgecolor="#BBBBBB"
                ),
            )

            ax.set_xticks(np.arange(len(SECTOR_ORDER)))
            ax.set_xticklabels(
                [SECTOR_LABELS_SHORT[s] for s in SECTOR_ORDER],
                rotation=12,
                fontsize=fs_xtick,
            )
            ax.tick_params(axis="x", labelsize=fs_xtick)
            ax.tick_params(axis="y", labelsize=fs_ytick)

            if r == 0:
                ax.set_title(
                    var, color=var_colors[var], fontweight="bold", fontsize=fs_col_title
                )
            if c == 0:
                ax.set_ylabel(
                    f"{window}\nCoefficient of variation (%)",
                    fontsize=fs_row_ylabel,
                    labelpad=20,
                )
            if r == 2:
                ax.set_xlabel("Sector", fontsize=fs_xlabel)

    hab_legend = [
        Patch(
            facecolor=core.HAB_COLORS[h],
            edgecolor="black",
            alpha=0.78,
            label=f"HAB {h}",
        )
        for h in HAB_ORDER
    ]
    hab_legend.append(
        Patch(facecolor="#A9A9A9", edgecolor="black", alpha=0.78, label="HAB pooled")
    )
    fig.legend(
        handles=hab_legend,
        title="Habitats",
        loc="upper right",
        bbox_to_anchor=(0.995, 0.995),
        fontsize=fs_legend,
        title_fontsize=fs_legend_title,
    )
    plt.tight_layout(rect=(0.08, 0, 0.96, 0.975))
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def generate_report_v2(
    metrics_df, arc_stat_df, inner_kw_df, sector_hab_df, output_path
):
    """Generate manuscript-style report emphasizing inner-sector differentiation."""

    def first_row(df, var, metric):
        s = df[(df["Variable"] == var) & (df["Metric"] == metric)]
        return s.iloc[0] if len(s) else None

    def sec_means(var, metric):
        sub = metrics_df[
            (metrics_df["VARIABLE"] == var)
            & (metrics_df["ARC"] == "inner")
            & (metrics_df["INNER_SECTOR"].isin(INNER_SECTOR_ORDER))
        ].copy()
        if metric in ["mean_duration_d", "mean_peak_z"]:
            sub = sub[sub["n_events"] > 0]
        g = sub.groupby("INNER_SECTOR", dropna=False)[metric].mean()
        return {sec: float(g.get(sec, np.nan)) for sec in INNER_SECTOR_ORDER}

    sst_arc = first_row(arc_stat_df, "SST", "n_events_per_year")
    dli_arc = first_row(arc_stat_df, "DLI", "n_events_per_year")
    chl_arc = first_row(arc_stat_df, "Chl-a", "n_events_per_year")

    sst_kw = first_row(inner_kw_df, "SST", "n_events_per_year")
    dli_kw = first_row(inner_kw_df, "DLI", "n_events_per_year")
    chl_kw = first_row(inner_kw_df, "Chl-a", "n_events_per_year")

    dli_hab = (
        metrics_df[metrics_df["VARIABLE"] == "DLI"]
        .groupby("HAB", dropna=False)["n_events_per_year"]
        .mean()
        .sort_values(ascending=False)
    )
    dli_hab_txt = (
        ", ".join([f"{h}={v:.2f}" for h, v in dli_hab.items()])
        if len(dli_hab)
        else "NA"
    )

    dli_cov = metrics_df[metrics_df["VARIABLE"] == "DLI"][
        "COVERAGE_PCT_PRODUCT"
    ].dropna()
    chl_cov = metrics_df[metrics_df["VARIABLE"] == "Chl-a"][
        "COVERAGE_PCT_PRODUCT"
    ].dropna()

    sst_sec = sec_means("SST", "n_events_per_year")
    dli_sec = sec_means("DLI", "n_events_per_year")
    chl_sec = sec_means("Chl-a", "n_events_per_year")

    lines = []
    lines.append("# Extreme Event Analysis Report (Option A v2; SITE x HAB)")
    lines.append("")
    lines.append("## Methodology")
    lines.append(
        "Daily SST, benthic DLI, and chlorophyll-a records (2002-2008) were analyzed at the profile level, "
        "with each profile defined as a unique SITE x HAB combination. This design preserves habitat-specific "
        "depth effects, particularly for DLI, which was computed as PAR * exp(-KdPAR * depth). For each profile "
        "and variable, a day-of-year climatology (mean and standard deviation) was estimated and smoothed using a "
        "30-day circular window. Standardized anomalies were computed as z-scores, and extreme events were defined "
        f"as z >= {core.THRESHOLD_Z}, with a minimum duration of {core.MIN_DURATION_DAYS} days and event merging across "
        f"gaps <= {core.MERGE_GAP_DAYS} day when no missing values occurred. Event frequency, duration, peak intensity, "
        "cumulative intensity, and percentage of extreme days were computed per profile. Arc-level contrasts (Inner vs Outer) "
        "were assessed using permutation and Mann-Whitney tests, while Inner-sector contrasts (North, Central, South) "
        "were assessed using Kruskal-Wallis and pairwise permutation contrasts."
    )
    lines.append("")
    lines.append("## Results")

    if sst_arc is not None:
        lines.append(
            "At arc scale, SST event frequency remained higher in the Inner Arc "
            f"({sst_arc['Inner_Mean']:.2f} events year-1) than in the Outer Arc ({sst_arc['Outer_Mean']:.2f} events year-1; "
            f"permutation p={core.format_p_value(sst_arc['P_Value_Perm'])})."
        )
    if dli_arc is not None:
        lines.append(
            "DLI frequency was comparable between arcs "
            f"(Inner {dli_arc['Inner_Mean']:.2f} vs Outer {dli_arc['Outer_Mean']:.2f} events year-1; "
            f"p={core.format_p_value(dli_arc['P_Value_Perm'])}), indicating that cross-shelf differences in light extremes "
            "are weaker than profile-level heterogeneity."
        )
    if chl_arc is not None:
        lines.append(
            "Chl-a frequency exhibited a moderate tendency toward higher values in the Outer Arc "
            f"(Inner {chl_arc['Inner_Mean']:.2f} vs Outer {chl_arc['Outer_Mean']:.2f} events year-1; "
            f"p={core.format_p_value(chl_arc['P_Value_Perm'])})."
        )

    lines.append(
        "Inner-sector differentiation showed distinct patterns in frequency: "
        f"SST (North {sst_sec['inner_north']:.2f}, Central {sst_sec['inner_central']:.2f}, South {sst_sec['inner_south']:.2f}), "
        f"DLI (North {dli_sec['inner_north']:.2f}, Central {dli_sec['inner_central']:.2f}, South {dli_sec['inner_south']:.2f}), "
        f"and Chl-a (North {chl_sec['inner_north']:.2f}, Central {chl_sec['inner_central']:.2f}, South {chl_sec['inner_south']:.2f}) "
        "events year-1."
    )

    dli_outer_hab = (
        metrics_df[(metrics_df["VARIABLE"] == "DLI") & (metrics_df["ARC"] == "outer")]
        .groupby("HAB", dropna=False)["n_events_per_year"]
        .mean()
    )
    if len(dli_outer_hab) > 0:
        dli_outer_text = ", ".join([f"{h}={v:.2f}" for h, v in dli_outer_hab.items()])
        lines.append(
            "Within the Outer Arc, habitat-resolved DLI frequency was: "
            f"{dli_outer_text} events year-1, demonstrating explicit HAB-level contrasts within the same arc."
        )

    if sst_kw is not None:
        lines.append(
            f"Kruskal-Wallis tests for Inner-sector frequency yielded p={core.format_p_value(sst_kw['Kruskal_p'])} (SST), "
            f"p={core.format_p_value(dli_kw['Kruskal_p']) if dli_kw is not None else 'NA'} (DLI), and "
            f"p={core.format_p_value(chl_kw['Kruskal_p']) if chl_kw is not None else 'NA'} (Chl-a)."
        )

    lines.append(
        "Across all profiles, DLI frequency by habitat followed: "
        f"{dli_hab_txt} events year-1, reinforcing that habitat-specific depth structure must be retained in extreme-event analyses."
    )

    lines.append("")
    lines.append("## Limitations")
    lines.append(
        "Ocean-color products remained affected by substantial cloud-related missingness. Median product coverage among retained "
        f"profiles was {np.nanmedian(dli_cov):.1f}% for DLI and {np.nanmedian(chl_cov):.1f}% for Chl-a, which can reduce sensitivity "
        "for short-lived events and may influence profile inclusion."
    )
    lines.append(
        "SST and chlorophyll extractions are sampled from the same satellite pixel for multiple habitats within a site, whereas "
        "DLI differs explicitly with depth; therefore, mechanistic habitat contrasts are strongest for DLI."
    )
    lines.append(
        "The 2002-2008 window captures relevant interannual variability but remains finite for tail inference. Future work should "
        "include sensitivity analyses for threshold and duration choices and hierarchical models that explicitly account for within-site "
        "dependence among habitats."
    )

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n\n".join(lines) + "\n")


def run():
    print("=" * 78)
    print("EXTREME EVENTS ANALYSIS (OPTION A v2) - INNER SECTOR ENHANCED")
    print(f"Period: {core.START_YEAR}-{core.END_YEAR} | Threshold z={core.THRESHOLD_Z}")
    print(f"Output folder: {OUTPUT_DIR_V2}")
    print("=" * 78)

    profiles_df = core.load_profile_metadata(core.INTEGRATED_DATA_PATH)
    print(f"\nLoaded profiles (SITE x HAB): {len(profiles_df)}")

    print("\nSTEP 1: Listing source files")
    sst_files = core.filter_files_by_period(
        glob.glob(os.path.join(core.SST_DIR, core.SST_PATTERN)),
        core.START_YEAR,
        core.END_YEAR,
    )
    kd_files = core.filter_files_by_period(
        glob.glob(os.path.join(core.MODIS_DIR, core.KD490_PATTERN)),
        core.START_YEAR,
        core.END_YEAR,
    )
    par_files = core.filter_files_by_period(
        glob.glob(os.path.join(core.MODIS_DIR, core.PAR_PATTERN)),
        core.START_YEAR,
        core.END_YEAR,
    )
    chl_files = core.filter_files_by_period(
        glob.glob(os.path.join(core.MODIS_DIR, core.CHL_PATTERN)),
        core.START_YEAR,
        core.END_YEAR,
    )
    print(f"  SST files: {len(sst_files)}")
    print(f"  Kd490 files: {len(kd_files)}")
    print(f"  PAR files: {len(par_files)}")
    print(f"  Chl-a files: {len(chl_files)}")

    if (
        len(sst_files) == 0
        and len(kd_files) == 0
        and len(par_files) == 0
        and len(chl_files) == 0
    ):
        raise SystemExit("No NetCDF files found. Check source directories.")

    print("\nSTEP 2: Extracting profile-level time series")
    sst_ts = core.load_timeseries_for_profiles(
        sst_files,
        ["analysed_sst", "sea_surface_temperature", "sst"],
        profiles_df,
        "SST",
    )
    chl_ts = core.load_timeseries_for_profiles(
        chl_files, ["chlor_a"], profiles_df, "Chl-a"
    )
    dli_ts = core.calculate_dli_for_profiles(kd_files, par_files, profiles_df)

    kd_dates = {core.parse_date_from_filename(f) for f in kd_files}
    par_dates = {core.parse_date_from_filename(f) for f in par_files}
    kd_dates = {d for d in kd_dates if d is not None}
    par_dates = {d for d in par_dates if d is not None}
    expected_days_by_var = {
        "SST": len(sst_files),
        "DLI": len(kd_dates.intersection(par_dates)),
        "Chl-a": len(chl_files),
    }

    analysis_start = pd.Timestamp(year=core.START_YEAR, month=1, day=1)
    analysis_end = pd.Timestamp(year=core.END_YEAR, month=12, day=31)

    print("\nSTEP 3: Climatology, anomalies, and event detection")
    ts_by_var = {"SST": sst_ts, "DLI": dli_ts, "Chl-a": chl_ts}
    anomaly_dict = {"SST": {}, "DLI": {}, "Chl-a": {}}
    event_dict = {"SST": {}, "DLI": {}, "Chl-a": {}}

    profile_lookup = profiles_df.set_index("PROFILE_ID").to_dict("index")
    metrics_rows = []

    for var, ts_map in ts_by_var.items():
        kept = 0
        excluded = 0
        for pid, series in ts_map.items():
            s = series.sort_index().loc[analysis_start:analysis_end]
            n_valid = int(s.notna().sum())
            if n_valid < core.MIN_DAYS_FOR_CLIMATOLOGY:
                excluded += 1
                continue

            clm = core.compute_climatology(s)
            anom = core.compute_anomaly(s, clm)
            if anom is None:
                excluded += 1
                continue

            ev = core.detect_events(
                anom,
                threshold_z=core.THRESHOLD_Z,
                min_duration=core.MIN_DURATION_DAYS,
                merge_gap=core.MERGE_GAP_DAYS,
                analysis_start=analysis_start,
                analysis_end=analysis_end,
            )

            anomaly_dict[var][pid] = anom
            event_dict[var][pid] = ev
            kept += 1

            meta = profile_lookup.get(pid)
            if meta is None:
                continue

            metrics = core.compute_profile_metrics(
                ev, n_valid, expected_days_by_var[var]
            )
            metrics_rows.append(
                {
                    "PROFILE_ID": pid,
                    "SITE_ID": meta["SITE_ID"],
                    "HAB": meta["HAB"],
                    "ARC": meta["ARC"],
                    "INNER_SECTOR": meta["INNER_SECTOR"],
                    "REEF_NAME": meta["REEF_NAME"],
                    "LAT": meta["LAT"],
                    "LON": meta["LON"],
                    "DEPTH_M": meta["DEPTH_M"],
                    "VARIABLE": var,
                    **metrics,
                    "THRESHOLD_Z": core.THRESHOLD_Z,
                    "MIN_DURATION_DAYS": core.MIN_DURATION_DAYS,
                    "MERGE_GAP_DAYS": core.MERGE_GAP_DAYS,
                    "METHODOLOGY": "OptionA_anomaly_profile_level_v2",
                }
            )

        print(
            f"  {var}: retained={kept}, excluded(<{core.MIN_DAYS_FOR_CLIMATOLOGY} days)={excluded}, "
            f"events={sum(len(v) for v in event_dict[var].values())}"
        )

    metrics_df = pd.DataFrame(metrics_rows)
    if metrics_df.empty:
        raise SystemExit("No output metrics generated.")

    print("\nSTEP 4: Statistical summaries")
    arc_stat_df, sector_hab_df, zero_df = core.build_stat_tables(metrics_df)
    inner_kw_df, inner_pair_df = compute_inner_sector_stats(metrics_df)
    event_sector_hab_long_df = summarize_event_metrics_by_sector_hab(metrics_df)

    print("\nSTEP 5: Figures")
    create_composite_event_metrics_figure(
        metrics_df,
        output_path=os.path.join(
            OUTPUT_DIR_V2, "Fig1_Composite_Event_Metrics_Sector_HAB.png"
        ),
    )

    # Fig4 series: three overlaid variables with three independent y-axes
    create_fig4_overlay_pair(
        anomaly_dict,
        event_dict,
        "ITA1_TP",
        "ITA1_PA",
        "Inner North (ITA1): Extreme-event anomalies (DLI, SST, Chl-a)",
        os.path.join(OUTPUT_DIR_V2, "Fig4A_Inner_North_ITA1_TP_PA.png"),
    )
    create_fig4_overlay_pair(
        anomaly_dict,
        event_dict,
        "TIM2_TP",
        "TIM2_PA",
        "Inner Central (TIM2): Extreme-event anomalies (DLI, SST, Chl-a)",
        os.path.join(OUTPUT_DIR_V2, "Fig4B_Inner_Central_TIM2_TP_PA.png"),
    )
    create_fig4_overlay_pair(
        anomaly_dict,
        event_dict,
        "SGO_TP",
        "SGO_PA",
        "Inner South (SGO): Extreme-event anomalies (DLI, SST, Chl-a)",
        os.path.join(OUTPUT_DIR_V2, "Fig4C_Inner_South_SGO_TP_PA.png"),
    )
    create_fig4_overlay_pair(
        anomaly_dict,
        event_dict,
        "SIR_RR",
        "FAR_RR",
        "Outer Rocky Reefs (SIR vs FAR): Extreme-event anomalies (DLI, SST, Chl-a)",
        os.path.join(OUTPUT_DIR_V2, "Fig4D_Outer_Rocky_SIR_FAR_RR.png"),
    )
    create_fig4_overlay_pair(
        anomaly_dict,
        event_dict,
        "PAB4_TP",
        "PAB4_PA",
        "Outer PAB4 (TP vs PA): Extreme-event anomalies (DLI, SST, Chl-a)",
        os.path.join(OUTPUT_DIR_V2, "Fig4E_Outer_PAB4_TP_PA.png"),
    )

    cv_long_multi = load_cv_long_table_multivar()
    create_cv_hab_sector_figure_multivar(
        cv_long_multi,
        os.path.join(OUTPUT_DIR_V2, "Fig5_DLI_CV_HAB_Sectors.png"),
    )

    print("\nSTEP 6: Save outputs")
    metrics_path = os.path.join(
        OUTPUT_DIR_V2, "extreme_events_metrics_per_site_hab.csv"
    )
    metrics_df.to_csv(metrics_path, index=False)

    events_rows = []
    for var in ["SST", "DLI", "Chl-a"]:
        for pid, ev in event_dict[var].items():
            if ev is None or len(ev) == 0:
                continue
            meta = profile_lookup.get(pid)
            if meta is None:
                continue
            e = ev.copy().reset_index(drop=True)
            e["EVENT_ID"] = np.arange(1, len(e) + 1)
            e["PROFILE_ID"] = pid
            e["SITE_ID"] = meta["SITE_ID"]
            e["HAB"] = meta["HAB"]
            e["ARC"] = meta["ARC"]
            e["INNER_SECTOR"] = meta["INNER_SECTOR"]
            e["VARIABLE"] = var
            events_rows.append(e)

    events_df = (
        pd.concat(events_rows, ignore_index=True) if events_rows else pd.DataFrame()
    )
    events_path = os.path.join(
        OUTPUT_DIR_V2, "extreme_events_details_per_event_site_hab.csv"
    )
    events_df.to_csv(events_path, index=False)

    arc_path = os.path.join(OUTPUT_DIR_V2, "statistical_comparison_results_arc.csv")
    arc_stat_df.to_csv(arc_path, index=False)

    inner_kw_path = os.path.join(
        OUTPUT_DIR_V2, "statistical_comparison_results_inner_sectors.csv"
    )
    inner_kw_df.to_csv(inner_kw_path, index=False)

    inner_pair_path = os.path.join(OUTPUT_DIR_V2, "pairwise_inner_sector_tests.csv")
    inner_pair_df.to_csv(inner_pair_path, index=False)

    metric_sector_hab_path = os.path.join(
        OUTPUT_DIR_V2, "event_metrics_summary_by_sector_hab_long.csv"
    )
    event_sector_hab_long_df.to_csv(metric_sector_hab_path, index=False)

    sector_hab_path = os.path.join(OUTPUT_DIR_V2, "summary_by_sector_hab.csv")
    sector_hab_df.to_csv(sector_hab_path, index=False)

    zero_path = os.path.join(OUTPUT_DIR_V2, "zero_inflation_summary_site_hab.csv")
    zero_df.to_csv(zero_path, index=False)

    report_path = os.path.join(OUTPUT_DIR_V2, "Extreme_Events_Report.md")
    generate_report_v2(metrics_df, arc_stat_df, inner_kw_df, sector_hab_df, report_path)

    print(f"  Saved: {os.path.basename(metrics_path)}")
    print(f"  Saved: {os.path.basename(events_path)}")
    print(f"  Saved: {os.path.basename(arc_path)}")
    print(f"  Saved: {os.path.basename(inner_kw_path)}")
    print(f"  Saved: {os.path.basename(inner_pair_path)}")
    print(f"  Saved: {os.path.basename(metric_sector_hab_path)}")
    print(f"  Saved: {os.path.basename(sector_hab_path)}")
    print(f"  Saved: {os.path.basename(zero_path)}")
    print(f"  Saved: {os.path.basename(report_path)}")

    print("\n" + "=" * 78)
    print("ANALYSIS COMPLETE (v2)")
    print(f"Outputs: {OUTPUT_DIR_V2}")
    print("=" * 78)


if __name__ == "__main__":
    run()
