import plotly.graph_objects as go

from dna_translator.core import ORF
from dna_translator.properties import hydropathy_profile

FRAME_ORDER = ["+1", "+2", "+3", "-1", "-2", "-3"]


def plot_orf_map(orfs: list[ORF], seq_length: int) -> go.Figure:
    """Render an interactive 6-frame genome-track map of discovered ORFs."""
    fig = go.Figure()

    # Draw subtle baseline rails for all 6 reading frames
    for frame in FRAME_ORDER:
        fig.add_trace(
            go.Scatter(
                x=[1, seq_length],
                y=[frame, frame],
                mode="lines",
                line={"color": "rgba(148, 163, 184, 0.25)", "width": 2},
                hoverinfo="skip",
                showlegend=False,
            )
        )

    # Overlay each ORF as a horizontal bar positioned at its exact genomic coordinates
    shown_strands: set[str] = set()
    for idx, orf in enumerate(orfs, start=1):
        is_fwd = orf.strand == "+"
        color = "#0EA5E9" if is_fwd else "#F97316"
        strand_name = "Forward Strand (+)" if is_fwd else "Reverse Strand (-)"
        show_legend = strand_name not in shown_strands
        shown_strands.add(strand_name)

        fig.add_trace(
            go.Bar(
                orientation="h",
                y=[orf.frame],
                x=[orf.nt_length],
                base=[orf.start],
                width=0.5,
                marker={"color": color, "line": {"color": "#0F172A", "width": 1}},
                name=strand_name,
                legendgroup=strand_name,
                showlegend=show_legend,
                hovertemplate=(
                    f"<b>ORF #{idx} ({orf.frame})</b><br>"
                    f"Coordinates: {orf.start} – {orf.end} nt<br>"
                    f"Length: {orf.aa_length} aa ({orf.nt_length} nt)<br>"
                    f"Preview: {orf.protein_seq[:24]}...<extra></extra>"
                ),
            )
        )

    fig.update_layout(
        title="6-Frame Open Reading Frame (ORF) Map",
        barmode="overlay",
        xaxis_title="Nucleotide Position (5' → 3' Forward Reference)",
        yaxis_title="Reading Frame",
        xaxis={"range": [0, seq_length + max(10, int(seq_length * 0.03))]},
        yaxis={"categoryorder": "array", "categoryarray": list(reversed(FRAME_ORDER))},
        height=360,
        margin={"l": 40, "r": 20, "t": 50, "b": 40},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1},
        template="plotly_dark",
    )
    return fig


def plot_hydropathy(
    wt_protein: str,
    mut_protein: str | None = None,
    window: int = 9,
    wt_label: str = "Wild-Type",
    mut_label: str = "Mutant",
) -> go.Figure:
    """Plot Kyte-Doolittle hydropathy curve, optionally overlaying a mutant for comparison."""
    wt_pos, wt_scores = hydropathy_profile(wt_protein, window=window)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=wt_pos,
            y=wt_scores,
            mode="lines",
            name=wt_label,
            line={"color": "#38BDF8", "width": 2.5},
            hovertemplate="Residue %{x}<br>Score: %{y:.2f}<extra>" + wt_label + "</extra>",
        )
    )

    if mut_protein:
        mut_pos, mut_scores = hydropathy_profile(mut_protein, window=window)
        fig.add_trace(
            go.Scatter(
                x=mut_pos,
                y=mut_scores,
                mode="lines",
                name=mut_label,
                line={"color": "#FB7185", "width": 2.5, "dash": "dash"},
                hovertemplate="Residue %{x}<br>Score: %{y:.2f}<extra>" + mut_label + "</extra>",
            )
        )

    # Reference line at 0: positive = hydrophobic, negative = hydrophilic
    fig.add_hline(
        y=0,
        line_dash="dot",
        line_color="rgba(226, 232, 240, 0.5)",
        annotation_text="Hydrophobic (+) / Hydrophilic (-)",
        annotation_position="top right",
    )

    fig.update_layout(
        title=f"Kyte-Doolittle Hydropathy Profile (Window = {window} aa)",
        xaxis_title="Amino Acid Residue Position",
        yaxis_title="Hydropathy Score",
        height=360,
        margin={"l": 40, "r": 20, "t": 50, "b": 40},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1},
        template="plotly_dark",
    )
    return fig