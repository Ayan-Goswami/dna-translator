from io import StringIO
from pathlib import Path

import pandas as pd
import streamlit as st
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

from dna_translator.core import clean_dna, find_orfs, translate_six_frames
from dna_translator.mutations import analyze_mutation
from dna_translator.properties import profile_protein
from dna_translator.viz import plot_hydropathy, plot_orf_map

st.set_page_config(
    page_title="DNA Translator v2.0",
    page_icon="🧬",
    layout="wide",
)

# --- Constants & Demo Presets ---
GENETIC_TABLES = {
    "Standard (Table 1)": 1,
    "Vertebrate Mitochondrial (Table 2)": 2,
    "Bacterial, Archaeal & Plant Plastid (Table 11)": 11,
}

DATA_DIR = Path(__file__).parent / "data"
FALLBACK_HBB_WT = (
    ">HBB_WT_Demo\n"
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGT"
    "GAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCT"
    "GATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCAC"
    "CTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTC"
    "AGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCC"
    "TATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
)


def load_file_or_fallback(filename: str, fallback: str) -> str:
    path = DATA_DIR / filename
    return path.read_text(encoding="utf-8") if path.exists() else fallback

@st.cache_data(show_spinner=False)
def cached_find_orfs(seq: str, min_aa: int, table_id: int, allow_alt: bool):
    return find_orfs(seq, min_aa_length=min_aa, table=table_id, allow_alt_starts=allow_alt)


@st.cache_data(show_spinner=False)
def cached_six_frames(seq: str, table_id: int):
    return translate_six_frames(seq, table=table_id)


@st.cache_data(show_spinner=False)
def cached_profile(protein_seq: str):
    return profile_protein(protein_seq)


@st.cache_data(show_spinner=False)
def cached_mutation(wt_seq: str, mut_seq: str, table_id: int):
    return analyze_mutation(wt_seq, mut_seq, table=table_id)

def parse_input_sequence(raw_text: str) -> tuple[str, str]:
    """Accept either FASTA format (>header + sequence) or plain raw DNA/RNA text."""
    text = raw_text.strip()
    if not text:
        raise ValueError("Please provide a DNA or RNA sequence.")
    if text.startswith(">"):
        records = list(SeqIO.parse(StringIO(text), "fasta"))
        if not records:
            raise ValueError("Could not parse FASTA header/sequence.")
        return records[0].id, clean_dna(str(records[0].seq))
    return "Custom_Sequence", clean_dna(text)


def make_sickle_demo(wt_fasta: str) -> str:
    """Derive Sickle Cell mutant (c.20A>T, p.Glu7Val) from the loaded HBB WT sequence."""
    _, seq = parse_input_sequence(wt_fasta)
    mut_seq = seq[:19] + "T" + seq[20:]
    return f">HBB_Sickle_Cell_p.Glu7Val\n{mut_seq}"


def make_nonsense_demo(wt_fasta: str) -> str:
    """Derive a premature stop mutation (c.19G>T, p.Glu7*) from HBB WT."""
    _, seq = parse_input_sequence(wt_fasta)
    mut_seq = seq[:18] + "T" + seq[19:]
    return f">HBB_Nonsense_p.Glu7Stop\n{mut_seq}"


def make_frameshift_demo(wt_fasta: str) -> str:
    """Derive a 1-bp deletion frameshift at position 20 from HBB WT."""
    _, seq = parse_input_sequence(wt_fasta)
    mut_seq = seq[:19] + seq[20:]
    return f">HBB_Frameshift_del20A\n{mut_seq}"


def build_alignment_html(wt_prot: str, mut_prot: str, chunk: int = 60) -> str:
    """Highlight mismatched amino acids between WT and Mutant proteins in monospace HTML."""
    max_len = max(len(wt_prot), len(mut_prot))
    wt_pad = wt_prot.ljust(max_len, "-")
    mut_pad = mut_prot.ljust(max_len, "-")

    lines = ["<div style='font-family: monospace; font-size: 13px; line-height: 1.6; overflow-x: auto;'>"]
    for start in range(0, max_len, chunk):
        wt_slice = wt_pad[start : start + chunk]
        mut_slice = mut_pad[start : start + chunk]
        mut_html = []
        for w, m in zip(wt_slice, mut_slice):
            if w == m:
                mut_html.append(m)
            else:
                mut_html.append(
                    f"<span style='background:#E11D48; color:#FFFFFF; font-weight:700; padding:0 2px; border-radius:3px;'>{m}</span>"
                )
        lines.append(f"<div><b>WT  {start + 1:>4}:</b> {wt_slice}</div>")
        lines.append(f"<div style='margin-bottom:8px;'><b>MUT {start + 1:>4}:</b> {''.join(mut_html)}</div>")
    lines.append("</div>")
    return "".join(lines)


# --- Header & Sidebar Controls ---
st.title("🧬 DNA Translator v2.0")
st.caption("6-Frame Open Reading Frame Finder • In Silico Variant Impact Simulator • Protein Physicochemical Profiler")

with st.sidebar:
    st.header("Analysis Parameters")
    table_label = st.selectbox("NCBI Genetic Code Table", list(GENETIC_TABLES.keys()))
    table_id = GENETIC_TABLES[table_label]
    min_aa = st.slider("Minimum ORF Length (aa)", min_value=10, max_value=200, value=30, step=5)
    allow_alt_starts = st.checkbox("Allow alternative start codons", value=False)
    kd_window = st.select_slider("Kyte-Doolittle Window (aa)", options=[5, 7, 9, 11, 13, 15], value=9)
    st.divider()
    st.markdown("**Built by Ayan Goswami**  \n[GitHub Repository](https://github.com/Ayan-Goswami/dna-translator)")

hbb_wt_default = load_file_or_fallback("hbb_coding.fasta", FALLBACK_HBB_WT)
hbb_sickle_default = load_file_or_fallback("hbb_sickle.fasta", make_sickle_demo(hbb_wt_default))

tab_orf, tab_mut = st.tabs(["🔎 6-Frame & ORF Explorer", "🧪 Mutation Impact Simulator"])

# ==========================================
# TAB 1: 6-FRAME TRANSLATION & ORF EXPLORER
# ==========================================
with tab_orf:
    col_in, col_upload = st.columns([3, 1])
    with col_upload:
        uploaded = st.file_uploader("Upload FASTA (.fasta, .fa, .fna)", type=["fasta", "fa", "fna", "txt"])
    with col_in:
        default_text = uploaded.getvalue().decode("utf-8") if uploaded else hbb_wt_default
        seq_input = st.text_area("DNA / RNA Sequence (FASTA or raw bases)", value=default_text, height=140)

    try:
        rec_id, cleaned_seq = parse_input_sequence(seq_input)
        orfs = cached_find_orfs(cleaned_seq, min_aa, table_id, allow_alt_starts)
        gc_pct = round(gc_fraction(cleaned_seq) * 100, 2)

        # Summary KPI row
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Sequence ID", rec_id[:22])
        m2.metric("Length", f"{len(cleaned_seq):,} nt")
        m3.metric("GC Content", f"{gc_pct}%")
        m4.metric(f"ORFs (≥ {min_aa} aa)", len(orfs))

        # Interactive 6-frame ORF track map
        st.plotly_chart(plot_orf_map(orfs, len(cleaned_seq)), use_container_width=True)

        if orfs:
            rows = []
            fasta_export_lines = []
            for idx, orf in enumerate(orfs, start=1):
                orf_name = f"ORF_{idx} ({orf.frame})"
                prof = cached_profile(orf.protein_seq)
                rows.append(
                    {
                        "ORF": orf_name,
                        "Frame": orf.frame,
                        "Start": orf.start,
                        "End": orf.end,
                        "Length (aa)": orf.aa_length,
                        "MW (kDa)": prof.molecular_weight_kda,
                        "pI": prof.isoelectric_point,
                        "GRAVY": prof.gravy,
                        "Instability": prof.instability_index,
                        "Stability": prof.stability_label,
                        "Protein": orf.protein_seq,
                    }
                )
                fasta_export_lines.append(
                    f">{rec_id}_ORF{idx} frame={orf.frame} coords={orf.start}-{orf.end} aa={orf.aa_length}\n{orf.protein_seq}"
                )

            df_orfs = pd.DataFrame(rows)
            st.subheader("Predicted ORFs & Physicochemical Summary")
            st.dataframe(df_orfs.drop(columns=["Protein"]), use_container_width=True, hide_index=True)

            # Export buttons
            dl1, dl2, _ = st.columns([1, 1, 3])
            dl1.download_button(
                "⬇️ Download CSV Report",
                data=df_orfs.to_csv(index=False).encode("utf-8"),
                file_name=f"{rec_id}_orfs.csv",
                mime="text/csv",
            )
            dl2.download_button(
                "⬇️ Download ORF Proteins (.fasta)",
                data="\n".join(fasta_export_lines).encode("utf-8"),
                file_name=f"{rec_id}_proteins.fasta",
                mime="text/plain",
            )

            # Single-ORF deep dive
            st.divider()
            selected_idx = st.selectbox(
                "Inspect Protein Profile & Hydropathy Curve",
                options=range(len(orfs)),
                format_func=lambda i: f"ORF #{i + 1} | Frame {orfs[i].frame} | {orfs[i].aa_length} aa ({orfs[i].start}–{orfs[i].end} nt)",
            )
            chosen_orf = orfs[selected_idx]
            chosen_prof = cached_profile(chosen_orf.protein_seq)

            p1, p2, p3, p4, p5 = st.columns(5)
            p1.metric("Molecular Weight", f"{chosen_prof.molecular_weight_kda} kDa")
            p2.metric("Isoelectric Point (pI)", chosen_prof.isoelectric_point)
            p3.metric("GRAVY", chosen_prof.gravy)
            p4.metric("Instability Index", f"{chosen_prof.instability_index} ({chosen_prof.stability_label})")
            p5.metric("Charge at pH 7.0", f"{chosen_prof.charge_at_ph7:+.2f}")

            if chosen_orf.aa_length >= kd_window:
                st.plotly_chart(
                    plot_hydropathy(chosen_orf.protein_seq, window=kd_window, wt_label=f"ORF #{selected_idx + 1}"),
                    use_container_width=True,
                )
            else:
                st.info(f"ORF length ({chosen_orf.aa_length} aa) is shorter than the selected hydropathy window ({kd_window} aa). Lower the window size in the sidebar to plot.")

            st.code(chosen_orf.protein_seq, language="text")
        else:
            st.warning(
                f"No ORFs of at least {min_aa} amino acids were found. Try lowering the minimum ORF length in the sidebar."
            )

        with st.expander("View Raw 6-Frame Translations"):
            frames = cached_six_frames(cleaned_seq, table_id)
            for label, seq_aa in frames.items():
                st.markdown(f"**Frame `{label}`** ({len(seq_aa):,} aa)")
                preview = (
                    seq_aa
                    if len(seq_aa) <= 5000
                    else f"{seq_aa[:5000]}\n... [truncated for display — {len(seq_aa):,} aa total]"
                )
                st.code(preview, language="text")

    except ValueError as err:
        st.error(f"Input Error: {err}")

# ==========================================
# TAB 2: IN SILICO MUTATION IMPACT SIMULATOR
# ==========================================
with tab_mut:
    preset = st.radio(
        "Quick Clinical & Functional Demos",
        options=[
            "Sickle Cell Disease (HBB c.20A>T, p.Glu7Val — Missense)",
            "Beta-Thalassemia Truncation (HBB c.19G>T, p.Glu7* — Nonsense)",
            "Frameshift Deletion (HBB del20A — Frameshift)",
        ],
        horizontal=True,
    )

    if "Sickle Cell" in preset:
        mut_preset_seq = hbb_sickle_default
    elif "Nonsense" in preset:
        mut_preset_seq = make_nonsense_demo(hbb_wt_default)
    else:
        mut_preset_seq = make_frameshift_demo(hbb_wt_default)

    c_wt, c_mut = st.columns(2)
    with c_wt:
        wt_input = st.text_area("Wild-Type CDS (FASTA or raw DNA)", value=hbb_wt_default, height=150)
    with c_mut:
        mut_input = st.text_area("Mutant CDS (FASTA or raw DNA)", value=mut_preset_seq, height=150)

    try:
        _, wt_clean = parse_input_sequence(wt_input)
        _, mut_clean = parse_input_sequence(mut_input)
        report = cached_mutation(wt_clean, mut_clean, table_id)

        if report.variant_type in {"Identical", "Synonymous (Silent)"}:
            st.success(f"**Variant Classification:** {report.variant_type}")
        elif report.variant_type in {"Missense", "In-frame Indel"}:
            st.warning(f"**Variant Classification:** {report.variant_type}")
        else:
            st.error(f"**Variant Classification:** {report.variant_type}")

        if report.changes:
            st.subheader("Detected Residue & Biochemical Shifts")
            df_changes = pd.DataFrame(
                [
                    {
                        "Position": ch.position,
                        "HGVS Protein": ch.hgvs_p,
                        "WT Residue": ch.wt_aa,
                        "Mutant Residue": ch.mut_aa,
                        "Biochemical Shift": ch.property_shift,
                    }
                    for ch in report.changes
                ]
            )
            st.dataframe(df_changes, use_container_width=True, hide_index=True)

        if report.wt_protein and report.mut_protein:
            wt_p = cached_profile(report.wt_protein)
            mut_p = cached_profile(report.mut_protein)

            st.subheader("Physicochemical Impact (WT → Mutant)")
            d1, d2, d3, d4, d5 = st.columns(5)
            d1.metric("Length", f"{mut_p.length} aa", delta=f"{mut_p.length - wt_p.length:+d} aa")
            d2.metric("MW (kDa)", f"{mut_p.molecular_weight_kda}", delta=f"{mut_p.molecular_weight_kda - wt_p.molecular_weight_kda:+.2f}")
            d3.metric("pI", f"{mut_p.isoelectric_point}", delta=f"{mut_p.isoelectric_point - wt_p.isoelectric_point:+.2f}")
            d4.metric("GRAVY", f"{mut_p.gravy}", delta=f"{mut_p.gravy - wt_p.gravy:+.3f}")
            d5.metric("Instability", f"{mut_p.instability_index}", delta=f"{mut_p.instability_index - wt_p.instability_index:+.2f}", delta_color="inverse")

            if len(report.wt_protein) >= kd_window and len(report.mut_protein) >= kd_window:
                st.plotly_chart(
                    plot_hydropathy(report.wt_protein, report.mut_protein, window=kd_window),
                    use_container_width=True,
                )
            else:
                st.info(
                    f"Mutant protein truncated to {len(report.mut_protein)} aa (shorter than window size {kd_window} aa), so hydropathy overlay is skipped."
                )

            st.subheader("Highlighted Protein Alignment (WT vs. Mutant)")
            st.markdown(build_alignment_html(report.wt_protein, report.mut_protein), unsafe_allow_html=True)

    except ValueError as err:
        st.error(f"Input Error: {err}")