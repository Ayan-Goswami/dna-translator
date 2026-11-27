import json
from pathlib import Path
from typing import Annotated, Optional

import pandas as pd
import typer
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from dna_translator.core import find_orfs, translate_six_frames
from dna_translator.mutations import analyze_mutation
from dna_translator.properties import profile_protein

app = typer.Typer(
    help="DNA Translator v2.0: 6-frame ORF finder, mutation impact simulator, and protein profiler.",
    no_args_is_help=True,
)


def _load_fasta_records(fasta_path: Path) -> list[SeqRecord]:
    """Validate and stream FASTA records from disk."""
    if not fasta_path.exists():
        typer.secho(f"Error: File not found -> {fasta_path}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        typer.secho(f"Error: No valid FASTA records found in {fasta_path}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)
    return records


@app.command()
def orfs(
    fasta_file: Annotated[Path, typer.Argument(help="Input DNA/RNA FASTA file (single or multi-FASTA).")],
    min_aa: Annotated[int, typer.Option("--min-aa", "-m", help="Minimum ORF length in amino acids.")] = 30,
    table: Annotated[int, typer.Option("--table", "-t", help="NCBI genetic code table ID (1=Standard, 2=Vert Mito, 11=Bacterial).")] = 1,
    alt_starts: Annotated[bool, typer.Option("--alt-starts", help="Allow organism-specific alternative start codons.")] = False,
    csv_out: Annotated[Optional[Path], typer.Option("--csv", help="Export ORF + physicochemical table to CSV.")] = None,
    json_out: Annotated[Optional[Path], typer.Option("--json", help="Export results to JSON.")] = None,
    fasta_out: Annotated[Optional[Path], typer.Option("--fasta", help="Export translated ORF proteins as a FASTA file.")] = None,
) -> None:
    """Scan all 6 reading frames, extract ORFs, and profile each predicted protein."""
    records = _load_fasta_records(fasta_file)
    rows: list[dict] = []
    protein_records: list[SeqRecord] = []

    for rec in records:
        found = find_orfs(
            str(rec.seq),
            min_aa_length=min_aa,
            table=table,
            allow_alt_starts=alt_starts,
        )
        for idx, orf in enumerate(found, start=1):
            orf_id = f"{rec.id}_ORF{idx}"
            try:
                prof = profile_protein(orf.protein_seq)
                mw_kda = prof.molecular_weight_kda
                pi = prof.isoelectric_point
                gravy = prof.gravy
                instability = prof.instability_index
                stability = prof.stability_label
            except ValueError:
                mw_kda = pi = gravy = instability = None
                stability = "N/A"

            rows.append(
                {
                    "record_id": rec.id,
                    "orf_id": orf_id,
                    "frame": orf.frame,
                    "start": orf.start,
                    "end": orf.end,
                    "aa_len": orf.aa_length,
                    "mw_kda": mw_kda,
                    "pI": pi,
                    "gravy": gravy,
                    "instability": instability,
                    "stability": stability,
                    "protein_seq": orf.protein_seq,
                }
            )
            protein_records.append(
                SeqRecord(
                    Seq(orf.protein_seq),
                    id=orf_id,
                    description=f"frame={orf.frame} coords={orf.start}-{orf.end} aa={orf.aa_length} MW={mw_kda}kDa",
                )
            )

    if not rows:
        typer.secho(f"No ORFs >= {min_aa} aa found in {fasta_file.name}.", fg=typer.colors.YELLOW)
        raise typer.Exit(code=0)

    df = pd.DataFrame(rows)
    display_cols = ["orf_id", "frame", "start", "end", "aa_len", "mw_kda", "pI", "gravy", "stability"]
    typer.secho(f"\nFound {len(rows)} ORF(s) across {len(records)} FASTA record(s):\n", fg=typer.colors.GREEN, bold=True)
    typer.echo(df[display_cols].to_string(index=False))

    if csv_out:
        df.to_csv(csv_out, index=False)
        typer.secho(f"\nSaved CSV -> {csv_out}", fg=typer.colors.CYAN)

    if json_out:
        json_out.write_text(json.dumps(rows, indent=2), encoding="utf-8")
        typer.secho(f"Saved JSON -> {json_out}", fg=typer.colors.CYAN)

    if fasta_out:
        SeqIO.write(protein_records, fasta_out, "fasta")
        typer.secho(f"Saved Protein FASTA -> {fasta_out}", fg=typer.colors.CYAN)


@app.command()
def translate(
    fasta_file: Annotated[Path, typer.Argument(help="Input DNA/RNA FASTA file.")],
    table: Annotated[int, typer.Option("--table", "-t", help="NCBI genetic code table ID.")] = 1,
    fasta_out: Annotated[Optional[Path], typer.Option("--fasta", help="Export all 6 frames to a protein FASTA file.")] = None,
) -> None:
    """Translate every sequence in a FASTA file across all 6 reading frames."""
    records = _load_fasta_records(fasta_file)
    out_records: list[SeqRecord] = []

    for rec in records:
        frames = translate_six_frames(str(rec.seq), table=table)
        typer.secho(f"\n=== {rec.id} ({len(rec.seq)} nt) ===", fg=typer.colors.GREEN, bold=True)
        for frame_label, prot in frames.items():
            preview = prot if len(prot) <= 70 else f"{prot[:70]}... ({len(prot)} aa)"
            typer.echo(f"  {frame_label:>2} | {preview}")
            out_records.append(
                SeqRecord(Seq(prot), id=f"{rec.id}_frame{frame_label}", description=f"table={table}")
            )

    if fasta_out:
        SeqIO.write(out_records, fasta_out, "fasta")
        typer.secho(f"\nSaved 6-frame FASTA -> {fasta_out}", fg=typer.colors.CYAN)


@app.command()
def mutate(
    wt_fasta: Annotated[Path, typer.Argument(help="Wild-type coding sequence (CDS) FASTA file.")],
    mut_fasta: Annotated[Path, typer.Argument(help="Mutant coding sequence (CDS) FASTA file.")],
    table: Annotated[int, typer.Option("--table", "-t", help="NCBI genetic code table ID.")] = 1,
) -> None:
    """Compare a wild-type CDS against a mutant CDS and report functional + biochemical shifts."""
    wt_rec = _load_fasta_records(wt_fasta)[0]
    mut_rec = _load_fasta_records(mut_fasta)[0]

    report = analyze_mutation(str(wt_rec.seq), str(mut_rec.seq), table=table)

    typer.secho(f"\nVariant Classification: {report.variant_type}", fg=typer.colors.CYAN, bold=True)
    typer.echo(f"WT Length:  {len(report.wt_protein)} aa  ({wt_rec.id})")
    typer.echo(f"Mut Length: {len(report.mut_protein)} aa  ({mut_rec.id})  [Delta: {report.length_delta_aa:+d} aa]")

    if report.changes:
        typer.secho("\nResidue Changes:", bold=True)
        for ch in report.changes:
            typer.echo(f"  • {ch.hgvs_p}  |  {ch.property_shift}")

    if report.wt_protein and report.mut_protein:
        wt_prof = profile_protein(report.wt_protein)
        mut_prof = profile_protein(report.mut_protein)
        typer.secho("\nPhysicochemical Delta (WT -> Mutant):", bold=True)
        typer.echo(f"  • MW (kDa):    {wt_prof.molecular_weight_kda} -> {mut_prof.molecular_weight_kda}")
        typer.echo(f"  • pI:          {wt_prof.isoelectric_point} -> {mut_prof.isoelectric_point}")
        typer.echo(f"  • GRAVY:       {wt_prof.gravy} -> {mut_prof.gravy}")
        typer.echo(f"  • Instability: {wt_prof.instability_index} ({wt_prof.stability_label}) -> {mut_prof.instability_index} ({mut_prof.stability_label})")


if __name__ == "__main__":
    app()