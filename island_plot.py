import sys
import argparse
import pandas as pd
import math
import os
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors

def get_gene_color(feature):
    """Determines the color of a gene based on its product or name."""
    product = " ".join(feature.qualifiers.get("product", [])).lower()
    gene_name = " ".join(feature.qualifiers.get("gene", [])).lower()
    note = " ".join(feature.qualifiers.get("note", [])).lower()
    
    full_text = f"{product} {gene_name} {note}"
    
    # 1. Check for mobile elements (Transposases, etc.) -> Red
    mobile_keywords = ["transposase", "integrase", "excisionase", "transposon", "insertion element", "phage", "recombinase", "integrations host factor", "plasmid"]
    if any(keyword in full_text for keyword in mobile_keywords):
        return colors.red
    
    # 2. Check for hypothetical proteins -> Gray
    if "hypothetical protein" in product:
        return colors.lightgrey
    
    # 3. Default for all other genes -> Lightblue (Strand-neutral)
    return colors.lightblue

def plot_single_island(island_id, island_start, island_end, island_genes, args, output_dir):
    """Creates the plot for a single genomic island."""
    gd_diagram = GenomeDiagram.Diagram(f"Island {island_id}")
    total_len = island_end - island_start
    chunk_size = math.ceil(total_len / args.rows)
    
    for i in range(args.rows):
        row_start = island_start + (i * chunk_size)
        row_end = min(row_start + chunk_size, island_end)
        
        track = gd_diagram.new_track(args.rows - i, name="", start=row_start, end=row_end, 
                                     greytrack=False, height=0.5)
        track.axis = True
        track.axis_labels = True
        feature_set = track.new_set()

        for gene in island_genes:
            if gene.location.end > row_start and gene.location.start < row_end:
                gene_color = get_gene_color(gene)
                label = gene.qualifiers.get("gene", gene.qualifiers.get("locus_tag", [""] ))[0]
                
                feature_set.add_feature(gene, sigil="ARROW", color=gene_color, 
                                        label=True, label_size=args.label_size, 
                                        label_angle=90, name=label,
                                        height=args.arrow_height,
                                        arrowshaft_height=args.shaft_height,
                                        arrowhead_length=args.head_length)

    base_out = os.path.join(output_dir, f"island_{island_id}")
    page_width, page_height = 1200, 200 * args.rows + 100
    
    gd_diagram.draw(format="linear", orientation="landscape", pagesize=(page_width, page_height),
                    fragments=1, start=island_start, end=island_end)

    if args.format in ['pdf', 'all']:
        gd_diagram.write(f"{base_out}.pdf", "PDF")
    if args.format in ['svg', 'all']:
        gd_diagram.write(f"{base_out}.svg", "SVG")

def main():
    parser = argparse.ArgumentParser(
        description="Visualize genomic islands with uniform gene coloring (except mobile/hypothetical).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-g", "--genbank", required=True, help="Path to the GenBank file.")
    parser.add_argument("-t", "--tsv", required=True, help="Path to the JSCB_output.tsv file.")
    parser.add_argument("-d", "--dir", default="island_plots", help="Directory to save the plots.")
    parser.add_argument("--rows", type=int, default=1, help="Number of rows to split each island into.")
    parser.add_argument("--format", choices=['pdf', 'svg', 'all'], default='all', help="Output format.")
    parser.add_argument("--label-size", type=int, default=8, help="Font size for labels.")
    
    # Sharp Arrow Defaults
    parser.add_argument("--arrow-height", type=float, default=0.3, help="Height of the arrow head.")
    parser.add_argument("--shaft-height", type=float, default=0.15, help="Height of the arrow tail.")
    parser.add_argument("--head-length", type=float, default=0.2, help="Proportional length of the arrow head.")

    args = parser.parse_args()

    if not os.path.exists(args.dir):
        os.makedirs(args.dir)

    try:
        df = pd.read_csv(args.tsv, sep='\t')
        id_col = [c for c in df.columns if 'Id' in c][0]
        start_col = [c for c in df.columns if 'Start' in c][0]
        end_col = [c for c in df.columns if 'End' in c][0]
    except Exception as e:
        print(f"Error loading TSV: {e}"); sys.exit(1)

    try:
        print("Reading GenBank file...")
        record = SeqIO.read(args.genbank, "genbank")
        all_cds = [f for f in record.features if f.type == "CDS"]
    except Exception as e:
        print(f"Error loading GenBank: {e}"); sys.exit(1)

    print(f"Processing {len(df)} islands...")
    for _, row in df.iterrows():
        island_id = str(row[id_col])
        island_start, island_end = int(row[start_col]), int(row[end_col])
        island_genes = [f for f in all_cds if f.location.start >= island_start and f.location.end <= island_end]
        
        print(f" -> Generating {island_id}")
        plot_single_island(island_id, island_start, island_end, island_genes, args, args.dir)

    print(f"\nDone! All plots saved in '{args.dir}'.")

if __name__ == "__main__":
    main()