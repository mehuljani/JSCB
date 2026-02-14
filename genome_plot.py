import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
from reportlab.graphics.shapes import String

def calculate_gc_stats(sequence, window_size=2000, step_size=2000):
    """Calculates GC Content and GC Skew in a sliding window."""
    gc_content, gc_skew = [], []
    for i in range(0, len(sequence), step_size):
        subseq = sequence[i:i+window_size]
        if not subseq: break
        g, c = subseq.count('G') + subseq.count('g'), subseq.count('C') + subseq.count('c')
        a, t = subseq.count('A') + subseq.count('a'), subseq.count('T') + subseq.count('t')
        total = g + c + a + t
        if total > 0:
            gc_content.append((i, (g + c) / total))
            gc_skew.append((i, (g - c) / (g + c) if (g + c) > 0 else 0))
    return gc_content, gc_skew

def main():
    parser = argparse.ArgumentParser(
        description="Visualize Genomic Islands from JSCB output on a circular genome map (PDF/SVG).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required Arguments
    parser.add_argument("-g", "--genbank", required=True, help="Path to the GenBank file.")
    parser.add_argument("-t", "--tsv", required=True, help="Path to the JSCB_output.tsv file.")
    
    # Optional Customization
    parser.add_argument("-o", "--output", default="genome_map", help="Base name for output files (saves .pdf and .svg).")
    parser.add_argument("-n", "--name", help="Override the genome name displayed in the center.")
    parser.add_argument("--color", default="firebrick", help="Color for Genomic Islands (e.g., 'blue', 'red', 'green').")
    parser.add_argument("--window", type=int, default=2000, help="Window size for GC calculation.")
    parser.add_argument("--core", type=float, default=0.5, help="Size of the inner circle core (0.1 to 0.9).")
    parser.add_argument("--no-labels", action="store_true", help="Disable island ID labels on the plot.")
    parser.add_argument("--label-size", type=int, default=25, help="Font size for island labels.")

    args = parser.parse_args()

    # 1. Load Data
    try:
        record = SeqIO.read(args.genbank, "genbank")
        genome_length = len(record)
        display_name = args.name if args.name else record.id
        print(f"Loaded {display_name}: {genome_length:,} bp")
        
        # Read TSV and identify columns robustly
        df = pd.read_csv(args.tsv, sep='\t')
        id_col = [c for c in df.columns if 'Id' in c][0]
        start_col = [c for c in df.columns if 'Start' in c][0]
        end_col = [c for c in df.columns if 'End' in c][0]
        print(f"Processing {len(df)} islands from {args.tsv}")
    except Exception as e:
        print(f"Error loading files: {e}")
        sys.exit(1)

    # 2. Plotting Logic
    print("Calculating GC statistics...")
    gc_data, skew_data = calculate_gc_stats(str(record.seq), window_size=args.window, step_size=args.window)
    avg_gc = sum(y for x, y in gc_data) / len(gc_data) if gc_data else 0.5

    gd_diagram = GenomeDiagram.Diagram(display_name)
    
    # Track 1: GC Skew (Inner)
    set_skew = gd_diagram.new_track(1, name="GC Skew", height=0.5).new_set(type="graph")
    set_skew.new_graph(skew_data, "GC Skew", style="bar", color=colors.purple, altcolor=colors.orange, center=0)

    # Track 2: GC Content (Middle)
    set_gc = gd_diagram.new_track(2, name="GC Content", height=0.5).new_set(type="graph")
    set_gc.new_graph(gc_data, "GC Content", style="line", color=colors.black, center=avg_gc)

    # Track 3: Genomic Islands (Outer)
    set_islands = gd_diagram.new_track(3, name="Islands", height=0.5).new_set()
    
    # Map the color string to a ReportLab color object
    island_color = getattr(colors, args.color.lower(), colors.firebrick)
    
    for _, row in df.iterrows():
        feature = SeqFeature(FeatureLocation(int(row[start_col]), int(row[end_col])), type="island")
        set_islands.add_feature(feature, color=island_color, label=(not args.no_labels), 
                                name=str(row[id_col]), label_size=args.label_size, sigil="BOX")

    # 3. Draw and Add Center Annotations
    size = 1500
    gd_diagram.draw(format="circular", circular=True, pagesize=(size, size), 
                    start=0, end=genome_length, circle_core=args.core)
    
    drawing = gd_diagram.drawing
    c_x, c_y = size / 2, size / 2
    
    # Add Genome Name
    drawing.add(String(c_x, c_y + 35, display_name, textAnchor="middle", fontSize=30, fontName="Helvetica-Bold"))
    # Fixed coordinate and text for genome length
    drawing.add(String(c_x, c_y, f"{genome_length:,} bp", textAnchor="middle", fontSize=22))
    # Add Average GC
    drawing.add(String(c_x, c_y - 30, f"Avg GC: {avg_gc*100:.1f}%", textAnchor="middle", fontSize=20))

    # 4. Save Outputs
    gd_diagram.write(f"{args.output}.pdf", "PDF")
    gd_diagram.write(f"{args.output}.svg", "SVG")
    print(f"Success! Map saved as:\n - {args.output}.pdf\n - {args.output}.svg")

if __name__ == "__main__":
    main()