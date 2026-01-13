"""Command line interface"""
import argparse
from pathlib import Path
import pandas as pd
from typing import List

from io_helpers import *
from model import CellModel, Location
from config import load_config, build_constraint

def parse_list(s: str) -> List[str]:
    """convert CLI Argument into a list"""    
    if not s:
        return None
    items = [x.strip() for x in s.split(",")]
    items = [x for x in items if x]
    return items or None

def run() -> int:
    p = argparse.ArgumentParser(
        prog="cell_layout",
        description=(
            "Constrained 3D layout of a whole-cell PPI graph by locations/organelles. "
            "Inputs: PPI edge list + localization table + JSON config. "
            "Outputs: per-location 3D coordinates + optional Plotly graph."
        ),
    )

    p.add_argument("--ppi", required=True, help="Path to 2-column edge list TSV/CSV")   #wprks
    p.add_argument("--localizations", required=True, help="Path to localization TSV")   #works
    p.add_argument("--config", required=True, help="Path to JSON config (see examples)")#works
    p.add_argument("--outdir", required=True, help="Output directory")                  #works

    p.add_argument("--only", help="Comma-separated list of location names to include")      #works
    p.add_argument("--exclude", help="Comma-separated list of location names to exclude")   #works

    p.add_argument("--seed", type=int, help="Override global seed")                     #does not work yet
    p.add_argument("--iterations", type=int, help="Override global iterations")         #does not work yet
    p.add_argument("--edge-strength", type=float, help="Override global edge_strength") #does not work yet, maybe get rid of it

    p.add_argument("--plot", action="store_true", help="Interactive Plotly graph") #works
    p.add_argument("--plot-title", default="3D cell", help="Plot title")           #works
    p.add_argument("--no-boundaries", action="store_true", help="Do not draw organelle boundaries") #does not work yet
    p.add_argument("--no-edges", action="store_true", help="Do not draw edges")                     #does not work

    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load config
    cfg = load_config(args.config)
    global_params = dict(cfg.global_params)
    if args.seed is not None:
        global_params["seed"] = args.seed
    if args.iterations is not None:
        global_params["iterations"] = args.iterations
    if args.edge_strength is not None:
        global_params["edge_strength"] = args.edge_strength
    
    only = parse_list(args.only)
    exclude = set(parse_list(args.exclude) or [])
    
    # Load data
    G = load_ppi_graph(args.ppi)
    print(G)
    loc_to_proteins, protein_to_locations, loc_df = load_localizations(args.localizations)
    init_locations_attribute(G)

    # Build model
    cell = CellModel(G)

    configured_names = []
    for spec in cfg.locations:
        if only is not None and spec.name not in only:
            continue
        if spec.name in exclude:
            continue

        raw = loc_to_proteins.get(spec.name, set())
        nodes = sorted(x for x in raw if isinstance(x, str) and x.strip() != "") # <- ignore NANs (occure because no loaction in locationfile)

        if not nodes:
            # Keep going, but warn into summary.
            pass

        constraint = build_constraint(spec.constraint_type, spec.constraint_params)
        loc = Location(
            name=spec.name,
            node_ids=nodes,
            min_degree=spec.min_degree,
            center=spec.center,
            constraint=constraint,
            repulsion_strength=spec.repulsion_strength,
        )
        cell.add_location(loc)
        configured_names.append(spec.name)

    if not configured_names:
        raise SystemExit("No locations selected. Check --only/--exclude and your config.")
    
    print(cell.locations)

    # Run layout
    cell.assign_all_positions()

    # Write positions (long format)
    rows = []
    for loc_name, loc in cell.locations.items():
        for node, (x, y, z) in loc.pos3d.items():
            #rows.append({"node": node, "location": loc_name, "x": x, "y": y, "z": z})
            rows.append({"node": node, "x": x, "y": y, "z": z})

    pos_path = outdir / "positions_per_location.tsv"
    print(pos_path)
    pd.DataFrame(rows).to_csv(pos_path, sep="\t", index=False)

    # Summary
    summary = cell.summary()
    
    '''
    # add selection warnings
    missing = [name for name in configured_names if name not in loc_to_proteins]
    summary["warnings"] = []
    if missing:
        summary["warnings"].append(
            {
                "type": "location_missing_in_localization_table",
                "locations": missing,
                "hint": "These location names were in the config but not found in the localization table.",
            }
        )
    '''
    #Plot
    if args.plot:
        #fig = cell.plot_all_locations_3d(title=args.plot_title)
        
        cell.plot_all_locations_3d(title=args.plot_title)
        
        #html_path = outdir / "cell_plot.html"
        #fig.write_html(html_path)

    return 0

def main() -> None:
    raise SystemExit(run())


if __name__ == "__main__":
    main()


