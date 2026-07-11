"""Download and cache central-region drive networks for the study cities.

Usage:
    python src/download.py [radius_m]

Graphs are saved as data/raw/{city}_{radius}m.graphml and re-downloads are
skipped if the file already exists.
"""

import sys
from pathlib import Path

import osmnx as ox

from cities import CITIES, DEFAULT_RADIUS_M

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"

ox.settings.use_cache = True
ox.settings.log_console = False


def graph_path(city: str, radius_m: int) -> Path:
    return DATA_DIR / f"{city}_{radius_m}m.graphml"


def download_city(city: str, radius_m: int = DEFAULT_RADIUS_M):
    """Download the drive network around the city's central point."""
    path = graph_path(city, radius_m)
    if path.exists():
        print(f"[skip] {city}: {path.name} already cached")
        return
    lat, lon = CITIES[city]
    print(f"[get ] {city}: r={radius_m} m around ({lat}, {lon}) ...", flush=True)
    G = ox.graph_from_point((lat, lon), dist=radius_m, network_type="drive",
                            simplify=True, retain_all=False)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    ox.save_graphml(G, path)
    print(f"[done] {city}: {len(G):,} nodes, {G.number_of_edges():,} edges")


def load_city(city: str, radius_m: int = DEFAULT_RADIUS_M):
    path = graph_path(city, radius_m)
    if not path.exists():
        download_city(city, radius_m)
    return ox.load_graphml(path)


if __name__ == "__main__":
    radius = int(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_RADIUS_M
    for name in CITIES:
        try:
            download_city(name, radius)
        except Exception as exc:  # keep going: Overpass hiccups on single cities
            print(f"[FAIL] {name}: {exc}")
