"""Top-10 most populous metro areas, with a hand-picked central point each.

Coordinates are (lat, lon) of a recognizable city-center landmark, so that
graph_from_point(radius) grabs the dense core rather than an administrative
centroid that may fall in a river or park.
"""

CITIES = {
    "tokyo": (35.6812, 139.7671),        # Tokyo Station
    "delhi": (28.6315, 77.2167),         # Connaught Place
    "shanghai": (31.2304, 121.4737),     # People's Square
    "dhaka": (23.7280, 90.4109),         # Motijheel / Gulistan
    "sao_paulo": (-23.5505, -46.6333),   # Praça da Sé
    "cairo": (30.0444, 31.2357),         # Tahrir Square
    "mexico_city": (19.4326, -99.1332),  # Zócalo
    "beijing": (39.9087, 116.3975),      # Tiananmen
    "mumbai": (18.9398, 72.8355),        # CST station
    "osaka": (34.7025, 135.4959),        # Osaka Station / Umeda
}

# 6 km download; with the 1 km boundary buffer the scored interior is a
# 5 km core, matching the ~5 km trip-locality radius from the literature
DEFAULT_RADIUS_M = 6000
