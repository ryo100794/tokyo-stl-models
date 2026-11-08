"""Explicit projected-metre to model-millimetre coordinates shared by models.

The public API always uses (easting, northing, height), independently of the
native axis order of a source CRS. ``mm_per_m`` is a unit conversion, not a
dimensionless map denominator: a 1:120000 model uses 1000 / 120000 mm/m.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class ModelCoordinates:
    origin_x: float
    origin_y: float
    origin_z: float = 0.0
    mm_per_m: float = 0.001
    vertical_exaggeration: float = 1.0
    z_offset_mm: float = 0.0
    crs: str = "EPSG:6677"

    def __post_init__(self):
        for name in ("origin_x", "origin_y", "origin_z", "mm_per_m",
                     "vertical_exaggeration", "z_offset_mm"):
            if not math.isfinite(float(getattr(self, name))):
                raise ValueError(f"{name} must be finite")
        if self.mm_per_m <= 0 or self.vertical_exaggeration <= 0:
            raise ValueError("mm_per_m and vertical_exaggeration must be positive")
        if not isinstance(self.crs, str) or not self.crs.strip():
            raise ValueError("crs must be a nonempty identifier")

    def point(self, point):
        """Convert one ENU point in projected metres into ENU millimetres."""
        x, y, z = map(float, point)
        return ((x - self.origin_x) * self.mm_per_m,
                (y - self.origin_y) * self.mm_per_m,
                (z - self.origin_z) * self.mm_per_m * self.vertical_exaggeration
                + self.z_offset_mm)

    def inverse(self, point):
        """Invert the model transform without changing its inherited height datum."""
        x, y, z = map(float, point)
        return (x / self.mm_per_m + self.origin_x,
                y / self.mm_per_m + self.origin_y,
                (z - self.z_offset_mm) / (self.mm_per_m * self.vertical_exaggeration)
                + self.origin_z)

    def points(self, points):
        import numpy as np
        array = np.asarray(points, dtype=np.float64)
        if array.ndim != 2 or array.shape[1] != 3 or not np.isfinite(array).all():
            raise ValueError("points must be a finite N by 3 array")
        return (array - [self.origin_x, self.origin_y, self.origin_z]) * [
            self.mm_per_m, self.mm_per_m, self.mm_per_m * self.vertical_exaggeration
        ] + [0, 0, self.z_offset_mm]

    def to_dict(self):
        data = asdict(self)
        data.update(schema_version="model_coordinates.1", projected_units="m",
                    model_units="mm", projected_axis_order=["easting", "northing", "height"],
                    model_axis_order=["east", "north", "up"], scale=self.mm_per_m)
        return data

    def write(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), ensure_ascii=False, indent=2) + "\n")

    @classmethod
    def from_dict(cls, data, *, mm_per_m=0.001, vertical_exaggeration=1.0,
                  z_offset_mm=0.0, crs="EPSG:6677"):
        """Read explicit fields or legacy origin_x/origin_y/scale JSON.

        Values stored in the file take precedence over fallback CLI defaults.
        Contradictory old/new scale fields are rejected, not silently selected.
        """
        if "mm_per_m" in data and "scale" in data and not math.isclose(
                float(data["mm_per_m"]), float(data["scale"]), rel_tol=1e-12, abs_tol=0):
            raise ValueError("origin JSON has inconsistent mm_per_m and legacy scale")
        if data.get("projected_units", "m") != "m" or data.get("model_units", "mm") != "mm":
            raise ValueError("origin JSON must use projected metres and model millimetres")
        if data.get("projected_axis_order", ["easting", "northing", "height"]) != ["easting", "northing", "height"]:
            raise ValueError("origin JSON public axes must be easting, northing, height")
        return cls(origin_x=float(data["origin_x"]), origin_y=float(data["origin_y"]),
            origin_z=float(data.get("origin_z", 0)),
            mm_per_m=float(data.get("mm_per_m", data.get("scale", mm_per_m))),
            vertical_exaggeration=float(data.get("vertical_exaggeration", vertical_exaggeration)),
            z_offset_mm=float(data.get("z_offset_mm", data.get("z_offset", z_offset_mm))),
            crs=data.get("crs", crs))

    @classmethod
    def read(cls, path, **fallbacks):
        return cls.from_dict(json.loads(Path(path).read_text()), **fallbacks)
