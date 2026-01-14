from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl


def _median_step(values: np.ndarray) -> Optional[float]:
    if values.size < 2:
        return None
    uniq = np.unique(values)
    if uniq.size < 2:
        return None
    diffs = np.diff(np.sort(uniq))
    diffs = diffs[diffs > 0]
    if diffs.size == 0:
        return None
    return float(np.median(diffs))


def estimate_grid_radius(df: pl.DataFrame) -> Optional[float]:
    if df.height < 2:
        return None
    x = df["x_mtx"].to_numpy()
    y = df["y_mtx"].to_numpy()
    dx = _median_step(x)
    dy = _median_step(y)
    steps = [v for v in [dx, dy] if v is not None and v > 0]
    if not steps:
        return None
    base = max(steps)
    return float(base * 2.1)


class NeighborIndex:
    def __init__(
        self,
        positions_df: pl.DataFrame,
        default_mode: str = "radius",
        radius: Optional[float] = None,
        k: Optional[int] = None,
    ):
        self.positions = positions_df.select(["barcode", "x_mtx", "y_mtx"])
        self.default_mode = default_mode
        self.default_radius = float(radius) if radius is not None else estimate_grid_radius(self.positions)
        self.default_k = int(k) if k is not None else None

        self._pairs = None
        self._pairs_radius = None
        self._neighbors_cache: Dict[Tuple[Optional[float], Optional[int], bool], Dict[str, List[str]]] = {}

    def _build_pairs(self, radius: float) -> pl.DataFrame:
        if radius <= 0:
            return pl.DataFrame({"barcode": [], "neighbor": [], "distance": []})
        cell_size = float(radius)
        left = self.positions.with_columns(
            (pl.col("x_mtx") / cell_size).floor().cast(pl.Int64).alias("cell_x"),
            (pl.col("y_mtx") / cell_size).floor().cast(pl.Int64).alias("cell_y"),
        )
        offsets = pl.DataFrame(
            {
                "dx": [-1, 0, 1, -1, 0, 1, -1, 0, 1],
                "dy": [-1, -1, -1, 0, 0, 0, 1, 1, 1],
            }
        )
        candidates = left.join(offsets, how="cross").with_columns(
            (pl.col("cell_x") + pl.col("dx")).alias("cell_x_n"),
            (pl.col("cell_y") + pl.col("dy")).alias("cell_y_n"),
        )
        right = left.select(
            pl.col("barcode").alias("neighbor"),
            pl.col("x_mtx").alias("nx"),
            pl.col("y_mtx").alias("ny"),
            pl.col("cell_x"),
            pl.col("cell_y"),
        )
        pairs = candidates.join(
            right, left_on=["cell_x_n", "cell_y_n"], right_on=["cell_x", "cell_y"], how="inner"
        ).with_columns(
            ((pl.col("x_mtx") - pl.col("nx")) ** 2 + (pl.col("y_mtx") - pl.col("ny")) ** 2)
            .sqrt()
            .alias("distance")
        )
        pairs = pairs.filter(pl.col("distance") <= radius)
        return pairs.select(["barcode", "neighbor", "distance"])

    def _ensure_pairs(self, radius: Optional[float]) -> None:
        if radius is None:
            return
        if self._pairs is None or self._pairs_radius is None or radius > self._pairs_radius:
            self._pairs = self._build_pairs(radius)
            self._pairs_radius = radius
            self._neighbors_cache.clear()

    def _neighbors_map(
        self,
        radius: Optional[float],
        k: Optional[int],
        include_center: bool,
    ) -> Dict[str, List[str]]:
        cache_key = (radius, k, include_center)
        if cache_key in self._neighbors_cache:
            return self._neighbors_cache[cache_key]

        if radius is None:
            return {}
        self._ensure_pairs(radius)
        pairs = self._pairs
        if not include_center:
            pairs = pairs.filter(pl.col("neighbor") != pl.col("barcode"))
        pairs = pairs.sort(["barcode", "distance"])
        if k is not None:
            grouped = pairs.group_by("barcode", maintain_order=True).agg(pl.col("neighbor").head(k))
        else:
            grouped = pairs.group_by("barcode", maintain_order=True).agg(pl.col("neighbor"))
        out = {b: n for b, n in zip(grouped["barcode"].to_list(), grouped["neighbor"].to_list())}
        self._neighbors_cache[cache_key] = out
        return out

    def neighbors_by_barcode(
        self,
        barcode: str,
        radius: Optional[float] = None,
        k: Optional[int] = None,
        include_center: bool = False,
    ) -> List[str]:
        if radius is None and k is None:
            mode = self.default_mode
            if mode == "k":
                k = self.default_k or 8
            else:
                radius = self.default_radius
        elif radius is None:
            radius = self.default_radius

        if radius is None:
            return []
        neighbors = self._neighbors_map(radius, k, include_center)
        return neighbors.get(barcode, [])
