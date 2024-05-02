import polars as pl
import polars.selectors as cs
from utils import ds2df
import numpy as np


def calc_pressure(ds, direction, breaks=None):
    df = pl.DataFrame(ds2df(ds))
    
    i = 0 if ds.dimensionality == 1 else 2
    breaks = breaks or np.linspace(
        ds.domain_left_edge[i], ds.domain_right_edge[i], ds.domain_dimensions[i] + 1
    )
    col = f"particle_position_{direction}"

    rename_map = {
        "brk": direction,
        "particle_momentum_x_var": "p_xx",
        "particle_momentum_y_var": "p_yy",
        "particle_momentum_z_var": "p_zz",
    }

    return (
        df.with_columns(
            # pl.col(f"particle_position_{direction}").qcut(n).alias(direction),
            pl.col(col).cut(breaks, include_breaks=True),
        )
        .with_columns(
            (cs.contains("momentum") - cs.contains("momentum").mean().over(col))
            .pow(2)
            .name.suffix("_var"),
            # cs.contains("momentum").mean().over(direction).name.suffix("_mean"),
        )
        .group_by(col)
        .agg(cs.contains("var").sum(), pl.col("time").first())
        .unnest(col)
        .rename(rename_map)
        .drop(cs.by_dtype(pl.Categorical))
    )

#TODO
def calc_pressure_parp_perp(df: pl.DataFrame):
    return (
        df.with_columns(Bmag=np.linalg.norm(b_df[Bfields], axis=1))
        .with_columns(
            p_parp=(
                pl.col("p_xx") * pl.col("Bx").abs()
                + pl.col("p_yy") * pl.col("By").abs()
                + pl.col("p_zz") * pl.col("Bz").abs()
            )
            / pl.col("Bmag")
        )
        .with_columns(
            p_perp=(pl.col("p_xx") + pl.col("p_yy") + pl.col("p_zz") - pl.col("p_parp"))
            / 2
        )
        .with_columns(anisotropy=pl.col("p_perp") / pl.col("p_parp"))
    )
    
def export_pressure_field(ts_part, ts_field, step=8):
    direction = "z"
    Bfields = ["Bx", "By", "Bz"]

    p_df = pl.concat(calc_pressure(ds, direction) for ds in ts_part[::step])
    b_df = pl.concat(pl.DataFrame(ds2df(ds, fields=Bfields)) for ds in ts_field[::step])

    index = f"{direction}_index"
    p_df = p_df.with_columns(pl.col(direction).rank("dense").alias(index))
    b_df = b_df.with_columns(pl.col(direction).rank("dense").alias(index))
    df = b_df.join(p_df, on=["time", index]).drop(f"{direction}_right", index).drop(Bfields)
    df.write_ipc("pressure.arrow")