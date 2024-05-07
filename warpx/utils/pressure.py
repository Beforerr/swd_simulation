import polars as pl
import polars.selectors as cs
from utils import ds2df
import numpy as np
from tqdm import tqdm


def momentum2velocity(
    df: pl.DataFrame,
    mass,
    pfields=["particle_momentum_x", "particle_momentum_y", "particle_momentum_z"],
    vfields=["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"],
    drop=True,
):
    new_df = df.with_columns(
        pl.col(pfield).alias(vfield) / mass for pfield, vfield in zip(pfields, vfields)
    )
    return new_df.drop(pfields) if drop else new_df


def deposit_part(
    part_df: pl.DataFrame, field_df: pl.DataFrame, direction="z", breaks=None
):
    col = f"particle_position_{direction}"
    rename_map = {"brk": direction}
    index = f"{direction}_index"

    p_df = (
        part_df.drop(["particle_cpu", "particle_id", "particle_weight"])
        .with_columns(pl.col(col).cut(breaks, include_breaks=True))
        .unnest(col)
        .rename(rename_map)
        .drop(cs.by_dtype(pl.Categorical))
    )

    p_df = p_df.with_columns(pl.col(direction).rank("dense").alias(index)).drop(
        direction
    )
    field_df = field_df.with_columns(pl.col(direction).rank("dense").alias(index))

    return p_df.join(field_df, on=["time", index]).drop(index)


def calc_pressure(df: pl.DataFrame, mass: float, direction):
    #: Note: pressure caculation needs to account for particle weight and hypercube volume

    return (
        df.with_columns(
            (
                cs.contains("velocity") - cs.contains("velocity").mean().over(direction)
            ).pow(2)
        )
        .group_by(direction)
        .agg(
            cs.contains("velocity").mean().sqrt().name.map(
                lambda x: x.replace("particle_velocity", "velocity_th")
            ),
            (mass * cs.contains("velocity").sum()).name.map(
                lambda x: x.replace("particle_velocity", "pressure")
            ),
            (~cs.contains("velocity")).first(),
        )
    )


def calc_pressure_parp_perp(
    part_df: pl.DataFrame,
    field_df: pl.DataFrame,
    mass: float,
    direction="z",
    breaks=None,
    Bfields=["Bx", "By", "Bz"],
    vfields=["particle_velocity_x", "particle_velocity_y", "particle_velocity_z"],
):

    field_df = field_df.with_columns(Bmag=np.linalg.norm(field_df[Bfields], axis=1))
    p_mesh_df = deposit_part(part_df, field_df, direction=direction, breaks=breaks)

    particle_velocity_parp_expr = sum(
        pl.col(v_comp) * pl.col(B_comp) for v_comp, B_comp in zip(vfields, Bfields)
    ) / pl.col("Bmag")

    pressure_perp_expr = (
        pl.col("pressure_x")
        + pl.col("pressure_y")
        + pl.col("pressure_z")
        - pl.col("pressure_parp")
    ) / 2
    
    velocity_th_perp_expr = ((
        pl.col("velocity_th_x").pow(2)
        + pl.col("velocity_th_y").pow(2)
        + pl.col("velocity_th_z").pow(2)
        - pl.col("velocity_th_parp").pow(2)
    ) / 2).sqrt()

    return (
        p_mesh_df.with_columns(particle_velocity_parp=particle_velocity_parp_expr)
        .pipe(calc_pressure, mass=mass, direction=direction)
        .with_columns(
            pressure_perp=pressure_perp_expr,
            velocity_th_perp=velocity_th_perp_expr,
        )
        .with_columns(anisotropy=pl.col("pressure_perp") / pl.col("pressure_parp"))
    )


def calc_pressure_parp_perp_ds(
    ds_part,
    ds_field,
    meta,
    direction="z",
    breaks=None,
    Bfields=["Bx", "By", "Bz"],
):
    mass = meta["m_ion"]
    field_df = pl.DataFrame(ds2df(ds_field, fields=Bfields))
    part_df = pl.DataFrame(ds2df(ds_part)).pipe(momentum2velocity, mass=mass)

    i = 0 if ds_part.dimensionality == 1 else 2
    breaks = breaks or np.linspace(
        ds_part.domain_left_edge[i],
        ds_part.domain_right_edge[i],
        ds_part.domain_dimensions[i] + 1,
    )
    return calc_pressure_parp_perp(
        part_df, field_df, mass=mass, direction=direction, breaks=breaks
    )


def calc_pressure_parp_perp_ts(ts_part, ts_field, step=1, **kwargs):
    return pl.concat(
        calc_pressure_parp_perp_ds(ds_part, ds_field, **kwargs)
        for ds_part, ds_field in zip(tqdm(ts_part[::step]), ts_field[::step])
    )

def export_pressure_field(ts_part, ts_field, file="pressure.arrow", **kwargs):
    return calc_pressure_parp_perp_ts(ts_part, ts_field, **kwargs).write_ipc(file)