import pandas as pd

def read_diag_file(path, rename_func=lambda x: x.split("]")[1].split("(")[0]):
    """
    read the diagnostic file as a dataframe from the path

    file format is like this:
    ```
    step time ...
    0 0.0 ...
    ```
    """
    df = pd.read_csv(path, sep=" ")
    # rename the columns from something like '[2]total(J)' to 'total'
    return df.rename(columns=rename_func)


def plot_energy_evolution(
    field_energy_diag_path="diags/reducedfiles/field_energy.txt",
    part_energy_diag_path="diags/reducedfiles/part_energy.txt",
    meta: dict = dict(),
):
    field_energy = read_diag_file(field_energy_diag_path).rename(
        columns={"total_lev0": "E_field"}
    )
    part_energy = read_diag_file(part_energy_diag_path).rename(
        columns={"total": "E_part"}
    )

    energy_df = (
        field_energy.merge(part_energy)
        .assign(E_total=lambda x: x.E_field + x.E_part)
        .assign(time_norm=lambda x: x.time / meta["t_ci"])
    )
    print(energy_df.columns)

    ax = energy_df.plot(x="time_norm", y=["E_total", "E_field", "E_part"])
    ax.set_xlabel(r"time ($T_{ci}$)")
    ax.set_ylabel(r"Energy ($J$)")
    return ax