using CSV, DataFrames, VegaLite
df = CSV.File("./Notch_params_complete.csv") |> DataFrame
df
