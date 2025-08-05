using DataFrames, CSV

# Test loading from db_complete.csv
println("Loading from Data/parameter_databases/db_complete.csv:")
db1 = CSV.File("Data/parameter_databases/db_complete.csv") |> DataFrame
println("Column names: ", names(db1))
println("Number of rows: ", nrow(db1))

# Test loading from the correct file
println("\nLoading from ../Notch_EMT_data/Notch_params_complete.csv:")
db2 = CSV.File("../Notch_EMT_data/Notch_params_complete.csv") |> DataFrame
println("Column names: ", names(db2))
println("Number of rows: ", nrow(db2))
println("\nFirst few parameter columns for row 49:")
println("k0 = ", db2[49, :k0])
println("k1 = ", db2[49, :k1])
println("k2 = ", db2[49, :k2])