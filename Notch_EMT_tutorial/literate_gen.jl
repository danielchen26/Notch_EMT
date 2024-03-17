using Revise
using Literate
current_folder_path = @__DIR__
# ===== Generate Markdown ======
Literate.markdown(
joinpath(current_folder_path, "Notch_EMT_main.jl"),
    current_folder_path,
    documenter=true)


# ===== Generate Notebook =======
Literate.notebook(
joinpath(current_folder_path, "Notch_EMT_main.jl"),
    current_folder_path,
    execute=true,
    documenter=true)