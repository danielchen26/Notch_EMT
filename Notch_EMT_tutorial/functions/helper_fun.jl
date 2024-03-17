function run_and_check_new_vars(script_path::String)
    # Get the list of variables before running the code
    vars_before = names(Main)

    # Run the code
    include(script_path)

    # Get the list of variables after running the code
    vars_after = names(Main)

    # Find the newly created variables
    new_vars = setdiff(vars_after, vars_before)

    println("Newly created variables: ", new_vars)
end