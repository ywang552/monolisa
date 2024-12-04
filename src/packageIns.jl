using Pkg

# Function to read dependencies.jl and extract package names from 'using' statements without using regex
function extract_packages_from_dependencies(dependencies_file)
    # Read the file line by line
    lines = readlines(dependencies_file)
    
    # List to store the package names
    packages = String[]
    
    # Iterate over each line
    for line in lines
        # Remove leading/trailing whitespace
        line = strip(line)
        
        # Check if the line starts with 'using' and has a package name
        if startswith(line, "using ")
            # Extract the part after 'using' (the package name)
            # Split by spaces and take the second part
            parts = split(line)
            push!(packages, parts[2])
        end
    end
    
    return packages
end

# Function to install missing packages
function install_packages_from_dependencies(dependencies_file)
    # Get the list of package names from dependencies.jl
    packages = extract_packages_from_dependencies(dependencies_file)

    # Loop through each package and install if not already installed
    for pkg in packages
        if !(pkg in keys(Pkg.installed()))
            println("Installing package: $pkg")
            try
                Pkg.add(pkg)
            catch e
                println("Failed to install $pkg: $e")
            end
        else
            println("Package $pkg is already installed")
        end
    end
end

# Run the function with the path to your dependencies.jl file
install_packages_from_dependencies(pwd()*"/src/dependencies.jl")
