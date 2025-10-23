using Pkg

# Function to read dependencies.jl and extract package names from 'using' statements
function extract_packages_from_dependencies(dependencies_file)
    lines = readlines(dependencies_file)
    packages = String[]

    for line in lines
        line = strip(line)
        if startswith(line, "using ")
            parts = split(line)
            pkg = parts[2]
            push!(packages, pkg)
        end
    end

    return packages
end

# Function to install missing packages (modern version)
function install_packages_from_dependencies(dependencies_file)
    packages = extract_packages_from_dependencies(dependencies_file)

    # Get the list of installed packages in the active environment
    installed_pkgs = keys(Pkg.dependencies())

    for pkg in packages
        if !(pkg in installed_pkgs)
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

# Run the function
install_packages_from_dependencies(pwd() * "/dependencies.jl")