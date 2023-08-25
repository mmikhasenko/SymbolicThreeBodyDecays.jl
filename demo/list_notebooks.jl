using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
Pkg.add("Pluto")
# 
using Pluto

# make sure the script is located in the notebooks folder
# the script assumes it is a subfolder of the project
notebooks_dir = "demo"
@assert splitpath(@__DIR__)[end] == notebooks_dir

files = readdir(@__DIR__)
notebooks = filter(files) do f
    Pluto.is_pluto_notebook(joinpath(@__DIR__, f))
end

function findindexstructure(file)
    index = []
    codelines = join(readlines(file), "\n")
    #
    basic_regex = r"md\"\"\"\nhashes (.*)\n"
    for i in 1:3
        string_regex = Regex(replace(basic_regex.pattern, "hashes" => repeat("\\#", i)))
        for t in eachmatch(string_regex, codelines)
            push!(index, (repeat("#", i) * " " * t.captures[1], t.offset))
        end
    end
    # 
    sort!(index, by=x -> getindex(x, 2))
    return getindex.(index, 1)
end

function markdown_string(notebook_name, notebooks_dir)
    notebook_relative_path = joinpath(notebooks_dir, notebook_name)
    notebook_path = joinpath(@__DIR__, notebook_name)
    structure = findindexstructure(notebook_path)
    title = structure[1][3:end]
    notebook_relative_path_html = replace(notebook_relative_path, ".jl" => ".html")
    "[$(notebook_name)]($(notebook_relative_path_html)) [[code]]($(notebook_relative_path)): $(title)"
end

notebook_list = markdown_string.(notebooks, Ref(notebooks_dir))
notebook_list_block = "## Notebooks\n\n" * prod(notebook_list) do l
    " - " * l * "\n"
end
readme_content = join(readlines(joinpath(@__DIR__, "..", "README.md")), "\n")
updated_readme_content = readme_content * "\n\n" * notebook_list_block


