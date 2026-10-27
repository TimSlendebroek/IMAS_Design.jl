# Included by ext/IMASGraphRecipesExt.jl (needs both Plots and GraphRecipes)
const Graphs = GraphRecipes.Graphs

@recipe function plot_circuit(circuit::IMAS.pf_active__circuit, supplies_labels::Vector{String}, coils_labels::Vector{String})
    title_str = circuit.name
    conn_matrix = circuit.connections
    num_supplies = length(supplies_labels)
    num_coils = length(coils_labels)

    conn_matrix = conn_matrix[:, 1:2:end] .+ conn_matrix[:, 2:2:end]
    num_nodes, num_elements = size(conn_matrix)
    @assert num_elements == num_supplies + num_coils "num_elements=$(num_elements) while (num_supplies + num_coils) * 2 = $((num_supplies + num_coils) * 2)"

    membership = [fill(2, num_nodes); fill(1, num_supplies); fill(3, num_coils)]

    names = []
    for k in 1:num_nodes
        push!(names, " $k ")
    end
    if isempty(supplies_labels)
        for k in 1:num_supplies
            push!(names, "S$k")
        end
    else
        append!(names, supplies_labels)
    end
    if isempty(coils_labels)
        for k in 1:num_coils
            push!(names, "C$k")
        end
    else
        append!(names, coils_labels)
    end
    names = collect(map(x -> replace(x, " " => "\n"), names))

    # hide unused supplies and coils
    index = (1:size(conn_matrix)[2])[dropdims(sum(conn_matrix; dims=1); dims=1).>0]
    conn_matrix = conn_matrix[:, index]
    membership = [membership[1:num_nodes]; membership[num_nodes.+index]]
    names = [names[1:num_nodes]; names[num_nodes.+index]]
    num_supplies = sum(membership .== 1)
    num_coils = sum(membership .== 3)

    # Pastel versions of red, green, and blue
    markercolor = [PlotUtils.Colors.colorant"#FFAAAA", PlotUtils.Colors.colorant"#AAFFAA", PlotUtils.Colors.colorant"#AAAAFF"]
    markercolor = markercolor[membership]
    node_weights = [1.0, 0.25, 1.0]
    node_weights = node_weights[membership]
    nodeshape = [:rect, :circle, :rect]
    nodeshape = nodeshape[membership]

    # Create an empty graph with nodes
    g = Graphs.SimpleDiGraph(num_nodes + num_supplies + num_coils)

    # Iterate through the matrix to add edges
    for i in 1:num_nodes
        connected_elements = findall(x -> x == 1, conn_matrix[i, :])
        for j in connected_elements
            Graphs.add_edge!(g, i, j + num_nodes)
        end
    end

    @series begin
        method --> :sfdp
        markercolor --> markercolor
        curves --> false
        arrows --> false
        names --> names
        nodeshape --> nodeshape
        node_weights --> node_weights
        titlefont --> font(12, "Arial", "bold")
        title --> title_str
        input_type := :sourcedestiny
        markerstrokewidth --> 0.0
        GraphRecipes.GraphPlot((g,))
    end
end
