
using CartesianGeneticProgramming
using LightGraphs
using MetaGraphs
using TikzGraphs
using TikzPictures
using LaTeXStrings


function to_graph(ind::CGPInd; active_outputs=trues(ind.n_out))
    actives = [n.active for n in ind.nodes]
    actives[1:ind.n_in] .= true
    vids = findall(actives)
    #pos = get_positions(c)
    mg = MetaDiGraph(SimpleDiGraph())
    add_vertices!(mg, length(vids)+sum(active_outputs))#c.nout)
    for i in 1:ind.n_in
        set_prop!(mg, i, :name, LaTeXString(string("\$in_{", i, "}\$")))
        set_prop!(mg, i, :type, 0)
    end
    for vi in (ind.n_in+1):length(vids)
        v = vids[vi]
        n = ind.nodes[v]
        f_name = split(split(repr(n.f), ".")[end], "_")[end]
        #if f_name == "const"
        #    set_prop!(mg, vi, :name, LaTeXString(@sprintf("%0.2f", n.p)))
        #else
        set_prop!(mg, vi, :name, LaTeXString(f_name))
        #end
        set_prop!(mg, vi, :function, n.f)
        set_prop!(mg, vi, :type, 2)
        #set_prop!(mg, vi, :param, n.p)
        cx = findfirst(x-> x==n.x,vids)
        cy = findfirst(x-> x==n.y,vids)
        if cx == cy
            add_edge!(mg, Edge(cx, vi))
            set_prop!(mg, cx, vi, :ci, 3)
        else
            add_edge!(mg, Edge(cx, vi))
            set_prop!(mg, cx, vi, :ci, 1)
            if typeof(cy)!=Nothing
                add_edge!(mg, Edge(cy, vi))
                set_prop!(mg, cy, vi, :ci, 2)
            end
        end
    end
    nid_count = 1
    for o in 1:ind.n_out
        if active_outputs[o]
            nid = length(vids)+nid_count
            set_prop!(mg, nid, :name, LaTeXString(string("\$out_{", o, "}\$")))
            set_prop!(mg, nid, :type, 1)
            oid = findfirst(x->x== ind.outputs[o],vids)
            add_edge!(mg, Edge(oid, nid))
            set_prop!(mg, nid, oid, :ci, 0)
            nid_count += 1
        end
    end
    mg
end

function chromo_draw(ind::CGPInd, file::String="graph.pdf"; active_outputs=trues(ind.n_out))
    mg = to_graph(ind, active_outputs=active_outputs)
    names = map(x->get_prop(mg, x, :name), 1:nv(mg))
    t = TikzGraphs.plot(mg.graph, names)
    TikzPictures.save(TikzPictures.PDF(file), t)
end
