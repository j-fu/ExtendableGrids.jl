
function makeplot(func,picdir; Plotter=nothing)
    f=getfield(Main,Symbol(func))
    p=gridplot(f(); Plotter)
    fname=joinpath(picdir,func*".svg")
    save(fname,p;Plotter)
    isfile(fname)
end

function makeplotx(func,picdir;Plotter=nothing)
    f=getfield(Main,Symbol(func))
    g,sg,sf=f()
    p=GridVisualizer(;Plotter, layout=(1,3),resolution=(600,200))
    gridplot!(p[1,1],g)
    gridplot!(p[1,2],sg)
    scalarplot!(p[1,3],sg,sf)
    fname=joinpath(picdir,func*".svg")
    save(fname,reveal(p);Plotter)
    isfile(fname)
end

function makeplots(picdir;Plotter=nothing)

    makeplot("interval_from_vector",picdir; Plotter)
    makeplot("interval_localref",picdir; Plotter)
    makeplot("interval_multiregion",picdir; Plotter)
    makeplot("interval_subgrid",picdir; Plotter)
    makeplot("rectangle",picdir; Plotter)
    makeplot("rectangle_localref",picdir; Plotter)
    makeplot("rectangle_multiregion",picdir; Plotter)
    makeplot("rectangle_subgrid",picdir; Plotter)
    makeplot("rect2d_bregion_function",picdir; Plotter)

    makeplot("quadrilateral",picdir; Plotter)
    makeplot("cross3d",picdir; Plotter)
    makeplotx("sorted_subgrid",picdir; Plotter)        
    
end
