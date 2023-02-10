
function makeplot(func,picdir)
    PyPlot.clf()
    f=getfield(Main,Symbol(func))
    gridplot(f(), Plotter=PyPlot)
    fname=joinpath(picdir,func*".svg")
    PyPlot.savefig(fname)
    isfile(fname)
end

function makeplotx(func,picdir)
    PyPlot.clf()
    f=getfield(Main,Symbol(func))
    g,sg,sf=f()
    vis=GridVisualizer(Plotter=PyPlot, layout=(1,3),resolution=(600,200))
    gridplot!(vis[1,1],g)
    gridplot!(vis[1,2],sg)
    scalarplot!(vis[1,3],sg,sf)
    fname=joinpath(picdir,func*".svg")
    PyPlot.savefig(fname)
    isfile(fname)
end

function makeplots(picdir)

    makeplot("interval_from_vector",picdir)
    makeplot("interval_localref",picdir)
    makeplot("interval_multiregion",picdir)
    makeplot("interval_subgrid",picdir)
    makeplot("rectangle",picdir)
    makeplot("rectangle_localref",picdir)
    makeplot("rectangle_multiregion",picdir)
    makeplot("rectangle_subgrid",picdir)
    makeplot("rect2d_bregion_function",picdir)

    makeplot("quadrilateral",picdir)
    makeplot("cross3d",picdir)
    
end
