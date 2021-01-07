
function makeplot(func,picdir)
    PyPlot.clf()
    f=getfield(Main,Symbol(func))
    gridplot(f(), Plotter=PyPlot)
    fname=joinpath(picdir,func*".svg")
    PyPlot.savefig(fname)
    isfile(fname)
end

function makeplot2(func,picdir)
    PyPlot.clf()
    f=getfield(Main,Symbol(func))
    f(Plotter=PyPlot)
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
    makeplot("quadrilateral",picdir)

    makeplot2("plotting_grid3d",picdir)
    makeplot2("plotting_func3d",picdir)
    makeplot2("plotting_multiscene",picdir)
    
end
