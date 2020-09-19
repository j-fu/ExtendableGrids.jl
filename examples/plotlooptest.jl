using ExtendableGrids
import Makie, PyPlot,Plots,VTKView


X=collect(0.0:0.05:1)

grid1d=simplexgrid(X)
grid2d=simplexgrid(X,X)

func1d=map(x->sin(5*x),grid1d)
func2d=map((x,y)->sin(5*x)*sin(5*y),grid2d)


function gridloop(Plotter;kwargs...)
    ctx=PlotterContext(Plotter,kwargs...)
    @time for i=1:20
        X=collect(0.0:0.05:1)
        Y=collect(0.0:0.05:1+0.1*i)
        grid=simplexgrid(X,Y)
        plot!(ctx,grid)
    end
end


function funcloop(Plotter;kwargs...)
    ctx=PlotterContext(Plotter,kwargs...)
    X=collect(0.0:0.05:1)
    grid=simplexgrid(X,X)
    @time for i=1:20
        func2d=map((x,y)->sin(0.5*i*x)*sin(0.5*i*y),grid2d)
        plot!(ctx,grid,func2d)
    end
end


