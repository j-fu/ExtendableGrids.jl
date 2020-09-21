### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ d432ad64-f91f-11ea-2e48-4bc7472ac64c
begin
	using Pkg
	Pkg.develop("ExtendableGrids")
	Pkg.add("PyPlot")
	Pkg.add("PlutoUI")	
end

# ╔═╡ 07ebc9c4-f920-11ea-2545-111b1a4b25b3
using PyPlot, ExtendableGrids,PlutoUI

# ╔═╡ 4926eca8-fb6d-11ea-15a0-6fbe372acc98
with_terminal() do
	Pkg.status()
end

# ╔═╡ 3e177bd8-fbeb-11ea-0c9c-c9161a2d36d4
ENV["LC_NUMERIC"]="C"; 

# ╔═╡ f32d9f04-f923-11ea-3a4a-53cc3df5642c
html"""<hr>"""

# ╔═╡ fd27b44a-f923-11ea-2afb-d79f7e62e214
md"""
### Using the GridFactory
"""

# ╔═╡ 0e4a4a78-f927-11ea-1870-3b97ad00f8cd
triangleflags()

# ╔═╡ 511b26c6-f920-11ea-1228-51c3750f495c
begin
	global factory=GridFactory(flags=triangleflags(:domain))
	p1=point!(factory,(0,0))
	p2=point!(factory,(1,0))
	p3=point!(factory,(1,1))
	p4=point!(factory,(0,1))
	facet!(factory,p1,p2)
	facet!(factory,p2,p3)
	facet!(factory,p3,p4)
	facet!(factory,p1,p4)
	facet!(factory,p1,p3)
	cellregion!(factory,(0.1,0.5),region=1)
	cellregion!(factory,(0.9,0.5),region=2,volume=0.01)
	factory
end

# ╔═╡ 8f0bd5c0-f920-11ea-3b1c-db90fc95f990
ExtendableGrids.plot(factory,Plotter=PyPlot,resolution=(600,600))

# ╔═╡ 46b7244c-fa91-11ea-3d98-73191a884235


# ╔═╡ 54677c86-fa91-11ea-3518-3b8aba2c8488


# ╔═╡ Cell order:
# ╠═d432ad64-f91f-11ea-2e48-4bc7472ac64c
# ╠═07ebc9c4-f920-11ea-2545-111b1a4b25b3
# ╟─4926eca8-fb6d-11ea-15a0-6fbe372acc98
# ╠═3e177bd8-fbeb-11ea-0c9c-c9161a2d36d4
# ╟─f32d9f04-f923-11ea-3a4a-53cc3df5642c
# ╟─fd27b44a-f923-11ea-2afb-d79f7e62e214
# ╠═0e4a4a78-f927-11ea-1870-3b97ad00f8cd
# ╠═511b26c6-f920-11ea-1228-51c3750f495c
# ╠═8f0bd5c0-f920-11ea-3b1c-db90fc95f990
# ╟─46b7244c-fa91-11ea-3d98-73191a884235
# ╟─54677c86-fa91-11ea-3518-3b8aba2c8488
