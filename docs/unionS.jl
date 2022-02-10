using BenchmarkTools
N=10_000
begin 
	struct S01 end
	struct S02 end
	struct S03 end
	struct S04 end
	struct S05 end
end

allS=(S01,S02,S03,S04,S05)

const UnionS=Union{allS...}

f(::S01,x)=1x
f(::S02,x)=2x
f(::S03,x)=3x
f(::S04,x)=4x
f(::S05,x)=5x

function sumup_f(collection)
	s=0.0
	for i=1:length(collection)
		s+=f(collection[i],1)
	end
	s
end

any_collection=[rand(allS)() for i=1:N]

unionS_collection=UnionS[s for s ∈ any_collection]

@btime sumup_f($any_collection)
@btime sumup_f($unionS_collection)
