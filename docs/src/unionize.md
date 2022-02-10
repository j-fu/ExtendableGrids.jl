# unionize
## Unionize collections

When working e.g with agent based models or finite elements with varying element geometries, a common pattern is the occurence of collections of objects of different types on which one wants to perform certain actions depending on their type.

We can observe at least three patterns related to this type of situation:

- __Collection of objects:__ This is e.g. a vector of instantations of structs of different types, and the most intuitive way of using this pattern. Calling methods of some function on a member of this collection leads to [dynamic dispatch](https://discourse.julialang.org/t/dynamic-dispatch/6963/2), which has certain runtime costs due to necessary allocations.
- __Collection of types:__  Julia allows to use types as variables. These can be stored  in a collection as well, and it is possible to  dispatch on the the type. No need to have an instantiation of the type, and one also can work with abstract type here.
- __Collection of functions:__ Instead of objects or types, one also can store functions in a collection. Accessing a function as a member of a collection once again will lead to dynamic dispatch.

[Union splitting](https://julialang.org/blog/2018/08/union-splitting/) allows to avoid dynamic dispatch and can save significant runtime costs. Below, we try to explain, how this works.

Here, we provide an example for a collection of objects.
```jldoctest
N=10_000

struct S01 end
struct S02 end
struct S03 end
struct S04 end
struct S05 end

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

function sumup_f_manual(collection)
	s=0.0
	for i=1:length(collection)
		c=collection[i]
		if isa(c,S01)
			s+=f(c,1)
		elseif isa(c,S02)
			s+=f(c,1)
		elseif isa(c,S03)
			s+=f(c,1)
		elseif isa(c,S04)
			s+=f(c,1)
		elseif isa(c,S05)
			s+=f(c,1)
		end
	end
	s
end


any_collection=[rand(allS)() for i=1:N]

unionS_collection=UnionS[s for s ∈ any_collection]

sumup_f(any_collection) # hide
sumup_f_manual(any_collection) #hide
sumup_f(unionS_collection) # hide

a_any=@allocated sumup_f(any_collection)
a_manual=@allocated sumup_f_manual(any_collection)
a_union=@allocated sumup_f(unionS_collection)
a_any, a_manual, a_union
# output
(160000, 16, 16)
```

As the benchmark shows, when defining the collection in the default way (resulting in a `Vector{Any}`, each access of an element is linked to an allocation with significant runtime overhead.

With "manual dispatch",  each time when c is accessed as a function parameter, due to the test via `isa`, the compiler knows the type of `c` and can choose the proper method of `f` at compile time,  avoiding the allocations.
While it is possible to generate the manual dispatch code with macros, another remedy of this situation appears to be more acessible.

"Unionizing" the collection means that one pins its element type to to the union of possible types of entries.
The compiler then  knows that the number of possible types of the elements of the collection is finite -- constrained by the list of types in the union. Consequently, it can automatically create code similar to the manual dispatch statement above.

This pattern works for all the cases mentioned above:
- vectors of instances of different types (as discussed in this thread), via Vector{Union{T1,T2}}
- for vectors of concrete types - one would use Vector{Union{Type{T1},Type{T2}}}
- for vectors of abstract types, same as for concrete types
- for vectors of functions, also with multiple methods via Vector{Union{typeof(f1),typeof(f2)}}

