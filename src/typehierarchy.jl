"""
$(TYPEDSIGNATURES)

Define children for types.
"""
AbstractTrees.children(T::Type)=InteractiveUtils.subtypes(T)

"""
$(TYPEDEF)

Apex type of all abstract types in this hierarchy.
"""
abstract type AbstractExtendableGridApexType end


"""
$(TYPEDSIGNATURES)

Print complete type hierachy for ExtendableGrids
"""
typehierarchy()=AbstractTrees.print_tree(AbstractExtendableGridApexType)


function leaftypes(TApex)
    function leaftypes!(leafs,t)
        st=subtypes(t)
        if length(st)==0
            push!(leafs,t)
        else
            for tsub in st
                leaftypes!(leafs,tsub)
            end
        end
        leafs
    end
    leaftypes!(Type[],TApex)
end
