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
typehierarchy()=AbstractTrees.print_tree(AbstractExtendableGridApexType,5,indicate_truncation=false)

