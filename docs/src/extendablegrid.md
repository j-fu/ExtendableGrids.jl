# Extendable grid

Grid data structure based on the idea of a type stable dictionary:
## TDict interface

Extendable  containers with type stable content access and lazy content creation via the Julia type system.

### Problem to be addressed

In certain contexts it is desirable to use containers with core components
which are user extendable and allow for type stable component acces. Moreover,
some components are necessary on demand only, so they should be created lazily.
Furthermore, there should be a kind of safety protocol which prevents errors
from typos in component names etc.

Julia default data structures do not provide these properties.

#### `struct` 
  - Julia structs with proper field type annotations guarantee type stability
  - Julia structs are not extendable, fields and their types are fixed upon definition
  - If we don't fix types of struct fields they become Any and a source 
    for type instability
  - The situation could be fixed if `getfield` could be overloaded but it cant't

#### `Dict`
  - Plain Dicts with flexible value types are a source of type instability
  - Typical use of Dicts with strings as keys needs a meta protocol to handle
    semantics of keys which at the end probably hinges on string comparison which
    will make things slow
  - Same for the implementation of a lazy evaluation protocol

### Proposed solution:

Harness the power of the Julia type system: 
- Use a struct containing a  Dict with DataType as keys. Every key is a type.
- Use the type system to dispatch between  `getindex`/`setindex!` methods for keys
- Use type hierarchies to manage different different value classes
- Extension requires declaring new types, keys can be only existing types almost removing
  typos as sources for errors
- Lazy extension is managed bye an  `instantiate` method called by `getindex` if necessary
- Component access is made type stable by type dependent `getindex` methods
- Component insertion is made safe by having  `setindex!`  calling a `veryform` method

#### Pros
See above ...

#### Cons
- Implemented using a Dict, so access is inherently slower than access to a component
  of a struct. Therefore it is not well suited for inner loops.




## Grid API

```@autodocs
Modules = [XGrid]
Pages = ["extendablegrid.jl"]
```

