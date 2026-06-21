
## A.63 library(ugraphs): Graph manipulation library

author  
\- R.A.O'Keefe  
- Vitor Santos Costa  
- Jan Wielemaker

license  
BSD-2 or Artistic 2.0

The S-representation of a graph is a list of (vertex-neighbours) pairs, where the pairs are in standard order (as produced by keysort) and the neighbours of each vertex are also in standard order (as produced by sort). This form is convenient for many calculations.

A new UGraph from raw data can be created using [vertices_edges_to_ugraph/3](ugraphs.html#vertices_edges_to_ugraph/3).

Adapted to support some of the functionality of the SICStus ugraphs library by Vitor Santos Costa.

Ported from YAP 5.0.1 to SWI-Prolog by Jan Wielemaker.

**vertices**(`+Graph, -Vertices`)  
Unify `Vertices` with all vertices appearing in `Graph`. Example:

``` code
?- vertices([1-[3,5],2-[4],3-[],4-[5],5-[]], L).
L = [1, 2, 3, 4, 5]
```

\[det\]**vertices_edges_to_ugraph**(`+Vertices:list, +Edges:pairs, -UGraph`)  
Create a `UGraph` from `Vertices` and `Edges`. `UGraph` must unify with the corresponding S-representation. Note that vertices that do not appear in any of the `Edges` appear in `UGraph` as `Vertice-[]`. The set of vertices in `UGraph` is the union of `Vertices` and all vertices that appear in the `Edges` pairs.

``` code
?- vertices_edges_to_ugraph([],[1-3,2-4,4-5,1-5], L).
L = [1-[3,5], 2-[4], 3-[], 4-[5], 5-[]]
```

In this case all vertices are defined implicitly. The next example shows three unconnected vertices:

``` code
?- vertices_edges_to_ugraph([1,2,6,7,8],[1-3,2-4,4-5,1-5], L).
L = [1-[3,5], 2-[4], 3-[], 4-[5], 5-[], 6-[], 7-[], 8-[]]
```

**add_vertices**(`+Graph, +Vertices, -NewGraph`)  
Unify `NewGraph` with a new graph obtained by adding the list of `Vertices` to `Graph`. Example:

``` code
?- add_vertices([1-[3,5],2-[]], [0,1,2,9], NG).
NG = [0-[], 1-[3,5], 2-[], 9-[]]
```

\[det\]**del_vertices**(`+Graph, +Vertices, -NewGraph`)  
Unify `NewGraph` with a new graph obtained by deleting the list of `Vertices` and all the edges that start from or go to a vertex in `Vertices` to the `Graph`. Example:

``` code
?- del_vertices([1-[3,5],2-[4],3-[],4-[5],5-[],6-[],7-[2,6],8-[]],
                [2,1],
                NL).
NL = [3-[],4-[5],5-[],6-[],7-[6],8-[]]
```

Compatibility  
Upto 5.6.48 the argument order was (+`Vertices`, +`Graph`, -`NewGraph`). Both YAP and SWI-Prolog have changed the argument order for compatibility with recent SICStus as well as consistency with [del_edges/3](ugraphs.html#del_edges/3).

**add_edges**(`+Graph, +Edges, -NewGraph`)  
Unify `NewGraph` with a new graph obtained by adding the list of `Edges` to `Graph`. Example:

``` code
?- add_edges([1-[3,5],2-[4],3-[],4-[5],
              5-[],6-[],7-[],8-[]],
             [1-6,2-3,3-2,5-7,3-2,4-5],
             NL).
NL = [1-[3,5,6], 2-[3,4], 3-[2], 4-[5],
      5-[7], 6-[], 7-[], 8-[]]
```

**ugraph_union**(`+Graph1, +Graph2, -NewGraph`)  
`NewGraph` is the union of `Graph1` and `Graph2`. Example:

``` code
?- ugraph_union([1-[2],2-[3]],[2-[4],3-[1,2,4]],L).
L = [1-[2], 2-[3,4], 3-[1,2,4]]
```

**del_edges**(`+Graph, +Edges, -NewGraph`)  
Unify `NewGraph` with a new graph obtained by removing the list of `Edges` from `Graph`. Notice that no vertices are deleted. Example:

``` code
?- del_edges([1-[3,5],2-[4],3-[],4-[5],5-[],6-[],7-[],8-[]],
             [1-6,2-3,3-2,5-7,3-2,4-5,1-3],
             NL).
NL = [1-[5],2-[4],3-[],4-[],5-[],6-[],7-[],8-[]]
```

**edges**(`+Graph, -Edges`)  
Unify `Edges` with all edges appearing in `Graph`. Example:

``` code
?- edges([1-[3,5],2-[4],3-[],4-[5],5-[]], L).
L = [1-3, 1-5, 2-4, 4-5]
```

**transitive_closure**(`+Graph, -Closure`)  
Generate the graph `Closure` as the transitive closure of `Graph`. Example:

``` code
?- transitive_closure([1-[2,3],2-[4,5],4-[6]],L).
L = [1-[2,3,4,5,6], 2-[4,5,6], 4-[6]]
```

\[det\]**transpose_ugraph**(`Graph, NewGraph`)  
Unify `NewGraph` with a new graph obtained from `Graph` by replacing all edges of the form V1-V2 by edges of the form V2-V1. The cost is O(`|`V`|`\*log(`|`V`|`)). Notice that an undirected graph is its own transpose. Example:

``` code
?- transpose([1-[3,5],2-[4],3-[],4-[5],
              5-[],6-[],7-[],8-[]], NL).
NL = [1-[],2-[],3-[1],4-[2],5-[1,4],6-[],7-[],8-[]]
```

Compatibility  
This predicate used to be known as [transpose/2](clpfd.html#transpose/2). Following SICStus 4, we reserve [transpose/2](clpfd.html#transpose/2) for matrix transposition and renamed ugraph transposition to [transpose_ugraph/2](ugraphs.html#transpose_ugraph/2).

**compose**(`+LeftGraph, +RightGraph, -NewGraph`)  
Compose `NewGraph` by connecting the *drains* of `LeftGraph` to the *sources* of `RightGraph`. Example:

``` code
?- compose([1-[2],2-[3]],[2-[4],3-[1,2,4]],L).
L = [1-[4], 2-[1,2,4], 3-[]]
```

\[semidet\]**ugraph_layers**(`Graph, -Layers`)  
\[semidet\]**top_sort**(`+Graph, -Sorted`)  
Sort vertices topologically. `Layers` is a list of lists of vertices where there are no edges from a layer to an earlier layer. The predicate [top_sort/2](ugraphs.html#top_sort/2) flattens the layers using [append/2](lists.html#append/2).

These predicates fail if `Graph` is cyclic. If `Graph` is not connected, the sub-graphs are individually sorted, where the root of each subgraph is in the first layer, the nodes connected to the roots in the second, etc.

``` code
?- top_sort([1-[2], 2-[3], 3-[]], L).
L = [1, 2, 3]
```

Compatibility  
\- The original version of this library provided top_sort/3 as a *difference list* version of [top_sort/2](ugraphs.html#top_sort/2). We removed this because the argument order was non-standard. Fixing causes hard to debug compatibility issues while we expect top_sort/3 was rarely used. A backward compatible top_sort/3 can be defined as

``` code
top_sort(Graph, Tail, Sorted) :-
    top_sort(Graph, Sorted0),
    append(Sorted0, Tail, Sorted).
```

The original version returned all vertices in a *layer* in reverse order. The current one returns them in standard order of terms, i.e., each layer is an *ordered set*.  
- [ugraph_layers/2](ugraphs.html#ugraph_layers/2) is a SWI-Prolog specific addition to this library.

\[det\]**neighbors**(`+Vertex, +Graph, -Neigbours`)  
\[det\]**neighbours**(`+Vertex, +Graph, -Neigbours`)  
`Neigbours` is a sorted list of the neighbours of `Vertex` in `Graph`. Example:

``` code
?- neighbours(4,[1-[3,5],2-[4],3-[],
                 4-[1,2,7,5],5-[],6-[],7-[],8-[]], NL).
NL = [1,2,7,5]
```

\[det\]**connect_ugraph**(`+UGraphIn, -Start, -UGraphOut`)  
Adds `Start` as an additional vertex that is connected to all vertices in `UGraphIn`. This can be used to create an topological sort for a not connected graph. `Start` is before any vertex in `UGraphIn` in the standard order of terms. No vertex in `UGraphIn` can be a variable.

Can be used to order a not-connected graph as follows:

``` code
top_sort_unconnected(Graph, Vertices) :-
    (   top_sort(Graph, Vertices)
    ->  true
    ;   connect_ugraph(Graph, Start, Connected),
        top_sort(Connected, Ordered0),
        Ordered0 = [Start|Vertices]
    ).
```

**complement**(`+UGraphIn, -UGraphOut`)  
`UGraphOut` is a ugraph with an edge between all vertices that are *not* connected in `UGraphIn` and all edges from `UGraphIn` removed. Example:

``` code
?- complement([1-[3,5],2-[4],3-[],
               4-[1,2,7,5],5-[],6-[],7-[],8-[]], NL).
NL = [1-[2,4,6,7,8],2-[1,3,5,6,7,8],3-[1,2,4,5,6,7,8],
      4-[3,5,6,8],5-[1,2,3,4,6,7,8],6-[1,2,3,4,5,7,8],
      7-[1,2,3,4,5,6,8],8-[1,2,3,4,5,6,7]]
```

To be done  
Simple two-step algorithm. You could be smarter, I suppose.

**reachable**(`+Vertex, +UGraph, -Vertices`)  
True when `Vertices` is an ordered set of vertices reachable in `UGraph`, including `Vertex`. Example:

``` code
?- reachable(1,[1-[3,5],2-[4],3-[],4-[5],5-[]],V).
V = [1, 3, 5]
```
