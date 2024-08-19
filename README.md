# PointSearchNumberFields

This repository contains [Magma](http://magma.maths.usyd.edu.au/magma/) code developed as part of a [URECA](https://ureca.wfu.edu/) project at Wake Forest by Flora Yi together with her advisors Jacob Mayle and Jeremy Rouse.

Given a projective curve `C` (defined over a number field *K*) and a specified `bound`, the main function `PointSearchNF` searches for *K*-rational points on `C` of height up to `bound`. Not every point on `C` of height up to `bound` is guaranteed to be found and points of higher height may be returned as well. The notion of height we use is defined as follows.

> Let $\mathcal{O}_K$ denote the ring of integers of $K$ and fix a integral basis $B_K = [b_1,...,b_d]$ of $\mathcal{O}_K$. Consider all representations of $P$ as $P = (x_0:x_1:\ldots:x_n)$ with $x_0,\ldots,x_n \in \mathcal{O}_K$. The *height* of $P$ is defined to be the minimum value (over all such representations) of $`\max_{i,j} |c_{i,j}|`$ where $`x_0 = c_{0,1} b_1 + ... + c_{0,d} b_d,\ldots, x_n = c_{n,1} b_1 + ... + c_{n,d} b_d`$, with each $`c_{i,j}`$ an integer and where $|·|$ denotes the absolute value in the usual sense. Note that the height of $P$ depends on the choice of basis $B_K$ for $\mathcal{O}_K$. The function `HeightByBasis`, available in `PointSearchNF.m`, computes the height of a point in $\mathbb{P}^n(K)$ in this sense.


The function `PointSearchNF` follows a lifting-based approach that involves searching integer lattices, building on the algorithm described in C.L. Turner's Ph.D. thesis (see also the notes of M. Watkins). A key improvement in our implementation is in choosing the ideal to work with. Our code searches for an ideal that will (heuristically) optimize the run time, resulting in significant runtime reductions in many cases compared to a more naive ideal choice.

Please contact us with any questions, comments, or suggestions.

References
1. Turner, Charlotte L. Lattice methods for finding rational points on varieties over number fields. Thesis (Ph.D.)–University of Warwick (United Kingdom) 2013.
2. Watkins, Mark. Searching for points with the Elkies ANTS-IV algorithm. Unpublished.
