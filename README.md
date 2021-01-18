# speig
A Matlab solver for **s**ym**p**lectic **eig**envalue problem via trace minimization and Riemannian optimization

## Problem
This solver aims to compute the symplectic eigenvalues and eigenvectors of a positive-definite matrix M.

The major step is based on the Riemannian optimization over the symplectic Stiefel manifold,

> min trace(X'AX), s.t.  X' J2n X = J2p.

where X is a 2n-by-2p matrix, J2n = [0 In; -In 0], and In is the n-by-n identity matrix.

## References
[Bin Gao](https://www.gaobin.cc/), [Nguyen Thanh Son](https://sites.google.com/view/ntson), [P.-A. Absil](https://sites.uclouvain.be/absil/), [Tatjana Stykel](https://www.uni-augsburg.de/en/fakultaet/mntf/math/prof/numa/team/tatjana-stykel/)
1. [Riemannian optimization on the symplectic Stiefel manifold](https://arxiv.org/abs/2006.15226)
2. [Symplectic eigenvalue problem via trace minimization and Riemannian optimization](https://arxiv.org/abs/2101.02618)

## Authors
+ [Nguyen Thanh Son](https://sites.google.com/view/ntson) (Thai Nguyen University of Sciences, Vietnam)
+ [Bin Gao](https://www.gaobin.cc/) (UCLouvain, Belgium)

## Copyright
Copyright (C) 2020, Nguyen Thanh Son, P.-A. Absil, Bin Gao, Tatjana Stykel

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/)
