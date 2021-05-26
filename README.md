vectormultibang
===============

This repository contains Matlab codes accompanying the paper [Convex relaxation of discrete vector-valued optimization problems](https://doi.org/10.1137/16M1104998) ([arXiv preprint](http://arxiv.org/abs/1611.07853)) by [Christian Clason](https://homepage.uni-graz.at/c.clason), Carla Tameling, and [Benedikt Wirth](https://www.uni-muenster.de/AMM/num/wirth/people/Wirth/index.html).

Contents
--------

##### directory `bloch`
contains the test scripts for the example concerning optimal control of the Bloch equation using discrete control vectors (run `test_bloch.m`)

##### directory `elasticity`
contains the test scripts for the example concerning optimal control of the equations of linearized elasticity (run `test_elasticity.m`)

##### directory `mbtransport`
contains the test scripts for the example concerning multimaterial branched transport (run `test_branchedTransport.m`)

If you find this approach useful, you can cite the paper as

    @article{vectormultibang,
        author = {Clason, Christian and Tameling, Carla and Wirth, Benedikt},
        title = {Convex relaxation of discrete vector-valued optimization problems},
        journal = {SIAM Review},
        volume = {63},
        number = {4},
        year = {2021}
    }
