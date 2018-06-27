vectormultibang
===============

This repository contains Matlab codes accompanying the paper [Vector-valued multibang control of differential equations](https://doi.org/10.1137/16M1104998) ([arXiv preprint](http://arxiv.org/abs/1611.07853)) by [Christian Clason](http://udue.de/clason), Carla Tameling, and Benedikt Wirth.

Contents
--------

##### directory `bloch` 
contains the test scripts for the example concerning optimal control of the Bloch equation using discrete control vectors (run `test_bloch.m`)

##### directory `elasticity` 
contains the test scripts for the example concerning optimal control of the equations of linearized elasticity (run `test_elasticity.m`)

If you find this approach useful, you can cite the paper as

    @article{vectormultibang,
        author = {Clason, Christian and Tameling, Carla and Wirth, Benedikt},
        title = {Vector-valued multibang control of differential equations},
        journal = {SIAM Journal on Control and Optimization},
        volume = {56},
        number = {3},
        pages = {2295--2326},
        doi = {10.1137/16M1104998},
        year = {2018}
    }
