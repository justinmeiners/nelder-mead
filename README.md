# Nelder-Mead

Fast C implementation of the [Nelder-Mead method](http://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)
multi-variable function optimization/minimization without derivatives.

This is a fork of [Matteo's work](https://github.com/matteotiziano/nelder-mead)  with the following improvements:

- improved performance and memory usage.
- optional parameters for configuring the initial simplex.
- improved interface especially focus on primitive types.
- single header style instead of two file.
- safe for C++ include.

## Include

Copy `nelder_mead.h` into your project. It is "single header style",
so whenever you want the forward definitions include it normally:

    #include "nm_optimizer.h"

Then it at least one `.c` file, include it again with a macro defined to get the implementation:

    #define NM_OPTIMIZER_IMPLEMENTATION
    #include "nm_optimizer.h`

## Example

See `main.c` for a simple example with the [Ackely function](http://www.sfu.ca/%7Essurjano/ackley.html).
The test can be compiled and tested from command line using the Makefile provided.

    make
    make test

The command line interface accepts a sequence of numbers as inputs:

    $ ./nm arg_1 arg_2 [...] arg_n

where `arg_i` is the i-th coordinate of n-dimensional initial position.
For example:

    $ ./nm -2.1 -3.04 4.5

    Iteration 0001     reflect         [ -2.10 -3.04 4.73 ]    11.00 
    Iteration 0002     expand          [ -1.96 -2.89 4.88 ]    10.46 
    Iteration 0003     reflect         [ -1.96 -2.89 4.88 ]    10.46 
    Iteration 0004     reflect         [ -1.90 -2.97 5.07 ]    10.45 
    Iteration 0005     contract out    [ -1.90 -2.97 5.07 ]    10.45 
    Iteration 0006     contract in     [ -1.92 -3.01 4.86 ]    10.44 
    Iteration 0007     reflect         [ -2.00 -2.96 4.89 ]    10.30 
    Iteration 0008     reflect         [ -2.00 -2.96 4.89 ]    10.30 
    Iteration 0009     contract in     [ -1.92 -3.00 4.99 ]    10.25 
    Iteration 0010     reflect         [ -1.92 -3.00 4.99 ]    10.25 
    Iteration 0011     contract out    [ -1.99 -2.94 4.97 ]    10.20 
    [...]
    Iteration 0040     contract in     [ -1.99 -2.98 4.97 ]    10.17 
    Initial point
    x = [ -2.10000000 -3.04000000 4.50000000 ], fx = 11.21198731 
    Solution
    x = [ -1.98940255 -2.98441227 4.97381887 ], fx = 10.16673927 

## Licence
MIT Licence. Copyright (c) 2017 Matteo Maggioni, 2021 Justin Meiners
