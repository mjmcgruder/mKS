# mKS

A Kuramoto-Sivashinsky equation solver using the discontinuous Galerkin finite element method.

## Building

The build system is a single Make file with supplementary build options specified in a machine-specific local.mk.

To start, install BLAS and LAPACK. For instance, on Ubuntu the following is sufficient.

```
sudo apt install liblapack-dev
```

Next, make a local.mk in the project root directory. It's easiest to copy the example if g++ is avilable in your PATH.

```
cp examples/local_example.mk local.mk
```

Finally

```
make
```
