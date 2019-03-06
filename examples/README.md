# Jupyter notebook examples

All the notebooks were create using [Jupyter notebooks](https://github.com/jupyter/jupyter) with different kernels:
* C++: [xeus-cling](https://github.com/QuantStack/xeus-cling)
  Install with conda. Note that you should use the same install steps used in their .travis.yml file:
    ```bash
    conda install xeus=0.18.1 cling=0.5.0 clang_variant=*=cling_6.14.06 clangdev=5.0 llvmdev=5 nlohmann_json=3.4.0 cppzmq=4.3.0 xtl=0.5.2 pugixml cxxopts=2.1.1 -c conda-forge
    conda install cmake -c conda-forge
    ```
  And then compile from source with cmake as shown in the project's README.md.  
* Perl: [iperl](https://metacpan.org/pod/Devel::IPerl)
  Install with perl's [cpan](https://metacpan.org/).
* Python: [ipython](https://github.com/ipython/ipython)
  Install with conda. It should be part of the default environment in anaconda.
* R: [IRKernel](https://irkernel.github.io)
  Follow the instructions [here](https://irkernel.github.io/installation/) to install

For all kernels just launch the Jupyter notebook server with `jupyter notebook` on the command line and choose the appropriate kernel.

For a complete list of available Jupyter notebook kernels go to https://github.com/jupyter/jupyter/wiki/Jupyter-kernels.
