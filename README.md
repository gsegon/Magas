![CMake-Build](https://github.com/gsegon/solver/actions/workflows/cmake.yml/badge.svg?event=push)


Magas
=====

Magas is an open source magnetostatic finite element solver featuring:
- Linear and Non-linear (experimental) 2D solver.
- Mesh support: Gmsh.
- Material definition as a constant permitivty or a CSV text file describing BH curve.
- VTU export for visualization with Paraview.
- JSON export for scalar results.
- User defined post processing Scalar and Field functions.
- Torque, force, flux linkage, energy,  Magnetic flux density (B), Magnetic field strength (H) 
- Force calculation based on Eggshell and Lorentz methods.
- Flux linkage calulcation.
- Periodic and anti-periodic boundaries.
- Virtual mesh rotation (experimental).


Installation from sources:
--------------------------
If building from sources, make sure all the dependencies are installed. See the dependecy list for more details.


Let's say you've unpacked the .tar.gz file into a directory /path/to/magas/sources. 
Then configure, compile, and install Magas with:

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/path/where/magas/should/be/installed/to /path/to/magas/sources
    sudo make install    (alternatively $ make -j<N> install)
    make test

To build from the repository, execute the following commands first:

    git clone https://github.com/gsegon/Magas.git
    cd Magas

Then continue as before.


Installation from package managers (Ubuntu/Debian):
-----------------------------------



Add magas and dependecy repositories:

    sudo add-apt-repository ppa:ginggs/deal.ii-9.4.0-backports
    sudo add-apt-repository ppa:gsegon/magas
    sudo apt update

Install magas:

    sudo apt-get install magas

    

Run an example:
--------------

Before running an example:

    sudo chown -R $USER ~/magas
    

Run an example:

    cd ~/magas/examples/motoric
    magas motoric.json

After a successful run two new files should be created in the current directory.

Results of user defined postprocessing section found in motoric.json
    
    results-motoric.json 

VTU file containing visualization of the mesh, solution and results of user defined postprocessor expessions found in motoric.json. Paraview may be used for visualization.

    motoric.vtu
    

Dependency list:
----------------

The only dependecy that needs to be installed spearately from magas is [deal.ii](https://www.dealii.org/). Instruction on how to install deal.ii can be found on the it's website. 
Other dependencies are included with Magas source code and can be found in _external_ directory.

Dependency  | Version | Bundled 
------------- | ------------- | ---------
[deal.ii](https://github.com/dealii/dealii)  | >=9.4.0 | No
[json](https://github.com/nlohmann/json)  |  | Yes
[cxxopts](https://github.com/jarro2783/cxxopts)  |  | Yes
[exprtk](https://github.com/ArashPartow/exprtk)  |  | Yes
[vtu11](https://github.com/phmkopp/vtu11)  |  | Yes



## License

See LICENSE file.

