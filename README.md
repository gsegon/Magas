![CMake-Build](https://github.com/gsegon/solver/actions/workflows/cmake.yml/badge.svg?event=push)


Magas
=====

Magas is an open source magnetostatic finite element solver. 


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


Installation from package managers:
-----------------------------------

Ubuntu/Debian

Add magas and dependecy repositories:

    sudo add-apt-repository ppa:ginggs/deal.ii-9.4.0-backports
    sudo add-apt-repository ppa:gsegon/magas
    sudo apt update

Install magas:

    sudo apt-get install magas

    

Run an example:
--------------

Before running an example:

    $ sudo chown -R $USER ~/magas
    

Run an example:

    $ cd ~/magas/examples/motoric
    $ magas motoric.json

To add a screenshot, create an `assets/images` folder in your repository and upload your screenshot to it. Then, using the relative filepath, add it to your README using the following syntax:

    ```md
    ![alt text](assets/images/screenshot.png)
    ```

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



## Credits

List your collaborators, if any, with links to their GitHub profiles.

If you used any third-party assets that require attribution, list the creators with links to their primary web presence in this section.

If you followed tutorials, include links to those here as well.

## License

The last section of a high-quality README file is the license. This lets other developers know what they can and cannot do with your project. If you need help choosing a license, refer to [https://choosealicense.com/](https://choosealicense.com/).

---

🏆 The previous sections are the bare minimum, and your project will ultimately determine the content of this document. You might also want to consider adding the following sections.

## Badges

![badmath](https://img.shields.io/github/languages/top/lernantino/badmath)

Badges aren't necessary, per se, but they demonstrate street cred. Badges let other developers know that you know what you're doing. Check out the badges hosted by [shields.io](https://shields.io/). You may not understand what they all represent now, but you will in time.

## Features

If your project has a lot of features, list them here.

## How to Contribute

If you created an application or package and would like other developers to contribute it, you can include guidelines for how to do so. The [Contributor Covenant](https://www.contributor-covenant.org/) is an industry standard, but you can always write your own if you'd prefer.

## Tests

Go the extra mile and write tests for your application. Then provide examples on how to run them here.
