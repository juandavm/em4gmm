# em4gmm

Extremely fast C implementation of the clustering Expectation Maximization (EM) algorithm for estimating Gaussian Mixture Models (GMM).

## Introduction

In statistics, the expectation–maximization (EM) algorithm is an iterative method for finding maximum likelihood or maximum a posteriori (MAP) estimates of parameters in statistical models, where the model depends on unobserved latent variables.

For this reason EM is frequently used for data clustering, verification and identification of the speaker (biometric tasks), author profiling based on his documents, automatic document categorization, and many more applications.

## Download

-   [Download the latest release of the source code on a zip file](https://github.com/juandavm/em4gmm/zipball/master).
-   Clone the repository: `git clone git://github.com/juandavm/em4gmm.git`

## Installing

In order to compile this program, you need first compile and install the [zlib library](http://www.zlib.net) from their website, or using your preferred software distribution channels (aptitude, yum, macports, etc) in order to install it (with the `dev` packages). In some systems this library is installed by default.

On MacOs and Linux distributions you can simple use the `make` command on the system shell to compile it, and then `sudo make install` to install it on yout system (by default on /usr/bin).

## Usage

-   [Download the last version of the user manual in pdf format](https://github.com/juandavm/em4gmm/blob/master/doc/UserManual.pdf?raw=true).

## Used by...

(does not imply any endorsement)

-   Steinbach, B. (Ed.). (2016). Problems and new solutions in the boolean domain. Cambridge Scholars Publishing.
-   Vandyke, D. J. (2014). Glottal Waveforms for Speaker Inference & A Regression Score Post-Processing Method Applicable to General Classification Problems (Doctoral dissertation, University of Canberra).
-   Schmidt, J., Blazek, R., & Fišer, P. (2014). On Probability Density Distribution of Randomized Algorithms Performance. In Proc. of 11th Int. Workshop on Boolean Problems (IWSBP). Freiberg (Germany) (pp. 67-74).
-   Schmidt, J., Fier, P., & Balcárek, J. (2014, August). On robustness of eda tools. In 2014 17th Euromicro Conference on Digital System Design (pp. 427-434). IEEE.

## Issues and Bugs

Do you have a bug or a feature request? Do not worry, [open a new issue](https://github.com/juandavm/em4gmm/issues). But please, before opening any new issue, search on existing the yours in order to avoid duplicates. And thanks you for your contribution!

## Original repository

-   [https://github.com/juandavm/em4gmm/](https://github.com/juandavm/em4gmm/)

## License

Expectation Maximization for Gaussian Mixture Models. Copyright (C) 2012-2022 Juan Valor.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
