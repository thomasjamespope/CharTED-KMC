<table align="center">
<tr><td align="center" width="10000">

# <strong> CharTED-KMC </strong>

<p>
    <a href="https://orcid.org/0000-0001-7552-9812">Dr. Thomas Pope </a>
    <br>
    <a href="https://ncl.ac.uk/nes/people/profile/tompenfold.html">Prof. Thomas Penfold </a>
    <br>
    <a href="https://orcid.org/0000-0002-8409-6702">Dr. Yvelin Giret </a>
</p>

<p>
    <a href="http://penfoldgroup.co.uk">Penfold Group </a> @ <a href="https://ncl.ac.uk">Newcastle University </a>
</p>

<p>
    <a href="#setup">Setup</a> • <a href="#getting">Quickstart</a> • <a href="#publications">Publications</a>
</p>

</td></tr></table>

#
CharTED-KMC is a Kinetic Monte-Carlo code used to model charge transport in organic semiconductors.

​The code has the following capabilities:

    -- Transport of electrons, holes and excitons.
    -- Interplay between multiple excited states on each site
    -- Miller-Abrahams and Marcus rates
    -- Order parameters to control relative orientation and order within thin film.

</a>

## SETUP

The quickest way to get started with CharTED-KMC is to clone this repository:

<!---
```
git clone https://github.com/thomasjamespope/CharTED-KMC/tree/main 
```
--->

This contains all the source files and an example input file - as well as the offline manual.

CharTED-KMC is entirely stand-alone and comes with its own Makefile. The compile, simply go to the code directory and type:

<!---
```
make all
```
--->

Now you're good to go!

## GETTING STARTING 
|---:|:---:|:---:|:---|
| DEBUG | 1 boolean | .false. | Sets the seed for the random number generator to zero for reproducable calculations |

## LICENSE

This project is licensed under the GPL-3.0 License - see the LICENSE.md file for details.

## PUBLICATIONS

[T. Pope](https://orcid.org/0000-0001-7552-9812), [Y. Giret](https://orcid.org/0000-0002-8409-6702), M. Fsadni, [P. Docampo](https://orcid.org/0000-0001-6164-4748), [C. Groves](https://orcid.org/0000-0003-2402-1618), & [TJ. Penfold](https://orcid.org/0000-0003-4490-5672) (2023). *[Modelling the effect of dipole ordering on charge-carrier mobility in organic semiconductors](https://doi.org/10.1016/j.orgel.2023.106760)*, *Organic Electronics*, **115**, 106760.
<table align="center">
<tr><td align="center" width="10000">

# <strong> CharTED-KMC </strong>

<p>
    <a href="https://orcid.org/0000-0001-7552-9812">Dr. Thomas Pope </a>
    <br>
    <a href="https://ncl.ac.uk/nes/people/profile/tompenfold.html">Prof. Thomas Penfold </a>
    <br>
    <a href="https://orcid.org/0000-0002-8409-6702">Dr. Yvelin Giret </a>
</p>

<p>
    <a href="http://penfoldgroup.co.uk">Penfold Group </a> @ <a href="https://ncl.ac.uk">Newcastle University </a>
</p>

<p>
    <a href="#setup">Setup</a> • <a href="#getting">Quickstart</a> • <a href="#publications">Publications</a>
</p>

</td></tr></table>

#
CharTED-KMC is a Kinetic Monte-Carlo code used to model charge transport in organic semiconductors.

​The code has the following capabilities:

    -- Transport of electrons, holes and excitons.
    -- Interplay between multiple excited states on each site
    -- Miller-Abrahams and Marcus rates
    -- Order parameters to control relative orientation and order within thin film.

</a>

## SETUP

The quickest way to get started with CharTED-KMC is to clone this repository:

<!---
```
git clone https://github.com/thomasjamespope/CharTED-KMC/tree/main 
```
--->

This contains all the source files and an example input file - as well as the offline manual.

CharTED-KMC is entirely stand-alone and comes with its own Makefile. The compile, simply go to the code directory and type:

<!---
```
make all
```
--->

Now you're good to go!

## GETTING STARTING 
|------------------|-----------|------------------|--------------------------------------------------------------------------------------------------------------------------------|
| DEBUG            | 1 boolean | .false.          | Sets the seed for the random number generator to zero for reproducable calculations                                            |
| NUMBER_STEPS     | 1 integer | 1,000,000        | Sets the number of events chosen before the calculation is terminated. This can be set to -1 for an ad infinitum calculation   |
| N_PRINT_STEP}    | 1 integer | 1,000}           | &Number of hopping events allowed to elapse between writing the various OUTPUT files, with the exception of the tajectory file |
| N_PRINT_STEP_XYZ | 1 integer | 10,000           | &Number of hopping events allowed to elapse between writing the trajectory file                                                |
| N_EQUILIBRATION  | 1 integer | 2 X N_PRINT_STEP | &The number of flux readings used in the calculation of the charge                                                             |
|------------------|-----------|------------------|--------------------------------------------------------------------------------------------------------------------------------|


## LICENSE

This project is licensed under the GPL-3.0 License - see the LICENSE.md file for details.

## PUBLICATIONS

[T. Pope](https://orcid.org/0000-0001-7552-9812), [Y. Giret](https://orcid.org/0000-0002-8409-6702), M. Fsadni, [P. Docampo](https://orcid.org/0000-0001-6164-4748), [C. Groves](https://orcid.org/0000-0003-2402-1618), & [TJ. Penfold](https://orcid.org/0000-0003-4490-5672) (2023). *[Modelling the effect of dipole ordering on charge-carrier mobility in organic semiconductors](https://doi.org/10.1016/j.orgel.2023.106760)*, *Organic Electronics*, **115**, 106760.
