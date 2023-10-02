
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

```
git clone https://github.com/thomasjamespope/CharTED-KMC/tree/main 
```

This contains all the source files and an example input file - as well as the offline manual.

CharTED-KMC is entirely stand-alone and comes with its own Makefile. The compile, simply go to the code directory and type:

```
make all
```

Now you're good to go!

## GETTING STARTING 
To run CharTED-KMC, type:

```
mpirun -np 4 CharTED-KMC.x input.file > output.o 2> error.o 
```

input.file contains the system parameters and is supplied by the user. output.o and error.o contain, respectively, the standard output and the sdandard error output for the calculation. Additional output files are also written. The MPI envirment is used to generate multiple independant simulations (in this case, 4) using the same input parameters. At the end of the calculation, statistical analysis is performed of the group of simulations generate physical predictions (like mobility, excitation rate, etc...). It is important to run as many simulations as possible or at the very least 2.

---


| Input parameter | Options | Default | Description |
| --- | --- | --- | --- |
| `DEBUG` | *1 boolean* | .false. | Sets the seed for the random number generator to zero for reproducable calculations |
| `NUMBER_STEPS` | *1 integer* | 1,000,000 | Sets the number of events chosen before the calculation is terminated. This can be set to -1 for an ad infinitum calculation |
| `N_PRINT_STEP` | *1 integer* | 1,000 | Number of hopping events allowed to elapse between writing the various OUTPUT files, with the exception of the tajectory file |
| `N_PRINT_STEP_XYZ` | *1 integer* | 10,000 | Number of hopping events allowed to elapse between writing the trajectory file |
| `N_EQUILIBRATION` | *1 integer* | 2 X N\_PRINT_STEP | The number of flux readings used in the calculation of the charge |
| `DIMENSIONS` | *3 integers* | no default | The number of grid points in x, y and z |
| `CORRELATED_NEIGHBORS` | *1 integer* | 7 | Number of nearest neighbours taken into account to calculate the electrostatic interaction |
| `COUL_CUTOFF` | *1 real number* | 15.0 | Cut-off radius for the Coulomb exclusion zone. The coulomb term is only recalculated is a hopping event ocurrs with the cutoff radius of a given site |
| `CALC_TYPE` | *HOLE*, *ELEC* or *FULL* | HOLE | Determines whether holes, electrons or both are considered in the system. If both are considered, singlet and triplet states are also considered |
| `CHARGE_RATE` | *MILLER-ABRAHAMS* or *MARCUS* | MARCUS | Determines type of rate equation you're using. The low temperature regime can be studied with Miller-Abrahams, and the high temperature by Marcus. Note, the high field regime is poorly described by both, so be careful! |
| `LATTICE` | *1 real number* | 1.0 nm | The distance between nearest-neighbour grid points in nm |
| `REORGANIZATION` | *1 real number* | 0.1 eV | Reorganization energy in eV for the marcus rate calculation. This parameter is ignored if you're using Miller-Abraham rates. |
| `TEMPERATURE` | *1 real number* | 300.0K | Temperature in Kelvin |
| `PERMITTIVITY` | *1 real number* | 3.5 AU | Relative permittivity in Atomic Units |
| `LOCAL` | *1 real number* | 0.1 eV | Charge Localization Factor in eV |
| `ELECTRIC_FIELD` | *1 real number* | 0.05 V/nm | Strength of the applied electric field in V/nm |
| `UNCORRELATED_DISORDER` | *1 real number* | 0.05 eV | Uncorrellated disorder strength in eV |
| `CORRELATED_TYPE` | *RANDOM*, *ORDERED* or *MIXED* | RANDOM | Determines whether the dipole moments are randomly distributed, or whether they are organized in an antiferromagnetic-type configuration. If MIXED is chosen, a scaling factor is used to determine how rigorous the ordering is (see: ORDER_PARAMETER) | 
| `ORDER_PARAMETER` | *1 real number* | 0.0 | Disorder from perfect antiferro symmetry (0 = no disorder) to completely random (1 = random) |
| `DIPOLE_MOMENT` | *1 real number* | 0.5 D | Strength of the average dipole moment in Debye |
| `FERMI_ENERGY` | *1 real number* | 0 eV | Fermi energy in eV |
| `HOLE_DENSITY` | *1 real number* | 0.001 $\text{site}^{-1}$ | Hole density per site |
| `HOLE_HOP` | *1 real number* | 1.0e12 $\text{s}^{-1}$ | Hopping attempt frequency for the holes. This is included in the prefactor to the rate equations and incorporates coupling between sites. It can be understood as an empirical parameter to tune. |
| `HOMO_ENERGY` | *1 real number* | 0 eV | Energy of the highest occupied atomic orbital in eV |
| `HOMO+n_ENERGY` | *1 real number* | inactive | Energy of additions HOMO levels in eV. Here, \emph{n} must incrementally increase from 1 for each additional level. |
| `ELEC_RELAX` | *1 real number* | 1.0e15 $\text{s}^{-1}$ | Electron relaxation rate. This is the rate at which a electron in a higher HOMO level will relax to the lowest. It is typically a few orders of magnetude faster than the hopping attempt frequency. |
| `ELEC_DENSITY` | *1 real number* | 0.001 $\text{site}^{-1}$ | Electron density per site |
| `ELEC_HOP` | *1 real number* | 1.0e12 $\text{s}^{-1}$ | Hopping attempt frequency for the electrons. This is included in the prefactor to the rate equations and incorporates coupling between sites. It can be understood as an empirical parameter to tune. |
| `LUMO_ENERGY` | *1 real number* | 0 eV | Energy of the lowest unoccupied atomic orbital in eV |
| `LUMO+n_ENERGY` | *1 real number* | inactive | Energy of additions LUMO levels in eV. Here, \emph{n} must incrementally increase from 1 for each additional level. |
| `ELEC_RELAX` | *1 real number* | 1.0e15 $\text{s}^{-1}$ | Electron relaxation rate. This is the rate at which a electron in a higher HOMO level will relax to the lowest. It is typically a few orders of magnetude faster than the hopping attempt frequency. |




## LICENSE

This project is licensed under the GPL-3.0 License - see the LICENSE.md file for details.

## PUBLICATIONS

[T. Pope](https://orcid.org/0000-0001-7552-9812), [Y. Giret](https://orcid.org/0000-0002-8409-6702), M. Fsadni, [P. Docampo](https://orcid.org/0000-0001-6164-4748), [C. Groves](https://orcid.org/0000-0003-2402-1618), & [TJ. Penfold](https://orcid.org/0000-0003-4490-5672) (2023). *[Modelling the effect of dipole ordering on charge-carrier mobility in organic semiconductors](https://doi.org/10.1016/j.orgel.2023.106760)*, *Organic Electronics*, **115**, 106760.
