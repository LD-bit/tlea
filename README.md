# Repository of the M.Sc. thesis, code name: tlea
M.Sc. thesis topic: Orbit determination of OPS-SAT-1 based on Doppler shifts observations

## Content
MATLAB scripts for the study of an alternative method to TLEs to obtain the position of ESA's first CubeSat, OPS-SAT, with the aid of Doppler shift measurements and Satellite Laser Ranging (SLR) measurements (range and directional angles of the antenna) from two and more ground stations: ESOC-1 (Doppler solely) and IZN-1 (ranging + angles).

The study is two-fold:
- a State of the Art implementation based on Labsir (2020) PhD thesis on Lie groups: the evolution model of the Extended Kalman Filter are based on his work (the observation model has been defined w.r.t. the observers, i.e. the ground stations).
- a contribution of the thesis with various scenario of observation (meaning different observation functions), detailed in Chapter 4 of the M.Sc. thesis.

## How to use

### Preliminary: preliminary and installation
- Developped on MATLABR2023a. I haven't tried myself but it should work on anterior versions as well, up to a certain point.
- Better if you use `git` in your MATLAB environment. I let you follow the instructions in this [tutorial](https://fr.mathworks.com/help/matlab/matlab_prog/set-up-git-source-control.html), it explains in a clear manner how to do.

### State of the Art solution
The scripts for this section are:
- `ClassicEKF.m`: main code
- `evolution_model/fk.m`: evolution model function (dynamic of the satellite)
- `evolution_model/Jfk.m`: Jacobian matrix function of the evolution model
- `observation_model/hk.m`: observation function
- `observation_model/Jhk.m`: Jacobian matrix function of the observation model

### M.Sc. contributions (heart of tlea project)
The scripts for this section are:
- `DopplerEKF.m`: main code of contribution
- `evolution_model/fk.m`: evolution model function (dynamic of the satellite)
- `evolution_model/Jfk.m`: Jacobian matrix function of the evolution model
- `observation_model/Doppler_hk.m`: observation function with Doppler functions (4 functions to select from in total)
- `observation_model/Doppler_Jhk.m`: Jacobian matrices function for Doppler scenarios (idem, 4 functions in total to select from)
- `Comparison.m`: the comparison between the state of the art solution from Labsir and tlea's implementation

A selection among the 4 different observation functions can be made. Then the code has to be modified in the observation function and its Jacobian expression and it is required to tailor the covariance matrices in the main code of this section (`DopplerEKF.m`).

### Additional content
Additional scripts were implemented to train myself on Kalman Filters (`LinearKF.m` and `BabyExtendedEKF.m`).

The scripts referred as `SensitivityClassicEKF.m`and `SensitivityDopplerEKF.m` were exploratory studies on the sensitivity of the evolution model (work on uncertainties). It is also mentioned in Chapter 4, briefly when studying the Classical approach to the problem.

## How to cite
Dubreil, L. (2023) Precise orbit determination of OPS-SAT-1 using Doppler shifts observations from ESOC-1.

The M.Sc. manuscript is available on ResearchGate: [Dubreil L., 2023](https://www.researchgate.net/publication/371139968_Precise_orbit_determination_of_OPS-SAT-1_using_Doppler_shifts_observations_from_ESOC-1)

## Key words
Extended-Kalman-Filter, TRANSIT, Doppler shift, Radar, OPS-SAT, OPS-SAT Space Lab, exp235

