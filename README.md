[![SCREAM Logo](https://portal.nersc.gov/cfs/e3sm/petercal/scream/SCREAM_Halloween_Logo.png)](https://e3sm.org)

Simple Cloud-Resolving E3SM Atmosphere Model (SCREAM)
================================================================================

SCREAM is a global atmosphere model targeted towards 3 km ("cloud resolving")
resolution. It is part of the Energy Exascale Earth System Model (E3SM) and
this git repo is forked from and follows https://github.com/E3SM-Project/scream

Please visit the [E3SM website](https://e3sm.org) for further details about E3SM. If you're part of the E3SM project, you can check in on SCREAM's internal happenings at [SCREAM confluence site](https://acme-climate.atlassian.net/wiki/spaces/NGDNA/overview).

Table of Contents 
--------------------------------------------------------------------------------
- [Quick Start](#quickstart)
- [Supported Machines](#supportedmachines)
- [Running](#running)
- [Contributing](#contributing)
- [Acknowledge](#acknowledge)
- [License](#license)

Quick Start
--------------------------------------------------------------------------------
The [Quick Start](https://e3sm.org/model/running-e3sm/e3sm-quick-start/) page 
includes instructions on obtaining the necessary code and input data for model 
setup and execution.

Supported Machines 
--------------------------------------------------------------------------------
E3SM is high-performance computing application and generally requires a cluster
with several hundred nodes to run a scientifically validated case at a useful
simulation speed.

To run E3SM, it is recommended that you obtain time on a 
[Supported Machine](https://e3sm.org/model/running-e3sm/supported-machines/).

Running
--------------------------------------------------------------------------------
Please refer to [Running E3SM](https://e3sm.org/model/running-e3sm/) 
 for instructions on running the model. 

Contributing
--------------------------------------------------------------------------------
Please refer to [Contributing](CONTRIBUTING.md) for details on our code development
process for submitting pull requests.

Acknowledgement
--------------------------------------------------------------------------------
The Energy Exascale Earth System Model (E3SM) Project should be acknowledged in
publications as the origin of the model using
[these guidelines](https://e3sm.org/resources/policies/acknowledge-e3sm/).

In addition, the software should be cited.  For your convenience,
the following BibTeX entry is provided.
```TeX
@misc{e3sm-model,
	title = {{Energy Exascale Earth System Model (E3SM)}},
	author = {{E3SM Project}},
	abstractNote = {{E3SM} is a state-of-the-art fully coupled model of the {E}arth's 
		climate including important biogeochemical and cryospheric processes.},
	howpublished = {[Computer Software] \url{https://dx.doi.org/10.11578/E3SM/dc.20180418.36}},
	url = {https://dx.doi.org/10.11578/E3SM/dc.20180418.36},
	doi = {10.11578/E3SM/dc.20180418.36},
	year = 2018,
	month = apr,
}
```

License
--------------------------------------------------------------------------------
The E3SM model became open development at the time of first model and data release.
Please see [LICENSE](LICENSE) for details.

