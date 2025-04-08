# DEWPython

[![License: GPL v3 License](https://img.shields.io/badge/License-GPL--3.0-blue.svg?style=flat-square)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![DOI](https://zenodo.org/badge/290678733.svg)](https://zenodo.org/badge/latestdoi/290678733)

## Overview

DEWPython is a Python implementation of the Deep Earth Water (DEW) model that allows users to compute thermodynamic and elastic properties of various aqueous species across a wide range of temperatures (100-1200°C) and pressures (1-60 kbar). The package is based on the [DEW spreadsheet](http://www.dewcommunity.org/) and provides similar functionality with enhanced features, including integrated SUPCRT support.

DEWPython calculates Gibbs energies of formation, equilibrium constants, and standard volume changes for reactions in extreme pressure-temperature conditions, making it particularly valuable for modeling geochemical processes in planetary interiors, hydrothermal systems, and deep-Earth environments.

## Key Features

- Calculate thermodynamic properties across extreme pressure-temperature ranges
- **Integrated support for SUPCRT96 and SUPCRTBL** directly within the package (executables included)
- Import and compare species between DEW and SUPCRT models
- Built-in database of minerals and aqueous species from slop16.dat
- Interactive and automated input options
- Plotting utilities for visualizing thermodynamic properties

## Applications

DEWPython has been successfully applied to studying:

- Thermodynamic constraints on the citric acid cycle and related reactions in ocean world interiors
- Chemical reaction networks in extreme environments
- Stability of organic compounds under high pressure-temperature conditions
- Biogeochemical processes relevant to astrobiology and prebiotic chemistry

## Installation

### Using pip
```bash
pip install DEWPython
```

### From source
```bash
git clone https://github.com/chandr3w/DEWPython.git
cd DEWPython
pip install -e .
```

## Quick Start

```python
import DEWPython
from DEWPython import DEWModel as dm

# Initialize a DEW object
reaction = dm.DEW()

# Define temperature and pressure arrays
pt_arr = [[100, 1000, 10], [1, 60, 1]]  # [T_min, T_max, T_step], [P_min, P_max, P_step] (°C, kbar)

# Run calculation with automatic inputs
reaction.run(pt_arr, 
             min_inp=[], 
             aq_inp=[['CO2,AQ', '4'], ['H2,aq', '5']], 
             g_inp=[], 
             h2o_inp=0,
             min_out=[], 
             aq_out=[['Oxaloacetate2-', '1'], ['H+', '2']], 
             g_out=[], 
             h2o_out=3,
             ptInp='Regular', 
             rhoWat='Z&D 2005', 
             forceBool=False, 
             dieEQ='Sverjensky', 
             forceSC=True, 
             WFEQ='D&H 1978', 
             dsV=False, 
             pdsV=False, 
             DV=False)

# Export results to CSV
reaction.export_to_csv()
```

## Using SUPCRT Within DEWPython

DEWPython includes built-in SUPCRT functionality with both SUPCRT96 and SUPCRTBL executables packaged within the distribution:

```python
import DEWPython
from DEWPython import DEWModel as dm

# Initialize a DEW object
reaction = dm.DEW()

# Run SUPCRT interactively
reaction.run_supcrt()  # Uses SUPCRT96 by default
# Or specify version
reaction.run_supcrt(version='96')  # SUPCRT96 - supports Maier-Kelly heat capacity and Powell & Holland 1990 equations
# reaction.run_supcrt(version='BL')  # SUPCRTBL - includes updated thermodynamic dataset with Holland & Powell (2011) properties

# Process SUPCRT output
reaction.calculate_supcrt(customFile="output_file.out")

# Export processed SUPCRT results
reaction.export_supcrt_to_csv()
```

## SUPCRT Variants

DEWPython includes multiple SUPCRT variants:

- **SUPCRT96**: An improved version of SUPCRT92 (Johnson et al., 1992) that supports both Maier-Kelly heat capacity expressions (Cp = a + b·T + c·T⁻²) and the four-term Powell & Holland (1990) equation (Cp = a + b·T + c·T⁻² + d·T⁻⁰·⁵)

- **SUPCRTBL**: A further enhanced version (Zimmer et al., 2016) that updates the thermodynamic dataset with:
  - Properties from Holland & Powell (2011) for mineral end-members
  - Added As-acid and As-metal aqueous species
  - Updated properties for Al-bearing species and other minerals relevant to geological carbon sequestration

## Complementary Software

While not directly integrated into DEWPython, these complementary tools are often used alongside it for comprehensive thermodynamic modeling:

- **SeaFreeze**: A separate package for calculating thermodynamic properties of ice polymorphs (Ih, II, III, V, VI and ice VII/ice X) and liquid water at conditions relevant to ocean worlds (130-500 K, 0-2300 MPa). SeaFreeze can be used to determine water phase boundaries when modeling planetary interiors. Available at [github.com/Bjournaux/SeaFreeze](https://github.com/Bjournaux/SeaFreeze).

## Detailed Usage

### Interactive Mode

DEWPython can be used interactively, prompting for inputs:

```python
reaction = dm.DEW()
reaction.set_inputs()   # Set input species
reaction.set_outputs()  # Set output species
reaction.set_preferences()  # Set calculation options
reaction.set_TPRho()    # Set temperature and pressure conditions
reaction.calculate()    # Perform calculations
reaction.make_plots()   # Generate plots
```

### Automated Mode

For automated or scripted use:

```python
reaction = dm.DEW()
reaction.run(pt_arr, min_inp, aq_inp, g_inp, h2o_inp, min_out, aq_out, g_out, h2o_out, 
             ptInp, rhoWat, forceBool, dieEQ, forceSC, WFEQ, dsV, pdsV, DV)
```

### Searching for Species

To find available species in the database:

```python
dm.search("string")  # Search for species containing "string"
```

## Example: Modeling Ocean World Reactions

The following example models the formation of oxaloacetate from CO₂ and H₂ under conditions relevant to ocean worlds:

```python
import DEWPython
from DEWPython import DEWModel as dm
import numpy as np

# Initialize model
reaction = dm.DEW()

# Define PT conditions relevant to ocean worlds (e.g., Enceladus)
pt_arr = [[100, 500, 20], [1, 10, 1]]  # [T_min, T_max, T_step], [P_min, P_max, P_step]

# Run calculation
reaction.run(pt_arr, 
             min_inp=[], 
             aq_inp=[['CO2,AQ', '4'], ['H2,aq', '7']], 
             g_inp=[], 
             h2o_inp=0,
             min_out=[], 
             aq_out=[['Oxaloacetate2-', '1'], ['H+', '2']], 
             g_out=[], 
             h2o_out=3,
             ptInp='Regular')

# Export results
reaction.export_to_csv()
```

## Theoretical Background

The DEW model extends the Helgeson-Kirkham-Flowers (HKF) equations of state to calculate thermodynamic properties at high temperatures (373-1473 K) and pressures (1-6 GPa). DEWPython implements:

1. Water density calculations using Zhang & Duan (2005, 2009)
2. Dielectric constant calculations using multiple equations (Johnson & Norton, Franck, Fernandez, Sverjensky)
3. Gibbs free energy calculations using Delaney & Helgeson (1978) or numerical integration
4. Born coefficients for aqueous species

These calculations enable accurate prediction of reaction properties under extreme conditions relevant to planetary interiors and deep-Earth environments.

## Range of Validity

| Model | Pressure Range | Temperature Range | Integration Status |
|-------|----------------|-------------------|-------------------|
| DEW   | 0.1-6 GPa | 373-1473 K | Core functionality |
| SUPCRT96 | 1-5000 bar | 273-873 K | Integrated in package |
| SUPCRTBL | 1-5000 bar | 273-1273 K | Integrated in package |
| SeaFreeze | 0.0001-2.3 GPa | 130-500 K | Separate complementary tool |

## Recent Applications

Recent work by Işık et al. (2025) has applied DEWPython to study thermodynamic constraints on the citric acid cycle and related reactions in ocean world interiors. This research demonstrated the feasibility of using DEWPython to:

1. Quantify the thermodynamic viability of metabolic reactions under pressure-temperature conditions in ocean worlds
2. Calculate equilibrium constants and Gibbs free energy changes for reactions in the TCA cycle
3. Model prebiotic reaction networks at high pressures
4. Identify potential energy bottlenecks in metabolic pathways on different planetary bodies

In this research, SeaFreeze was used as a complementary tool to determine water phases along PT profiles, helping to exclude icy regions from calculations for aqueous species.

## Dependencies

* Pandas
* NumPy
* Matplotlib
* Collections
* JSON

## References

* Huang, F., & Sverjensky, D. A. (2019). Extended Deep Earth Water Model for predicting major element mantle metasomatism. Geochimica et Cosmochimica Acta, 254, 192-230.
* Sverjensky, D. A. (2019). Thermodynamic modelling of fluids from surficial to mantle conditions. Journal of the Geological Society, 176(2), 348-374. https://doi.org/10.1144/jgs2018-105
* Facq, S., Daniel, I., Montagnac, G., Cardon, H., & Sverjensky, D. A. (2016). Carbon speciation in saline solutions in equilibrium with aragonite at high pressure. Chemical Geology, 431, 44-53.
* Johnson, J.W., Oelkers, E.H. & Helgeson, H.C. (1992). SUPCRT92 - A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 bar to 5000 bar and 0°C to 1000°C. Computers and Geosciences 18:899-947.
* Sverjensky, D. A., Harrison, B., & Azzolini, D. (2014). Water in the deep Earth: The dielectric constant and the solubilities of quartz and corundum to 60 kbar and 1200°C. Geochimica et Cosmochimica Acta, 129, 125-145.
* Işık, S., Melwani Daswani, M., Işık, E., Weber, J., & Olgun Kıyak, N. (2025). Thermodynamic constraints on the citric acid cycle and related reactions in ocean world interiors.
* Chan, A., Melwani Daswani, M., Vance, S. (2020). DEWPython: A Python Implementation of the Deep Earth Water Model and Application to Ocean Worlds.
* Zimmer, K., Zhang, Y.L., Lu, P., Chen, Y.Y., Zhang, G.R., Dalkilic, M. & Zhu, C. (2016). SUPCRTBL: A revised and extended thermodynamic dataset and software package of SUPCRT92. Computer and Geosciences 90:97-111.
* Journaux, B., Brown, J.M., Pakhomova, A., Collings, I.E., Petitgirard, S., Espinoza, P., Boffa Ballaran, T., Vance, S.D., Ott, J., Cova, F., Garbarino, G., Hanfland, M. (2020). Holistic Approach for Studying Planetary Hydrospheres: Gibbs Representation of Ices Thermodynamics, Elasticity, and the Water Phase Diagram to 2,300 MPa. Journal of Geophysical Research: Planets, 125, e2019JE006176.
* Holland, T.J.B. & Powell, R. (2011). An improved and extended internally consistent thermodynamic dataset for phases of petrological interest, involving a new equation of state for solids. Journal of Metamorphic Geology, 29, 333-383.

## Authors

* **Andrew Chan** - *Div. of Geological and Planetary Sciences, California Institute of Technology, Pasadena, CA, USA 91125* 
* **Mohit Melwani Daswani** - *Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA 91109*
* **Steven Vance** - *Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA 91109*
* **Seda Işık** - *Eurasia Institute of Earth Sciences, Istanbul Technical University, 34469 Istanbul, Türkiye*
* **Emre Işık** - *Max Planck Institute for Solar System Research, 37077 Göttingen, Germany*

## License

Copyright (C) 2025. Jet Propulsion Laboratory, California Institute of Technology. Government sponsorship acknowledged.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Acknowledgments

M.M.D. was supported by the NASA Planetary Science Early Career Award Program NNH19ZDA001N-ECA to proposal #19-ECA19-0032. A part of this research was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (80NM0018D0004). Financial support for A.C. was provided by the Jet Propulsion Laboratory, California Institute of Technology Summer Undergraduate Research Fellowship program, and the Caltech Associates. S.I. acknowledges funding by the Scientific and Technological Research Council of Türkiye (TÜBITAK), under grant 122F287. 

## Previous Version History

V 2.0.0

The DEW package allows the user to compute the thermodynamic and elastic properties of various aqueous inputs for a general range of 100 - 1200 C and a pressure range of 1 - 60 kbar. It is based on the [DEW spreadsheet](http://www.dewcommunity.org/) and behaves similarly. The DEW package additionally provides integrated support for [SUPCRTBL](https://models.earth.indiana.edu/supcrtbl.php) and can be used to directly import and compare species between the two models.

V 2.0.1 (this version; not yet on pip)

Automatic input function updated (2022/10/30). Now you can use it in the form as shown in `tools/Example-dew-calc.ipynb`. Beware: this is a beta release. You can use it only by cloning this repository into your local disk, then importing DEWPython according to the location of the cloned repository. 



