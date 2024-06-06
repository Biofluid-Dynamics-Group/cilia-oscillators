# Filament Oscillator

This repository contains the source code and documentation for the filament oscillator model. The project aims to study the dynamics and coordination of cilia.

## Main simulation files

||Number of cilia|Beat plane angle|Background flow|Sync|Bending stiffness|Runs?
|---|---|---|---|---|---|---|
|01|1|0|0|0|1000|Yes|
|02|1|0|0|0|0.001|Yes|
|03|1|0|uniform|0|1|Yes|
|04|1|0|oscillatory|0|1|Yes|
|05|1|0|oscillatory|flow out of phase|1|Yes|
|06|1|12.5°|0|0|1|Yes|
|07|2|0|0|0|1000|Slightly wobbly $\psi_1$, $\psi_2$ is very oscillatory|
|08|2|0|0|0|1|Breaks at the end of the ninth period, and $\psi_2$ varies wildly|
|09|2|0|0|random cilia start|1|No|
|10|2|12.5°|0|random cilia start|1|No|
