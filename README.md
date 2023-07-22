# radar-vad
### Quick summary ###

This is a C module for retrieve horizontal wind profile using a Velocity-Azimuth Display
(VAD) technique.

This wind retrieval code is licensed under GNU GENERAL PUBLIC LICENSE (GNU GPL).

If you use this software in publication,
an acknowledgment and citation of the following paper(s) would be appreciated.
If you have any comments, suggestions for improvements, bug fixes or you need help to
interface the software with your radar data, please contact us.

Please cite the following paper:

Oue, M., B. A. Colle, S. E. Yuter, P. Kollias, P. Yeh, L. M. Tomkins, 2023: Microscale Updrafts Within Northeast U.S. Coastal Snowstorms Using High-Resolution Cloud Radar Measurements, Mon. Wea. Rev. in review.

### Related publications ###

Kalesse, H., G. de Boer, A. Solomon, M. Oue, M. Ahlgrimm, D. Zhang, M. Shupe, E. Luke, and A. Protat, 2016: Understanding rapid changes in phase partitioning between cloud liquid and ice in stratiform mixed-phase clouds: An Arctic Case Study. Mon. Wea. Rev., 144, 4805-4826, doi: 10.1175/MWR-D-16-0155.1.

### Summary of set up ###

This module requres NetCDF4. Then, excute:
make

VAD2 is the excutable file.

Example of excute:
VAD2 [PPI data file path] [Velocity sign convention] [name for velocity] [name of SNR] [SNR threshold] [Azimuth offset]

A sample script is available in ./examples


### Contacts ###

* [Mariko Oue](mailto:mariko.oue@stonybrook.edu)