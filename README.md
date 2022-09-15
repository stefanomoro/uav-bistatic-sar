# UAV Bistatic SAR

uav-bistatic-sar is the MATLAB repository for the PoliMi project of a demonstrator for a UAV-based bistatic low cost SAR system.

## Requirements

MATLAB with installed Phased Array System Toolbox and Parallel Computing Toolbox. The current version of TDBP focusing implements `gpuArray` and it s being developed a version using CUDA acceleration.

## Folders explanation

### SDR_scripts

Scripts used to run the SDR board for acquisition of the data. The most important files are inside `gnuradio` and are used to start transmission and reception of samples.

### lib

Functions often needed in the other scripts, mainly for storing and loading binary complex data in an interleaved format.

### mat_files

Here should be copied the dataset to be downloaded independently.

### range compression

Scripts needed to perform the FFT based range compression on the raw binary data collected by the SDR.

### tx_waveforms

Waveforms transmitted in the various experiment, in different formats.

### uav processing

Main script used to elaborate the RangeCompressed data. First it is time and phase corrected using navigation data, then it is focused with TDBP algorithm on a specified grid.
