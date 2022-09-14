# UAV Bistatic SAR

uav-bistatic-sar is the MATLAB repository for the PoliMi project of a demonstrator for a UAV-based bistatic low cost SAR system.

## Requirements

MATLAB with installed Phased Array System Toolbox and Parallel Computing Toolbox. The current version of TDBP focusing implements gpuArray and it s being developed a version using CUDA acceleration.

## Folders explanation

### SDR_scripts

Here are presents all the scripts used to run the SDR board for acquisition of the data. The most important files are inside `gnuradio` and are used to start transmission and reception of samples.

### lib

Here are all the functions often needed in the other scripts, mainly for storing and loading binary complex data in an interleaved format.

### mat_files

Here should be copied the dataset to be downloaded independently from [here].

### range compression

Here are presents the script needed to perform the FFT based range compression on the raw binary data collected by the SDR.

### tx_waveforms

Here are presents all waveform transmitted in the various experiment

### uav processing

Here is presentthe main script used to elaborate the RangeCompressed data. First it is time and phase corrected using navigation data, then it is focused with TDBP algorithm on a specified grid.
