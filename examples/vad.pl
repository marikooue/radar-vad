#! /usr/bin/perl
# Example of VAD2

# If the data include an offset of azimuth measurement, use this option. Otherwise, use "0"
$AZ_OFFSET =  0;

# Convention of velocity sign
$sign = -1; #1: positive is going away; -1: negative is going away

# Parameter names
$para_v = "corrected_velocity"; #name of Doppler velocity
$para_s = "snr_copol_vtx"; #name of SNR
$s_thre = 5; # minimum value of para_s

# Input data file path
$data = "KASPR_PP_MOMENTS_20210201-141333.ppi.nc";
# Output file path
$odata = "sbukasprppivad.c1.20210201.141333.nc";

system("../src/VAD2 $data $sign $para_v $para_s $s_thre $AZ_OFFSET > tmp.txt");
print "$data ->\n$odata \n";
system("mv vad.nc $odata ");


