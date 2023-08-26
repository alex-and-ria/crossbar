































#!/bin/bash
#put in the same directory as hspice_gen.cpp and hspice_parse.cpp
#set -x #debugging mode; print each command executed;
set -e #exit on error;

M=5
N=7

gen_file=hspice_gen #.cpp
parse_file=hspice_parse #.cpp
lang_ext=.cpp

g++ $gen_file$lang_ext -o $gen_file #compile hspice file generator;
gen_out=(`./$gen_file $M $N`) #run hspice file generator and store output (with generated file name); store output in array (white space is a separator);

echo -e "\narr=${gen_out[@]} with ${#gen_out[@]} elements";
echo ${gen_out[1]} #file path prom current directory;

IFS='/' # Set '/' as delimiter to parse path;
read -a gen_out_arr <<< ${gen_out[1]} #split path in array;
echo ${gen_out_arr[@]}
IFS=' ' #reset IFS to white space;


cd ${gen_out_arr[1]} #cd to the folder containing generated hspice file, to keep files outed by hspice itself in the same folder;
hspice_out=`/share/reconfig/synopsys/hspice/P-2019.06-2/hspice/bin/hspice ${gen_out_arr[2]}`
hs_tot_time_arr=(`echo "$hspice_out" | grep -wn "job total runtime"`)

echo "$hspice_out"
echo -e "\nhspice tot_time=(${hs_tot_time_arr[@]})" # with ${#hs_tot_time_arr[@]} elements"


cd ../ #cd to the genrator and parser source files;




echo -e "\nInput:"; echo "m=$M; n=$N; ${gen_out[2]} ${gen_out[3]} ${gen_out[4]} ${gen_out[5]}"
echo -e "\nOutput:";
g++ $parse_file$lang_ext -o $parse_file
./$parse_file ${gen_out[1]}.ic0 $M $N

#echo "${gen_out[1]} with ${#gen_out[@]} elements"



