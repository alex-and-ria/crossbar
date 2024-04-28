































#!/bin/bash
#put in the same directory as hspice_gen_dc.cpp and hspice_parse_dc.cpp
set -x #debugging mode; print each command executed;
set -e #exit on error;

M=3
N=2
N_swp=2
n_stat=1

gen_file=hspice_gen_dc #.cpp
parse_file=hspice_parse_dc_ #.cpp
lang_ext=.cpp
hsp_gen_f_ext=.sw0

parse_ifn=q.lis
matlab_ifn=q.m #matlab template file name; should be in the folder;
matlab_ofn=q1.m #generated matlab source file;
output_log_fn=q_out.txt #file to output results;
log_fl=log_fl.txt
echo "" > "$output_log_fn" #clear log_file;
fl=0

#for((N=2;N<=2;N*=2)) do
echo -e "" >> "$output_log_fn" #output just new line to indicate new record in the csv file;
#for((M=2;M<=2;M*=2)) do
#echo -e "*****\n" >> "$output_log_fn"
hspice_stat=0;
lu_stat_d=0; lu_stat_sol=0;
lib_lu_stat_d=0; lib_lu_stat_sol=0;
for((ii=0;ii<n_stat;ii++)) do

g++ $gen_file$lang_ext -o $gen_file #compile hspice file generator;
gen_out=$(./$gen_file $M $N $N_swp) #command substitution; save command output to variable;



#echo "$gen_out"


IFS=$'\n'
#printf "%d" "'$IFS"

read -d '' -a gen_out_arr < <(echo "$gen_out") || true #probably can be used directly with ./$gen_file $M $N $N_swp instead, but for better readability it is done with process substitution; stores data in gen_out_arr array, splited with $IFS; continue untill all strings are processed; to avoid exiting on error || true is used; returns 1 presumably because EOF is reached: -d '' is used;
#gen_out_arr[0] should store conductances; gen_out_arr[1] should store file path of generated file;
#read  -d '' -a array <<< "$gen_out"






IFS='/' #used to parse file path in gen_out_arr[1];
read -a gen_out_fn_arr <<<  "${gen_out_arr[1]}" #split path in array;
#echo "${gen_out_fn_arr[2]}"
#gen_out_fn_arr[0] should be .
#gen_out_fn_arr[1] should be hspice
#gen_out_fn_arr[2] should be filename of hspice file generated by gen_file;

cd "${gen_out_fn_arr[1]}"
# Synopsys
export SNPSLMD_LICENSE_FILE=27020@saratoga.cse.sc.edu
#export PATH="/share/reconfig/synopsys/hspice/P-2019.06-2/hspice/bin:$PATH" #hspice
#export PATH="/share/reconfig/synopsys/cscope64/P-2019.06/ai_bin:$PATH" #scope;
/share/reconfig/synopsys/hspice/P-2019.06-2/hspice/bin/hspice "${gen_out_fn_arr[2]}" > "$parse_ifn"
echo -e "\n\n$m $n ">>"$log_fl"
if [ "$fl" -eq "0" ] 
then cat "$parse_ifn" >"$log_fl"
 else cat "$parse_ifn" >>"$log_fl"; let "$fl=1";
 fi
hspice_jtr=$(cat "$parse_ifn" | grep "job total runtime")
cd ..

#pwd

g++ $parse_file$lang_ext -o $parse_file #compile hspice file generator;
parse_out=$(./$parse_file $M $N $N_swp "./${gen_out_fn_arr[1]}/${gen_out_fn_arr[2]}$hsp_gen_f_ext") #command substitution; save command output to variable;

echo -e "m=$M; n=$N; ${gen_out_arr[0]}\n$parse_out" > "$matlab_ofn"
cat "$matlab_ifn" >> "$matlab_ofn"

matlab_out=$(~/abc/matlab/bin/matlab -batch "${matlab_ofn%.m}") #remove extension ".m" from $matlab_ofn file name and run matlab for this file in batch mode;

#echo -e "ii=$ii m=$M; n=$N;\nN_swp=$N_swp;\n$hspice_jtr\n$matlab_out\n\n" >> "$output_log_fn"
echo -e "$hspice_jtr"
#echo -e "$hspice_jtr" | grep -o -E "[0-90+\.[0-9]*" #-o output only maching part; -E extended regular expression;
hspice_stat="$hspice_stat"+$(echo -e "$hspice_jtr" | grep -o -E "[0-90+\.[0-9]*") #concatenate strings with '+' as delimiter;
#echo $(echo "$matlab_out" |awk '{print $1}') #parse matlab data (delimiter is whitespace);
lu_stat_d="$lu_stat_d"+$(echo "$matlab_out" |awk '{print $1}')
lu_stat_sol="$lu_stat_sol"+$(echo "$matlab_out" |awk '{print $2}')
lib_lu_stat_d="$lib_lu_stat_d"+$(echo "$matlab_out" |awk '{print $3}')
lib_lu_stat_sol="$lib_lu_stat_sol"+$(echo "$matlab_out" |awk '{print $4}')

done
#echo -e "hspice_stat=$hspice_stat"
hspice_mean=$(echo "($hspice_stat)/($n_stat.0)" | bc -l)
#hspice_var=$(echo "($hspice_stat)/($n_stat.0)" | bc -l)
lu_mean_d=$(echo "($lu_stat_d)/($n_stat.0)" | bc -l)
lu_mean_sol=$(echo "($lu_stat_sol)/($n_stat.0)" | bc -l)
lib_lu_mean_d=$(echo "($lib_lu_stat_d)/($n_stat.0)" | bc -l)
lib_lu_mean_sol=$(echo "($lib_lu_stat_sol)/($n_stat.0)" | bc -l)
echo -n "$hspice_mean,$lu_mean_d,$lu_mean_sol,$lib_lu_mean_d,$lib_lu_mean_sol,|" >> "$output_log_fn" #do not append '\n' to the output;
#done
#done
