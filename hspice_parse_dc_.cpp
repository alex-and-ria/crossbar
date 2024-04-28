
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
#include<iostream>//cout; cin
#include<fstream>//fstream
#include<cstring>//strcmp
#include <stdlib.h>/* strtod */

//#define M 3
//#define N 2
//#define N_swp 2

//#define ifn "hspice/3x2_2_dc_sweep_gen.sw0"
#define buff_sz 128

int main(int argc, char** argv){
	unsigned int M,N,N_swp; char* ifn;
	if(argc<5){
		std::cout<<"\nusage: [exec_name] M N N_swp ifn";
	
	}
	else{
		M=atoi(argv[1]); N=atoi(argv[2]); N_swp=atoi(argv[3]); ifn=argv[4];

	std::fstream ifs(ifn,std::ios::in|std::ios::binary);
	char ifs_buff[buff_sz]={}, f_num_buff[buff_sz]; char q;
	char* p_sub_str;
	double **voltages=new double*[2*M*N];//define array on heap to avoid using all stack; parse only node voltages for memristors (2*M*N); Vin should be provided by hspice_gen_dc.cpp
	for(unsigned int i=0;i<(2*M*N);i++){
		voltages[i]=new double[N_swp];
	
	}
	if(ifs.is_open()){
		ifs>>ifs_buff;//ifs.getline(ifs_buff,buff_sz,' ');
		while(strcmp(ifs_buff,"$&%#")!=0){//data starts after "$&%#" delimiter;
			//std::cout<<"\nifs_buff="<<ifs_buff;
			ifs>>ifs_buff;
		
		}
		unsigned int nd_cnt=0,swp_cnt=0;
		while(1){
			ifs>>ifs_buff;
			if(ifs.eof()) break;
			
			
			
			
			
			
			p_sub_str=strtok(ifs_buff,"."); char frst_sgn[3];//should have first character which represents the sign ('0' for '+' or '-' for '-') of the number parsed, dot ('.') separator and terminating '\0' character for convenient usage of strcat;
			frst_sgn[0]=p_sub_str[0]; frst_sgn[1]='.'; frst_sgn[2]='\0';
			p_sub_str=strtok(NULL,".");//get the rest of the number and new sign character;
			//std::cout<<"\nfrst_sgn="<<frst_sgn<<"\tp_sub_str="<<p_sub_str;
			while(p_sub_str!=NULL){
				f_num_buff[0]='\0';//clear the buffer;
				strcat(f_num_buff,frst_sgn);//set the sign;
				strcat(f_num_buff,p_sub_str);//copy the rest of the number;
				frst_sgn[0]=p_sub_str[strlen(p_sub_str)-1];//get new frst_sgn of the next number;
				p_sub_str=strtok(NULL,".");
				if(p_sub_str!=NULL){//if it was not the last number in the line
					f_num_buff[strlen(f_num_buff)-1]='\0';//clear the last character (which is the new frst_sgn of the new number);
				
				}
				
				
				//std::cout<<"\nf_num_buff="<<f_num_buff<<"\tfrst_sgn="<<frst_sgn<<"\tp_sub_str="<<p_sub_str;
				
				//std::cout<<"\nf_num_buff="<<f_num_buff<<"\tfloat_val="<<strtod(f_num_buff,NULL);
				
				if(nd_cnt>=2 && nd_cnt<2*M*N+2){//first two numbers are skipped, corresponding node voltages start from the third number; after required voltages recorded, skip the rest of values until new sweep values;
					voltages[nd_cnt-2][swp_cnt]=strtod(f_num_buff,NULL);
				
				}
				if(nd_cnt==2+2*M*N+2*M-1){//first two values are skipped, 2*M*N node voltages should be recoreded, M voltages (input voltages) and M currents are skipped;
					swp_cnt++; nd_cnt=0;
				
				}
				else{
					nd_cnt++;
				
				}
				
				
				
				//std::cin>>q;
			
			}
			
		
		}
		
		
		ifs.close();
		
		std::cout<<"\nHspcOut=[";
		for(unsigned int i=0;i<N_swp;i++){
			std::cout<<'[';
			for(unsigned int j=0;j<2*M*N-1;j++){//nodes 1 to 2*M*N-1;
				std::cout<<voltages[j][i]<<',';
				
			
			}
			if(i==N_swp-1){//i==N_swp-1; j==2*M*N-2;
				std::cout<<voltages[2*M*N-1][i]<<"]];";
			
			}
			else{
				std::cout<<voltages[2*M*N-1][i]<<"];";
			
			}
		
		}
		
	
	}
	else{
	
	std::cout<<"\nifs.is_open()";
	}

	for(unsigned int i=0;i<2*M*N;i++){
		delete[] voltages[i];
		
	}
	delete[] voltages;
	

}
return 0;

}
