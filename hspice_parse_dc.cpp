































						  
#include<iostream>
#include<cstring>
#include<fstream>

//#define M 3
//#define N 3
//#define N_swp 2
#define buff_sz 256

//#define ifn "./hspice/q.lis"

int main(int argc, char** argv){
	unsigned int M,N,N_swp; char* ifn;
	if(argc<5){
		std::cout<<"\nusage: [exec_name] M N N_swp ifn";
	
	}
	else{
		M=atoi(argv[1]); N=atoi(argv[2]); N_swp=atoi(argv[3]); ifn=argv[4];


	std::fstream ifs(ifn,std::fstream::in);
	char ifs_buff[buff_sz]; char q;
	double **voltages=new double*[2*M*N+M+1];//define array on heap to avoid using all stack;
	for(unsigned int i=0;i<(2*M*N+M+1);i++){
		voltages[i]=new double[N_swp];
	
	}
	//double voltages[2*M*N+M+1][N_swp];
	if(ifs.is_open()){
		ifs.getline(ifs_buff,buff_sz);
		while(!ifs.eof()&&strstr(ifs_buff,"dc transfer curves")==NULL){//skip title untill node voltages is listed;
			ifs.getline(ifs_buff,buff_sz);

		}
		while(!ifs.eof()&&strstr(ifs_buff,"***** job concluded")==NULL){
			char* curr_ch=NULL;
			while((curr_ch=strstr(ifs_buff,"voltage"))==NULL){//skip lines until "voltage" header is reached;
				ifs.getline(ifs_buff,buff_sz);
				if(strstr(ifs_buff,"***** job concluded")!=NULL){goto fin_parse;} //if all nodes are parsed, close the file;

			}
			
			unsigned int n_clns=0;
			while(curr_ch!=NULL){//count how many volage columns in a row;
				n_clns++;
				curr_ch=strstr(curr_ch+strlen("voltage"),"voltage");
			}
			unsigned int* node_indxs=new unsigned int [n_clns];
			for(unsigned int i=0;i<n_clns;i++){//get node indexes;
				ifs>>node_indxs[i];
			}
			ifs.getline(ifs_buff,buff_sz);//read the rest of the line;
			
			
			char* strtk_ch; unsigned int curr_nd=0; double nv_curr=0; char mdf_tmp;
			for(unsigned int i=0;i<N_swp;i++){
				curr_nd=0;
				ifs.getline(ifs_buff,buff_sz);
				strtk_ch=strtok(ifs_buff," ");//ignoring first number; v_0 is provided at the beginning of the column;
				strtk_ch=strtok(NULL," ");
				while(strtk_ch!=NULL){
					mdf_tmp='\0';
					sscanf(strtk_ch,"%lf%c",&nv_curr,&mdf_tmp);
					if(mdf_tmp=='m'){//millivolts;
						nv_curr=nv_curr/1000.;
					
					}
					voltages[node_indxs[curr_nd]][i]=nv_curr;
					//std::cout<<" voltages["<<node_indxs[curr_nd]<<"]["<<i<<"]="<<voltages[node_indxs[curr_nd]][i];
					curr_nd++;
					strtk_ch=strtok(NULL," ");
				
				}
			}
			ifs.getline(ifs_buff,buff_sz);

		}
		
fin_parse:
		ifs.close();
		
		
		
		
		
		
		
		std::cout<<"\nHspcOut=[";
		for(unsigned int i=0;i<N_swp;i++){
			std::cout<<'[';
			for(unsigned int j=1;j<2*M*N;j++){//nodes 1 to 2*M*N;
				std::cout<<voltages[j][i]<<',';
				
			
			}
			if(i==N_swp-1){//i==N_swp; j==2*M*N;
				std::cout<<voltages[2*M*N][i]<<"]];";
			
			}
			else{
				std::cout<<voltages[2*M*N][i]<<"];";
			
			}
		
		}
		
		

	
	
	}
	else{
		std::cout<<"\nifs.is_open()";
	
	}
	for(unsigned int i=0;i<(2*M*N+M+1);i++){
		delete[] voltages[i];
	
	}
	delete[] voltages;
	}
	


	return 0;
}
