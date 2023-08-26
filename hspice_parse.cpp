































#include<iostream>//cout;
#include<fstream>//ifstream;
#include<string.h>//strtok;

//#define ifn "./hspice/3x3_crossbar_gen.ic0"
//#define M 3
//#define N 3
#define ibuff_size 256
#define debug 0
	

int main(int argc, char** argv){
	char* ifn;
	unsigned int M,N;
	if(argc!=4){
		std::cout<<"\nusage: parse_obj input_file_path M N";	
	}
	else{
		ifn=argv[1];
		M=atoi(argv[2]);
		N=atoi(argv[3]);
	
	}
	std::ifstream ifh(ifn,std::ios::in|std::ios::binary);
	unsigned char ibuff[ibuff_size]; char get_nl; unsigned char* uc_arr_tkn;
	bool strt_fl=0; unsigned int num_nodes=2*M*N+M;//2*m*n nodes for crossbar itself, m nodes for input sources + ground node (which is 0 V);
	double node_voltages[num_nodes];
	if(ifh.is_open()){
		while(!ifh.eof()){
			ifh.get((char*)ibuff,ibuff_size,'\n');//get line of size at most ibuff_size, or until '\n' symbol;
			//std::cout<<'\n'<<"line:"<<ibuff;
			ifh.get(get_nl);//read '\n' character from stream to allow getting new line on next iteration;
			if(strt_fl==0 && strstr((char*)ibuff,".nodeset")!=NULL){//skip title and comment strings;
				strt_fl=1; if(debug) std::cout<<"\nhere!";
			}
			else if(strt_fl==1 && strstr((char*)ibuff,"END")==NULL){
				if(debug) std::cout<<'\n'<<"line:"<<ibuff;
				uc_arr_tkn=(unsigned char*)strtok((char*)ibuff,"+= ");
				//std::cout<<'\n'<<"\nf_token="<<uc_arr_tkn<<" ("<<atoi(uc_arr_tkn)<<")";
				unsigned int node_idx=atoi((char*)uc_arr_tkn);
				if(debug) std::cout<<"\nnode_idx="<<node_idx;
				uc_arr_tkn=(unsigned char*)strtok(NULL,"+= ");
				double d_scl=1.;
				if(uc_arr_tkn[strlen((char*)uc_arr_tkn)-1]=='m') d_scl=0.001;//milli volts;
				else if(uc_arr_tkn[strlen((char*)uc_arr_tkn)-1]=='u') d_scl=0.000001;//micro volts;
				node_voltages[node_idx-1]=strtof((char*)uc_arr_tkn,NULL)*d_scl;//convert to float,
				
				/*while(uc_arr_tkn!=NULL){
					
					std::cout<<"\ntoken="<<uc_arr_tkn;
					uc_arr_tkn=(unsigned char*)strtok(NULL,"+= ");
					//std::cout<<"\ntoken="<<uc_arr_tkn;
				
				}*/
				
			
			}
			else if(strstr((char*)ibuff,"END")!=NULL){//"END" string; if this condition not specified, Segmentation fault arises;
				strt_fl=0;
				//std::cout<<"\nelse"; ifh.close(); break;
			}
		
		}
		
		
		ifh.close();
	
	}
	else{
		std::cout<<"\nifh.is_open()";
	
	}
	for(unsigned int i=0;i<num_nodes;i++) if(debug) std::cout<<"\nnode["<<i+1<<"]="<<node_voltages[i];
	
	
	
	
	//cout node voltages for Mathematica; skip first M nodes, because it is input (voltage) sources nodes;
	std::cout<<"\nHspcOut=";
	std::cout<<"{"<<node_voltages[M]; for(unsigned int i=M+1;i<num_nodes;i++) std::cout<<','<<node_voltages[i]; std::cout<<"};\n";

	return 0;
}
