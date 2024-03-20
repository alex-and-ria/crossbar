































								
#include<iostream>//cout
#include<fstream>//std::fstream;

//#define M 3
//#define N 3

//#define N_swp 2

#define R_wl 5
#define R_bl 5
#define R_min 1000
#define R_max 100000
#define V_min 1
#define V_max 100

#define buff_sz 256
#define gen_ofn(pchar_buf,m,n,n_swp) sprintf(pchar_buf,"./hspice/%dx%d_%d_dc_sweep_gen",m,n,n_swp);

int main(int argc, char** argv){
	unsigned int M,N,N_swp;
	if(argc<4){
		std::cout<<"\nusage: [exec_name] M N N_swp";
	
	}
	else{
		M=atoi(argv[1]); N=atoi(argv[2]); N_swp=atoi(argv[3]);
	unsigned int** Rstrs=new unsigned int*[M*N];//can be large array, so allocate memory for it in heap;
	for(unsigned int i=0;i<M*N;i++){
		Rstrs[i]=new unsigned int[M*N];
	
	}
	int** Vin=new int*[N_swp];//can be large array, so allocate memory for it in heap;
	for(unsigned int i=0;i<N_swp;i++){
		Vin[i]=new int[M];
	
	}
	std::fstream ofs;
	char ofn_buff[buff_sz], data_buf[buff_sz];
	srand(time(NULL));
	gen_ofn(ofn_buff,M,N,N_swp)
	ofs.open(ofn_buff,std::ios_base::out);
	if(ofs.is_open()){
		sprintf(data_buf,"\n*%dx%d cd seep hspice file generated by ../hpice_gen_dc.cpp",M,N);
		ofs<<data_buf;
		for(unsigned int i=0;i<M;i++){
			for(unsigned int j=0;j<N;j++){
				Rstrs[i][j]=rand()%(R_max-R_min+1)+R_min;
				sprintf(data_buf,"\nR%d_%d %d %d %u",i,j,2*i*N+2*j+1,2*i*N+2*j+2,Rstrs[i][j]);
				ofs<<data_buf;
				if(j==0){
					sprintf(data_buf,"\nRwl%d_%d %d %d %d",i,j,2*M*N+i+1,2*i*N+2*j+1, R_wl);
					ofs<<data_buf;
				
				}
				else{
					sprintf(data_buf,"\nRwl%d_%d %d %d %d",i,j,2*i*N+2*(j-1)+1,2*i*N+2*j+1, R_wl);
					ofs<<data_buf;
				
				}
				if(i==M-1){
					sprintf(data_buf,"\nRbl%d_%d %d %d %d",i,j,2*i*N+2*j+2,0, R_bl);
					ofs<<data_buf;
				
				}
				else{
					sprintf(data_buf,"\nRbl%d_%d %d %d %d",i,j,2*i*N+2*j+2,2*(i+1)*N+2*j+2, R_bl);
					ofs<<data_buf;
				
				}
				
				
			
			}
		
		}
		for(unsigned int i=0;i<M;i++){
			sprintf(data_buf,"\nV%d %d 0 V_%d",i,2*M*N+i+1,i);
			ofs<<data_buf;
		
		}
		ofs<<"\n\n*analysis";
		ofs<<"\n.DC DATA=voltages";
		ofs<<"\n.DATA voltages\n\t";
		for(unsigned int i=0;i<M;i++){
			sprintf(data_buf,"V_%d\t",i);
			ofs<<data_buf;
		
		}
		for(unsigned int i=0;i<N_swp;i++){
			ofs<<'\n';
			for(unsigned int j=0;j<M;j++){
				Vin[i][j]=(rand()%(V_max-V_min+1)+V_min);
				ofs<<'\t'<<Vin[i][j];
			
			}
		
		}
		ofs<<"\n.ENDDATA";
		ofs<<"\n.option post=2 POST_VERSION=2001";//set otput .sw0 in ASCII format;
		ofs<<"\n.probe DC V(*)";//save probe results to .sw0 file;
//		ofs<<"\n.print DC V(*)";
		ofs<<"\n.END";
		
		
		
		ofs<<"\n33";
		ofs.close();
	
	}
	else{
		std::cout<<"\nofn.is_open()";
	
	}
	
	
	std::cout<<"\nCnds=[";
	for(unsigned int i=0;i<M;i++){
		std::cout<<'[';
		for(unsigned int j=0;j<N-1;j++){
			std::cout<<1./(Rstrs[i][j]+0.)<<',';
		
		}
		if(i==M-1){//i==M-1; j==N-1;
			std::cout<<1./(Rstrs[i][N-1]+0.)<<"]];";
		
		}
		else{
			std::cout<<1./(Rstrs[i][N-1]+0.)<<"];";
		
		}
	}
	std::cout<<" Gwl="<<1./(R_wl+0.)<<"; Gbl="<<1./(R_bl+0.)<<';';
	
	std::cout<<" Vin=[";
	for(unsigned int i=0;i<N_swp;i++){
		std::cout<<'[';
		for(unsigned int j=0;j<M-1;j++){
			std::cout<<Vin[i][j]<<',';
		
		}
		if(i==N_swp-1){//i==N_swp-1; j=M-1;
			std::cout<<Vin[i][M-1]<<"]];";
		
		}
		else{
			std::cout<<Vin[i][M-1]<<"];";
		
		}
	
	}
	std::cout<<'\n'<<ofn_buff;
	
	for(unsigned int i=0;i<M*N;i++){//free memory allocated in heap to avoind memory leaks;
		delete[] Rstrs[i];
	
	}
	delete[] Rstrs;
	for(unsigned int i=0;i<N_swp;i++){
		delete[] Vin[i];
	
	}
	delete[] Vin;
	
	}
	
	return 0;







}								
