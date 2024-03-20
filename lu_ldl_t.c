































#include<stdio.h>
#include<math.h>


void show_arr(double* d_arr,unsigned int m){
	for(unsigned int i=0;i<m;i++){
		for(unsigned int j=0;j<m;j++){
			printf("%f ",*(d_arr+i*m+j));
		
		}
		printf("\n");
	
	
	
	}
		
		

}





unsigned int gen_lu(double* M,double* L,double* U,unsigned int m, double eps,unsigned int* P){//square matrixes M;
	//LU for square matrixes; Doolitle algorithm; partial pivotin;q
	//Input: M, m,esp: M is matrix to be decomposed (square m by m matrix); eps is pivotng sensitivity;
	//Inout: L,U,P: L and U are resulting decomposition matrixes, expected to be m by m zero initialized matrixes; P is a pivoting vector expected to be 1:m vector; L*U(:,P)==M(:,P);
	for(unsigned int i=0;i<m;i++){
		*(L+i*m+i)=1;//set ones on the diagonal;
	}
	for(unsigned int i=0;i<m;i++){
		*(U+i)=*(M+i);//first row of U is same as first row of M;
	}
	//arr[i][j] (i=1:m;j=:n) == *(arr+n*(i-1)+j)
	unsigned int pkk=0,pjj=0;
	double dot_pr=0.;
	for(int i=0;i<m;i++){
		
		if(fabs(*(U+i*m+*(P+i)-1))<eps){//pivoting;
		
		
			pkk=*(P+i)-1;
			double U_max_curr=fabs(*(U+i*m+pkk));
			unsigned int max_idx=i;
			for(int idx_curr=i;idx_curr<m;idx_curr++){
				if(fabs(*(U+i*m+*(P+idx_curr)-1))>U_max_curr){
					U_max_curr=fabs(*(U+i*m+*(P+idx_curr)-1));
					max_idx=idx_curr;
				
				}
			
			}
			*(P+i)=*(P+max_idx);
			*(P+max_idx)=pkk+1;
			
			
		
		}
		pkk=*(P+i)-1;
		for(unsigned int j=i+1;j<m;j++){
		
			dot_pr=0.;
			for(int idx=0;idx<=i-1;idx++){//should set idx to be of type int instead of unsigned int because (unsigned int)0<(-1)=1; '<=' is because i is run from 0 to m-1;
				dot_pr+=(*(L+j*m+idx))*(*(U+idx*m+pkk));
			
			}
			*(L+j*m+i)=(*(M+j*m+pkk)-dot_pr)/(*(U+i*m+pkk));
		
		}
		for(unsigned int j=i+1;j<m;j++){
			pjj=*(P+j)-1;
			dot_pr=0;
			for(unsigned int idx=0;idx<=i;idx++){//can keep idx as unsigned int because (unsigned int) 0<0=0;
			
			
			
				dot_pr+=(*(L+(i+1)*m+idx))*(*(U+idx*m+pjj));
			}
			*(U+(i+1)*m+pjj)=*(M+(i+1)*m+pjj)-dot_pr;
		
		}
	
	
	
	}
	return 0;




}
