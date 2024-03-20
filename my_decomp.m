                                                                                                                                                                                                                            































function varargout=my_decomp(M,type,eps)
varargout = cell(1,nargout-1);
    switch(type)
        case 'LU'
            disp('LU decomp');
            [L,U,P]=gen_lu(M,eps);
            varargout{1} = L;
            varargout{2} = U;
            varargout{3} = P;
            max(max(abs(M-L*U)))
%             max(max(abs(M-L1*U1)))
            


        case 'LDL^T'
            disp('LDL^T decomp');
        case 'LU_Dltl'
            disp('dltl')
            LU = LUDecompDoolittle(M);
            varargout{1} = LU;
        otherwise
            disp('unknown decomp');
    end
end

function [L,U,P]=gen_lu(M,eps)%LU for square matrices; Doolittle algorithm; partial pivoting;
dim_m=size(M,1);
L=eye(dim_m,dim_m); %diagonal entries of L is onesl
U=zeros(dim_m,dim_m);
%   L1=eye(dim_m,dim_m); %diagonal entries of L is onesl
%   U1=zeros(dim_m,dim_m);
%   U1(1,:)=M(1,:); %first row of U is same as first row of M;
U(1,:)=M(1,:); %first row of U is same as first row of M;
P=1:dim_m; %initially permetation is just identity;
%v_swp=zeros(dim_m,1);

for(kk=1:dim_m)
% for(jj=kk:dim_m)%prepare next row of U for next iteration;
%         U1(kk,jj)=M(kk,jj)-dot(L1(kk,1:kk-1),U1(1:kk-1,jj));
% end    
    if(abs(U(kk,P(kk)))<eps && (P(kk))<dim_m)%do partial pivoting (use columns);
        
        [~,U_max_idx]=max(abs(U(kk,P(kk):dim_m)));%find index of element with biggest absolute value in the reminder (including current element in case the rest of the raw is zeros) of the raw;
        P_kk_old=P(kk);%remember what should be swaped;
        P(P==(U_max_idx+P_kk_old-1))=P_kk_old;%find index of new column value and replace (swap) it with old column value;
        P(kk)=U_max_idx+P_kk_old-1;%replace old value with value of new column;
        %idx_upd=find(P==P_kk_old);
%         for(ii=kk+1:dim_m)
%             if(U(kk,ii)>eps)
%                 v_swp=U(:,kk); U(:,kk)=U(:,ii); U(:,ii)=v_swp;
%                 v_swp=M(:,kk); M(:,kk)=M(:,ii); M(:,ii)=v_swp;
%                 P(ii)=kk; P(kk)=ii;
%                 break;
%             end
%         end
    end
    for(jj=kk+1:dim_m)
        L(jj,kk)=(M(jj,P(kk))-dot(L(jj,1:kk-1),U(1:kk-1,P(kk))))/U(kk,P(kk));
%             L1(jj,kk)=(M(jj,kk)-dot(L1(jj,1:kk-1),U1(1:kk-1,kk)))/U1(kk,kk);
    end
    for(jj=kk+1:dim_m)%prepare next row of U for next iteration;
        U(kk+1,P(jj))=M(kk+1,P(jj))-dot(L(kk+1,1:kk),U(1:kk,P(jj)));
    end
end
                

end


function LU = LUDecompDoolittle(A)
    n = length(A);
    LU = A;
    % decomposition of matrix, Doolittle's Method
    for i = 1:1:n
        for j = 1:(i - 1)
            LU(i,j) = (LU(i,j) - LU(i,1:(j - 1))*LU(1:(j - 1),j)) / LU(j,j);
        end
        j = i:n;
        LU(i,j) = LU(i,j) - LU(i,1:(i - 1))*LU(1:(i - 1),j);
    end
    %LU = L+U-I
end


