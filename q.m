




















%dim=16; %M=N=dim;%Cnds=sym('G',dim); Gwl=sym('Gwl'); Gbl=sym('Gbl');
%Ivec=sym(zeros(1,2*dim*dim));
%Vin=sym('Vin',[1,dim]);
%Gm=sym(zeros(2*dim*dim));
fl_sprs=1;

% m=3; n=3;
% Cnds=[[1.41445e-05,2.74296e-05,2.10009e-05];[9.58681e-05,1.17854e-05,0.000392311];
%     [1.99132e-05,1.96282e-05,1.04814e-05]]; Gwl=0.2; Gbl=0.2;
%  Vin=[[50,31,56];[21,18,7]];
% HspcOut=[[49.9844,0.0457723,49.9723,0.0296753,49.9671,0.138938,30.9232,0.0422405,30.8613,0.0228258,30.8011,0.133706,55.986,0.0239062,55.9776,0.0141591,55.9747,0.0683178];[20.9934,0.0223386,20.9884,0.011428,20.9862,0.0768214,17.9554,0.0208555,17.9194,0.008551,17.8845,0.0746259,6.9983,0.0107756,6.9972,0.0046186,6.9968,0.0374953]];
N_swp=size(HspcOut,1);
Ivec=zeros(size(HspcOut));
if (fl_sprs==0)
    Gm=zeros(2*m*n);
else
    Gm=sparse(1,1,0,2*m*n,2*m*n,m*(6+(n-2)*4)+n*(6+(m-2)*4));
end



for ii=1:m
    for jj=1:n
        if(jj==1)
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+1)=2*Gwl+Cnds(ii,jj);
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+3)=-Gwl;
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+2)=-Cnds(ii,jj);
            
            

            for(kk=1:N_swp)
                Ivec(kk,2*(ii-1)*n+2*(jj-1)+1)=Gwl*Vin(kk,ii);
            end

            %Ivec1(2*(ii-1)*n+2*(jj-1)+1)=Gwl*Vin(1,ii);
            %Ivec2(2*(ii-1)*n+2*(jj-1)+1)=Gwl*Vin(2,ii);
        elseif(jj==n)
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+2*(jj-1)+1)=Gwl+Cnds(ii,jj);
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+2*(jj-2)+1)=-Gwl;
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+2*(jj-1)+2)=-Cnds(ii,jj);
        else
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+2*(jj-1)+1)...
            =2*Gwl+Cnds(ii,jj);
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+2*(jj-2)+1)=-Gwl;
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+2*(jj)+1)=-Gwl;
            Gm(2*(ii-1)*n+2*(jj-1)+1,2*n*(ii-1)+2*(jj-1)+2)=-Cnds(ii,jj);
        end
        if(ii==1)
            Gm(2*(ii-1)*n+2*jj,2*(jj-1)+2)=Gbl+Cnds(ii,jj);
            Gm(2*(ii-1)*n+2*jj,2*n+2*(jj-1)+2)=-Gbl;
            Gm(2*(ii-1)*n+2*jj,2*(jj-1)+1)=-Cnds(ii,jj);
        elseif(ii==m)
            Gm(2*(ii-1)*n+2*jj,2*n*(ii-1)+2*(jj-1)+2)=2*Gbl+Cnds(ii,jj);
            Gm(2*(ii-1)*n+2*jj,2*n*(ii-2)+2*(jj-1)+2)=-Gbl;
            Gm(2*(ii-1)*n+2*jj,2*n*(ii-1)+2*(jj-1)+1)=-Cnds(ii,jj);
        else
            Gm(2*(ii-1)*n+2*jj,2*n*(ii-1)+2*(jj-1)+2)=2*Gbl+Cnds(ii,jj);
            Gm(2*(ii-1)*n+2*jj,2*n*(ii)+2*(jj-1)+2)=-Gbl;
            Gm(2*(ii-1)*n+2*jj,2*n*(ii-2)+2*(jj-1)+2)=-Gbl;
            Gm(2*(ii-1)*n+2*jj,2*n*(ii-1)+2*(jj-1)+1)=-Cnds(ii,jj);
        end
        
    end
end
% P_lu=1:2*m*n;
% for(ii=1:m)
%     for(jj=1:n)
%         %val_tmp=P_lu(ii+1);
%         %P_lu(ii+1)=P_lu()
% 
%     end
% end

%figure(1); spy(Gm);
% Gpm=zeros(m,2*n,2*m*n);
% for k=1:m
%     for ii=1:2*n
%         for jj=1:2*m*n
%             if((ii+(k-1)*2*n)==jj)
%                 Gpm(k,ii,jj)=1;
%             end
%         end
%     end
%  end
% Pm=zeros(2*m*n);
% for ii=1:m
%     for jj=1:n
%         Pm(2*m*(jj-1)+ii,:)=permute(Gpm(ii,2*(jj-1)+1,:),[3,1,2]);
%         Pm(2*m*(jj-1)+m+ii,:)=permute(Gpm(ii,2*(jj-1)+2,:),[3,1,2]);
%     end
%  end
% tic(); V_nd=linsolve(Gm,Ivec'); toc();
% [arr_P,arr_R, arr_C]=equilibrate(Gm);
% tic(); V_nd1=cgs(arr_R*arr_P*Gm,arr_R*arr_P*(Ivec')); toc();
% fh_gm=@(x) (Gm*x);
% tic(); V_nd2=newtonraphson(fh_gm,zeros(1,2*dim*dim)'); toc();

%tic(); y= Gm\(Ivec'); toc();
%times_v=zeros(1,4);
%fprintf("m=%d; n=%d; N_swp=%d;\n",m,n,N_swp);
lu_decomp_t=0; lu_solv_t=0; lu_max_abs_err=0;
%disp('lu'); %tic();
tic();[L,U,P] = lu(Gm); lu_decomp_t=toc();
%disp('runs'); times_lu=0;
for(kk=1:N_swp)
    tic(); lu_y = L\(P*Ivec(kk,:)');
    lu_x = U\lu_y; lu_solv_t=lu_solv_t+toc();
    if(lu_max_abs_err<max(abs(lu_x-HspcOut(kk,:)')))
        lu_max_abs_err=max(abs(lu_x-HspcOut(kk,:)'));
    end
end
%disp(times_lu);
%fprintf("lu_decomp_t=%g, lu_solv_t=%g, lu_max_abs_err=%g\n", lu_decomp_t,lu_solv_t,lu_max_abs_err);
fprintf("%f %f ", lu_decomp_t,lu_solv_t);


% disp('ldl');
%     tic(); dA = decomposition(Gm,'ldl','upper'); toc();
% disp('runs'); times_ldl=0; max_abs_err=0.;
% for(kk=1:N_swp)
%     %disp(kk);
%     tic(); ldl_x_curr=dA\(Ivec(kk,:)'); times_ldl=times_ldl+toc();
%     if(max_abs_err<max(abs(ldl_x_curr'-HspcOut(kk,:))))
%         max_abs_err=max(abs(ldl_x_curr'-HspcOut(kk,:)));
%     end
%     %figure(kk); plot(1:2*m*n,abs(ldl_x_curr'-HspcOut(kk,:)));
% 
% end
% disp(times_ldl);
% fprintf("max_abs_err=%g\n", max_abs_err);
if(2*m*n<100)
lib_lu_decomp_t=0; lib_lu_solv_t=0; lib_lu_max_abs_err=0;
tst_L=zeros(2*m*n,2*m*n); tst_U=zeros(2*m*n,2*m*n); tst_P=1:2*m*n;
ptst_M=libpointer('doublePtr',full(Gm)'); ptst_L=libpointer('doublePtr',tst_L); ptst_U=libpointer('doublePtr',tst_U); ptst_P=libpointer('uint32Ptr',tst_P');
loadlibrary('lu_ldl_t.so','lu_ldl_t.h');
tic();
    calllib('lu_ldl_t','gen_lu',ptst_M,ptst_L,ptst_U,2*m*n,0.00001,ptst_P);
lib_lu_decomp_t=toc();
L0=ptst_L.Value';U0=ptst_U.Value'; P0=ptst_P.Value';
unloadlibrary('lu_ldl_t');
for(kk=1:N_swp)
    tic(); lib_lu_y = L0\(Ivec(kk,P0)');
    lib_lu_x = U0(:,P0)\lib_lu_y(P0); lib_lu_solv_t=lib_lu_solv_t+toc();
    if(lib_lu_max_abs_err<max(abs(lib_lu_x-HspcOut(kk,P0)')))
        lib_lu_max_abs_err=max(abs(lib_lu_x-HspcOut(kk,P0)'));
    end
end
%disp(times_lu);
%fprintf("lib_lu_decomp_t=%g, lib_lu_solv_t=%g, lib_lu_max_abs_err=%g\n", lib_lu_decomp_t,lib_lu_solv_t,lib_lu_max_abs_err);
fprintf("%f %f\n", lib_lu_decomp_t,lib_lu_solv_t);

else
fprintf("%f %f\n", 0,0);
end

%tic(); ldl_x2=dA\Ivec(2)'; toc();
%figure(2); plot(1:2*m*n,abs(ldl_x2'-HspcOut(2,:)));

%disp('4'); tic();  gmres(Gm,Ivec'); toc();
%tic(); disp('prec cg'); L = ichol(Gm,struct('michol','on')); [x1,fl1,rr1,it1,rv1] = pcg(Gm,Ivec',1e-6,100,L,L'); toc();
%fh_gm=@(x) (Gm*x);
%disp('NR'); tic(); V_nd2=newtonraphson(fh_gm,zeros(1,2*m*n)'); toc();
%figure(2); plot(HspcOut'-Pm*lu_x);
%waitforbuttonpress()



% 
% 
% 
% 
% 
% 
% 
% function [sol, error] = newtonraphson(F, x0, maxiter)
% % Newton-Raphson solver for a system of nonlinear equations, where F is the
% % handle to the system of equations and x0 is the starting point for the
% % solver. The output of F is a column vector, x0 is a column vector.
% % Maxiter species the maximum number of iterations, default is 1000.
% % initialise values
% error1 = 1e8;
% x = x0;
% iter = 0;
% if nargin < 3
%     maxiter = 1e3;
% end
% error = zeros(maxiter, 1);
% sol = F(x0); % compute test solution
% nvar = length(sol); % compute size of the system
% h = 1e-4 .* ones(nvar, 1); % size of step for derivative computation
% while error1 > 1e-12
%     iter = iter+1; % update iteration
%     f = F(x); % evaluate function at current point
%     
%     J = jacobiannum(F, x, h); % compute the Jacobian
%     y = -J\f;  % solve the linear equations
%     x = x + y; % move the solution
%     
%     % calculate errors
%     error1 = sqrt(sum(y.^2));
%     error(iter) = sqrt(sum(f.^2));
%     
%     % break computation if maximum number of iterations is exceeded
%     if iter == maxiter
%         warning('Maximum number of iterations (%i) exceeded, solver may not have converged.', maxiter);
%         break;
%     end 
% end
% % return solution and error for each iteration
% sol = x;
% error = error(1:iter);
% end
% function J = jacobiannum(F,x,h)
% % Computes the Jacobian matrix of the function F at the point x, where h is
% % the step size to take on each dimension (has to be small enough for
% % precision). Note that the result of f and vectors x and h must be column
% % vectors of the same dimensions.
% J = (F(repmat(x,size(x'))+diag(h))-F(repmat(x,size(x'))))./h';
% end