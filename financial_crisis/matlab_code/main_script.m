prompt = {'Enter dimension:','Enter number of observations'};
answer = inputdlg(prompt);
d=str2double(answer{1}); nrows=str2double(answer{2});
%Filename
[C,b,c0]=getEllMat(array,d,60);
B=inv(C);
A=transpose(chol(B));  %get 
V=eye(d);
V(:,d+1)=zeros(d,1);
%V
%repmat(b,[1 d+1])
V=V-repmat(b,[1 d+1]);
%V
V=linsolve(A,V);
V=transpose(V);
F=[];
%V
for i=1:(d+1)
    R=V;
    R(i,:)=[];
    %i
    %R
    F=[F; linsolve(R,ones(d,1))'];
    %F
end
%F
P = polytope(F,ones(d+1,1));

xc = sdpvar(d,1); %center of Chebychev ball
r = sdpvar(1); %radius of Chebychev ball
xc1 = zeros(d,1); %center of input ball
r1 = c0; %radius of input ball

Con = [F*xc+r*sqrt(sum(F.^2,2)) <= ones(d+1,1), norm(xc-xc1,2) <= r1-r];
sol = solvesdp(Con,-r);

xc=value(xc);
%value(r);

eig_min=eig(C);
%eig_min
radius=min(eig_min);
center=A*xc+b;

fileID2 = fopen('mat_ellips_ve.txt','w');
%fileID2 = fopen('mat_ellips3.txt','w');
for i=1:d
    for j=1:d
        fprintf(fileID2,'%f',C(i,j));
        fprintf(fileID2,'%s',' ');
    end
end
for j=1:d
    fprintf(fileID2,'%f',b(j));
    fprintf(fileID2,'%s',' ');
end
fprintf(fileID2,'%f',c0);
fclose(fileID2);
%fclose(fileID2);


fileID = fopen('center_radius_ve.txt','w');
for i=1:(d+1)
    if i<(d+1)
        fprintf(fileID,'%f',center(i));
        fprintf(fileID,'%s',' ');
    else
        fprintf(fileID,'%f',radius);
    end
end     
fclose(fileID);

movefile center_radius_ve.txt ~/fullpathto/volume_approximation/financial_crisis/project_main
movefile mat_ellips_ve.txt ~/fullpathto/volume_approximation/financial_crisis/project_main


cd ~/fullpathto/volume_approximation/financial_crisis/project_main
[status,cmdout] = system('-ve 1 2 [dim] [k] [W] [e] mat_ellips_ve.txt center_radius_ve.txt');
cmdout
cd ~/fullpathto/volume_approximation/financial_crisis/matlab_code
