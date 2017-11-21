function finalx=manifold_alignment(refImg,DataCubeNorm,correspondence,Wref,Wspectral);

% M_source=imread('shiyiA.bmp');
% M_target=imread('shiyiB_nonlinear_power.bmp');
M_source=double(refImg);
%M_target=imread('shiyiB_nonlinear_serious.bmp');

M_target=double(DataCubeNorm);


source_pixel_rows=size(M_source,1);
source_pixel_cols=size(M_source,2);
target_pixel_rows=size(M_target,1);
target_pixel_cols=size(M_target,2);
channels=size(M_target,3);

x_source=reshape(M_source,[source_pixel_rows*source_pixel_cols,3]);
x_target=reshape(M_target,[target_pixel_rows*target_pixel_cols,channels]);
x_source=x_source';
x_target=x_target';

sample_source=size(x_source,2);
sample_target=size(x_target,2);

dimension1=size(x_source,1);
dimension2=size(x_target,1);



W_source=Wref;

W_target=Wspectral;

alpha=100;
C_target=correspondence*alpha;
C_source=C_target';

W=sparse([W_source,C_target;C_source,W_target]);

RowIndex=1:(sample_source+sample_target); 
Value=sum(W,2)';

D=sparse(RowIndex,RowIndex,Value,sample_source+sample_target,sample_source+sample_target);
% RD=D^(-1/2);
L=D-W;




Z=[x_source,sparse(size(x_source,1),size(x_target,2));sparse(size(x_target,1),size(x_source,2)),x_target];


EigL=Z*L*Z';
EigD=Z*D*Z';

k=4;

[vec,val]=eigs(EigL,EigD,k,'sm');

[~,idx] = sort(diag(val));
count=sum(diag(val)==0);
Vec=vec(:,idx(count+1:count+k));
% Vec=vec(:,idx(3));

F_source=Vec(1:dimension1,:);
F_target=Vec(dimension1+1:end,:);

Y_source=F_source'*x_source;
Y_target=F_target'*x_target;

Y_target_last=pinv(F_source')*Y_target;
y_target_last=zeros(target_pixel_rows,target_pixel_cols,k);
for i=1:k
    y_target_last(:,:,i)=reshape(Y_target_last(i,:),target_pixel_rows,target_pixel_cols);
end

finalx=uint8(y_target_last);


