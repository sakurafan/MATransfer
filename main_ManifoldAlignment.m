function main_ManifoldAlignment(DataCube,refImg)
refImg=double(refImg);
DataCube=double(DataCube);
[rows,columns,channels]=size(DataCube);
pixel=rows*columns;
DataCubeNorm=DataCube;
for i=1:channels
    temp=DataCubeNorm(:,:,i);
    DataCubeNorm(:,:,i)=DataCubeNorm(:,:,i)/max(temp(:));
end

%%% compute W for both images
refImgMatrix=reshape(refImg,[rows*columns,3]);
DataCubeMatrix=reshape(DataCube,[rows*columns,channels]);
Wref=calculateW(refImgMatrix);
Wspectral=calculateW(DataCubeMatrix);

%%%% Compute correspondence matrix  (pixel-wise matched by default)
input_i=1:10:rows;
input_j=1:10:columns;
[RowIndex,ColumnIndex]=meshgrid(input_i,input_j);
input_i=RowIndex(:);
input_j=ColumnIndex(:);
index=(input_j-1)*rows+input_i;
correspondence=sparse(index,index,ones(1,length(input_j)),pixel,pixel);

%%%% Call manifold alignment
outImage=manifold_alignment(refImg,DataCubeNorm,correspondence,Wref,Wspectral);