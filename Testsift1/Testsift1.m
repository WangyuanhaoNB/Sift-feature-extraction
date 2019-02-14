%%% Test sift Dog
clear;close all;
Img=imread('Macintosh HD/user/george/picture/pap.er/D5NmSl8IQig.jpg');
img=double(rgb2gray(Img));
[M,N]=size(img);
S=3;
sigma=0.5;thresh=0.03;gama=10;  % param
Ht=cell(S+3,1);
Gauss=cell(S+3,1);
Dog=cell(S+2,1);
for i=1:S+3
    Ht{i}=fspecial('gaussian', 6, sigma*2^((i-1)/S));
    Gauss{i}=imfilter(img,Ht{i},'replicate');
end
for i=1:S+2
    Dog{i}=Gauss{i+1}-Gauss{i};
end
%%  recursion

Improc=cell(S+2,1);   %   variable process
for i=1:S+2
    Improc{i}=zeros(M+2,N+2);
    Improc{i}(2:M+1,2:N+1)=Dog{i};
    Improc{i}(1:M+2,1)=Improc{i}(1:M+2,2);
    Improc{i}(1:M+2,end)=Improc{i}(1:M+2,end-1);
    Improc{i}(1,2:N+1)=Improc{i}(2,2:N+1);
    Improc{i}(end,2:N+1)=Improc{i}(end-1,2:N+1);
end
indic=zeros(M,N,S);
for i=2:S+1
    for j=2:M+1
        for k=2:N+1
            if Improc{i}(j,k)>Improc{i}(j+1,k) && Improc{i}(j,k)>Improc{i}(j-1,k) && Improc{i}(j,k)>Improc{i}(j,k+1) && Improc{i}(j,k)>Improc{i}(j,k-1) && Improc{i}(j,k)>Improc{i}(j+1,k+1) && Improc{i}(j,k)>Improc{i}(j-1,k-1) && Improc{i}(j,k)>Improc{i}(j+1,k-1) && Improc{i}(j,k)>Improc{i}(j-1,k+1)...
                && Improc{i}(j,k)>Improc{i+1}(j,k) && Improc{i}(j,k)>Improc{i+1}(j+1,k) && Improc{i}(j,k)>Improc{i+1}(j-1,k) && Improc{i}(j,k)>Improc{i+1}(j,k+1) && Improc{i}(j,k)>Improc{i+1}(j,k-1) && Improc{i}(j,k)>Improc{i+1}(j+1,k+1) && Improc{i}(j,k)>Improc{i+1}(j-1,k-1) && Improc{i}(j,k)>Improc{i+1}(j+1,k-1) && Improc{i}(j,k)>Improc{i+1}(j-1,k+1)...
                && Improc{i}(j,k)>Improc{i-1}(j,k) && Improc{i}(j,k)>Improc{i-1}(j+1,k) && Improc{i}(j,k)>Improc{i-1}(j-1,k) && Improc{i}(j,k)>Improc{i-1}(j,k+1) && Improc{i}(j,k)>Improc{i-1}(j,k-1) && Improc{i}(j,k)>Improc{i-1}(j+1,k+1) && Improc{i}(j,k)>Improc{i-1}(j-1,k-1) && Improc{i}(j,k)>Improc{i-1}(j+1,k-1) && Improc{i}(j,k)>Improc{i-1}(j-1,k+1)
            indic(j-1,k-1,i-1)=1;
            end
            if Improc{i}(j,k)<Improc{i}(j+1,k) && Improc{i}(j,k)<Improc{i}(j-1,k) && Improc{i}(j,k)<Improc{i}(j,k+1) && Improc{i}(j,k)<Improc{i}(j,k-1) && Improc{i}(j,k)<Improc{i}(j+1,k+1) && Improc{i}(j,k)<Improc{i}(j-1,k-1) && Improc{i}(j,k)<Improc{i}(j+1,k-1) && Improc{i}(j,k)<Improc{i}(j-1,k+1)...
                && Improc{i}(j,k)<Improc{i+1}(j,k) && Improc{i}(j,k)<Improc{i+1}(j+1,k) && Improc{i}(j,k)<Improc{i+1}(j-1,k) && Improc{i}(j,k)<Improc{i+1}(j,k+1) && Improc{i}(j,k)<Improc{i+1}(j,k-1) && Improc{i}(j,k)<Improc{i+1}(j+1,k+1) && Improc{i}(j,k)<Improc{i+1}(j-1,k-1) && Improc{i}(j,k)<Improc{i+1}(j+1,k-1) && Improc{i}(j,k)<Improc{i+1}(j-1,k+1)...
                && Improc{i}(j,k)<Improc{i-1}(j,k) && Improc{i}(j,k)<Improc{i-1}(j+1,k) && Improc{i}(j,k)<Improc{i-1}(j-1,k) && Improc{i}(j,k)<Improc{i-1}(j,k+1) && Improc{i}(j,k)<Improc{i-1}(j,k-1) && Improc{i}(j,k)<Improc{i-1}(j+1,k+1) && Improc{i}(j,k)<Improc{i-1}(j-1,k-1) && Improc{i}(j,k)<Improc{i-1}(j+1,k-1) && Improc{i}(j,k)<Improc{i-1}(j-1,k+1)
            indic(j-1,k-1,i-1)=1;
            end
        end
    end
end

hy=[-0.5;0;0.5];
hx=[-0.5,0,0.5];
fx=zeros(M+2,N+2,S);
fy=zeros(M+2,N+2,S);
fxx=zeros(M+2,N+2,S);
fyy=zeros(M+2,N+2,S);
fxy=zeros(M+2,N+2,S);
fyx=zeros(M+2,N+2,S);
fz=zeros(M+2,N+2,S);
fzy=zeros(M+2,N+2,S);
fzx=zeros(M+2,N+2,S);
fzz=zeros(M+2,N+2,S);
fxz=zeros(M+2,N+2,S);
fyz=zeros(M+2,N+2,S);
for i=1:S
    fx(:,:,i)=imfilter(Improc{i+1},hx,'replicate');
    fy(:,:,i)=imfilter(Improc{i+1},hy,'replicate');
    fxx(:,:,i)=imfilter(fx(:,:,i),hx,'replicate');
    fxy(:,:,i)=imfilter(fx(:,:,i),hy,'replicate');
    fyx(:,:,i)=imfilter(fy(:,:,i),hx,'replicate');
    fyy(:,:,i)=imfilter(fy(:,:,i),hy,'replicate');
    fz(:,:,i)=(Improc{i+2}-Improc{i})/2;
    fzy(:,:,i)=imfilter(fz(:,:,i),hy,'replicate');
    fzx(:,:,i)=imfilter(fz(:,:,i),hx,'replicate');
    if i==1
        fzz(:,:,i)=(fz(:,:,i+1)-Improc{i+1}/2)/2;
        fxz(:,:,i)=(fx(:,:,i+1)-imfilter(Improc{i},hx,'replicate'))/2;
        fyz(:,:,i)=(fy(:,:,i+1)-imfilter(Improc{i},hy,'replicate'))/2;
    elseif i==S
        fzz(:,:,i)=(-Improc{i+1}/2-fz(:,:,i-1))/2;
        fxz(:,:,i)=(imfilter(Improc{i+2},hx,'replicate')-fx(:,:,i-1))/2;
        fyz(:,:,i)=(imfilter(Improc{i+2},hy,'replicate')-fy(:,:,i-1))/2;
    else
        fzz(:,:,i)=(fz(:,:,i+1)-fz(:,:,i-1))/2;
        fxz(:,:,i)=(fx(:,:,i+1)-fx(:,:,i-1))/2;
        fyz(:,:,i)=(fy(:,:,i+1)-fy(:,:,i-1))/2;
    end
end

Indic=zeros(size(indic,1),size(indic,2));
for i=1:S
    Ind=find(indic(:,:,i)==1);
    [I,J]=ind2sub(size(indic(:,:,i)),Ind);
    D_x=zeros(length(Ind),1);
    for j=1:length(Ind)
        F=[fxx(I(j)+1,J(j)+1,i),fxy(I(j)+1,J(j)+1,i),fxz(I(j)+1,J(j)+1,i);
            fyx(I(j)+1,J(j)+1,i),fyy(I(j)+1,J(j)+1,i),fyz(I(j)+1,J(j)+1,i);
            fzx(I(j)+1,J(j)+1,i),fzy(I(j)+1,J(j)+1,i),fzz(I(j)+1,J(j)+1,i)];
        Finv=inv(F);
        x_x=-Finv*[fx(I(j)+1,J(j)+1,i);fy(I(j)+1,J(j)+1,i);fz(I(j)+1,J(j)+1,i)];
        D_x(j)=Improc{i+1}(I(j)+1,J(j)+1)+0.5*[fx(I(j)+1,J(j)+1,i),fy(I(j)+1,J(j)+1,i),fz(I(j)+1,J(j)+1,i)]*x_x;
        if max(x_x)>0.5 || D_x(j)/max(max(img))<thresh
            indic(I(j),J(j),i)=0;
        end
    end
    Indic=Indic | indic(:,:,i);
end
disp(sum(sum(Indic)));
%   fxx¡¢fxy¡¢fyy
for i=1:S
    Ind=find(indic(:,:,i)==1);
    [I,J]=ind2sub(size(indic(:,:,i)),Ind);
    for j=1:length(Ind)
        H=[fxx(I(j)+1,J(j)+1,i),fxy(I(j)+1,J(j)+1,i);fxy(I(j)+1,J(j)+1,i),fyy(I(j)+1,J(j)+1,i)];
        Tr=trace(H);
        Det=det(H);
        if Tr^2/Det >= (gama+1)^2/gama
            indic(I(j),J(j),i)=0;
        end
    end
    Indic=Indic | indic(:,:,i);
end
disp(sum(sum(Indic)));

imgshow=drawkeypoints(img,J+1,I+1,15);
imgshow=uint8(imgshow);
imshow(imgshow);
