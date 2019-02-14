function Imgout=drawkeypoints(Img, pointx,pointy,siz)
if length(pointx) ~=length(pointy)
    error('Length error!');
end
if max(size(Img)) < max(max([pointx,pointy]))
    error('Size error!');
end
imagesize=size(Img);
if numel(imagesize)==3
    Img=rgb2gray(Img);
end
[M,N]=size(Img);
Imgout=zeros(M,N,3);
Imgout(:,:,1)=Img;
Imgout(:,:,2)=Img;
Imgout(:,:,3)=Img;

for i=1:length(pointx)
    if pointx(i)>=siz && N-pointx(i)+1>=siz && pointy(i)>=siz && M-pointy(i)+1>=siz
        Imgout(pointy(i),pointx(i),:)=0;
        Imgout(pointy(i)-1,pointx(i),:)=0;
        Imgout(pointy(i)-2,pointx(i),:)=0;
        Imgout(pointy(i)+1,pointx(i),:)=0;
        Imgout(pointy(i)+2,pointx(i),:)=0;
        Imgout(pointy(i),pointx(i)-1,:)=0;
        Imgout(pointy(i),pointx(i)-2,:)=0;
        Imgout(pointy(i),pointx(i)+1,:)=0;
        Imgout(pointy(i),pointx(i)+2,:)=0;
    elseif pointx(i)>=siz && N-pointx(i)+1>=siz && pointy(i)>=siz && M-pointy(i)+1<siz
        Imgout(pointy(i):M,pointx(i),:)=0;
        Imgout(pointy(i)-1,pointx(i),:)=0;
        Imgout(pointy(i)-2,pointx(i),:)=0;
        Imgout(pointy(i),pointx(i)-1,:)=0;
        Imgout(pointy(i),pointx(i)-2,:)=0;
        Imgout(pointy(i),pointx(i)+1,:)=0;
        Imgout(pointy(i),pointx(i)+2,:)=0;
    elseif pointx(i)>=siz && N-pointx(i)+1<siz && pointy(i)>=siz && M-pointy(i)+1>=siz
        Imgout(pointy(i),pointx(i):N,:)=0;
        Imgout(pointy(i)-1,pointx(i),:)=0;
        Imgout(pointy(i)-2,pointx(i),:)=0;
        Imgout(pointy(i)+1,pointx(i),:)=0;
        Imgout(pointy(i)+2,pointx(i),:)=0;
        Imgout(pointy(i),pointx(i)-1,:)=0;
        Imgout(pointy(i),pointx(i)-2,:)=0;
    elseif pointx(i)<siz && N-pointx(i)+1>=siz && pointy(i)>=siz && M-pointy(i)+1>=siz
        Imgout(pointy(i),1:pointx(i),:)=0;
        Imgout(pointy(i)-1,pointx(i),:)=0;
        Imgout(pointy(i)-2,pointx(i),:)=0;
        Imgout(pointy(i)+1,pointx(i),:)=0;
        Imgout(pointy(i)+2,pointx(i),:)=0;
        Imgout(pointy(i),pointx(i)+1,:)=0;
        Imgout(pointy(i),pointx(i)+2,:)=0;
    elseif pointx(i)>=siz && N-pointx(i)+1>=siz && pointy(i)<siz && M-pointy(i)+1>=siz
        Imgout(1:pointy(i),pointx(i),:)=0;
        Imgout(pointy(i)+1,pointx(i),:)=0;
        Imgout(pointy(i)+2,pointx(i),:)=0;
        Imgout(pointy(i),pointx(i)-1,:)=0;
        Imgout(pointy(i),pointx(i)-2,:)=0;
        Imgout(pointy(i),pointx(i)+1,:)=0;
        Imgout(pointy(i),pointx(i)+2,:)=0;
    end
end
