function num = match(image1, image2)
[im1, des1, loc1] = sift(image1);
[im2, des2, loc2] = sift(image2);
distRatio = 0.6;   
des2t = des2';                          
for i = 1 : size(des1,1)
   dotprods = des1(i,:) * des2t;        
   [vals,indx] = sort(acos(dotprods));  
   if (vals(1) < distRatio * vals(2))
      match(i) = indx(1);
   else
      match(i) = 0;
   end
end
im3 = appendimages(im1,im2);
figure('Position', [100 100 size(im3,2) size(im3,1)]);
colormap('gray');
imagesc(im3);
hold on;
cols1 = size(im1,2);
for i = 1: size(des1,1)
  if (match(i) > 0)
    line([loc1(i,2) loc2(match(i),2)+cols1], ...
         [loc1(i,1) loc2(match(i),1)], 'Color', 'c');
  end
end
hold off;
num = sum(match > 0);
fprintf('Found %d matches.\n', num);




