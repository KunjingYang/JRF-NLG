function  ouput_image = barrel_distortion (A, k1)
    image = A; % ¶ÁÈ¡ÄúµÄÍ¼Ïñ

%focal_length = 100;
% k1 = 0.000001;

[h,w,L] = size(image);
cx = (w+1)/2;
cy = (h+1)/2;
[x,y] = meshgrid(1:w,1:h);
x = x-cx;
y = y-cy;
r2 = x.^2+y.^2;
x = x.*(1+k1*r2);
y = y.*(1+k1*r2);
x = x+cy;
y = y+cy;

ouput_image = zeros(size(image));
for i = 1:L
    ouput_image(:,:,i)=interp2(double(image(:,:,i)),x,y);
end