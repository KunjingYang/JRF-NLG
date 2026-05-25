function SS = remark(S,a,c,b,d)
%[M,N] = size(S);
% t1 = b-a;
% t2 = d-c;
% E = zeros(time*t1,time*t2,3);
% for k=1:3
%     for i = 0:t1-1
%         for j = 0:t2-1
%             E(1+i*time:time*(i+1),1+time*j:time*(j+1),k)=S(a+i,c+j,k);
%         end
%     end
% end
% S(1:time*t1,1:time*t2,1:3)=E;
%% green
for i = a:1:b
    S(i,c-3,2:3)=0; S(i,c-3,1)=1e9;
    S(i,c-2,2:3)=0; S(i,c-2,1)=1e9;
    S(i,c-1,2:3)=0; S(i,c-1,1)=1e9;
    S(i,c,2:3)=0; S(i,c,1)=1e9;
    
    S(i,d,2:3)=0; S(i,d,1)=1e9;
    S(i,d+1,2:3)=0; S(i,d+1,1)=1e9;
    S(i,d+2,2:3)=0; S(i,d+2,1)=1e9;
    S(i,d+3,2:3)=0; S(i,d+3,1)=1e9;
end
for j = c:1:d
    S(a-3,j,2:3)=0; S(a-3,j,1)=1e9;
    S(a-2,j,2:3)=0; S(a-2,j,1)=1e9;
    S(a-1,j,2:3)=0; S(a-1,j,1)=1e9;
    S(a,j,2:3)=0; S(a,j,1)=1e9;
    S(b,j,2:3)=0; S(b,j,1)=1e9;
    S(b+1,j,2:3)=0; S(b+1,j,1)=1e9;
    S(b+2,j,2:3)=0; S(b+2,j,1)=1e9;
    S(b+3,j,2:3)=0; S(b+3,j,1)=1e9;
end
% %% blue
% for i = 1:time*t1
%     S(i,time*t2,1:2)=0; S(i,time*t2,3)=1e9;
%     S(i,time*t2+1,1:2)=0; S(i,time*t2+1,3)=1e9;
% end
% for j = 1:time*t2
%     S(time*t1,j,1:2)=0; S(time*t1,j,3)=1e9;
%     S(time*t1+1,j,1:2)=0; S(time*t1+1,j,3)=1e9;
% end

SS=S;
