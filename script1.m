
A =load('output_1000.dat');
q=A(:,4);
C =zeros(100,100,100);


for i=1:100
   for j=1:100
      for k=1:100
          C(i,j,k)=q((k)+10*(j+10*i));
       end
   end
end
[x,y,z]=meshgrid(1:100,1:100,1:100);
scatter3(x(:),y(:),z(:),5,C(:),'filled')
