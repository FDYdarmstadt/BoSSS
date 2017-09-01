clear

xy = dlmread('RegLHS.txt');
z  = dlmread('RegRHS.txt');

a = xy\z

x = xy(:,1);
y = xy(:,2);

I = max(x) + 1;
J = max(y) + 1;


A = zeros(I,J);
for i=1:I
  for j=1:J
    A(i,j) = (i-1)*a(1) + (j-1)*a(2);
  end
end

for l=1:length(z)
  A(x(l) + 1,y(l) + 1) = z(l);  
end

surf(A);
