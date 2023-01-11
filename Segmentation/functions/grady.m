function M=grady(I)

%Calcul le gradient en y d'une image I
%Syntaxe: grady(I)

[m,n]=size(I);
M=zeros(m,n);

M(:,1:end-1)=I(:,2:end)-I(:,1:end-1);


