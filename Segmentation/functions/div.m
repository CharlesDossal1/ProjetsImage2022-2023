function M=div(px,py)

%calcul de la divergence d'un vecteur (px,py)
%px et py de meme taille
%(de sorte que div=-(grad)^*)
%Syntaxe: div(px,py)

[m,n]=size(px);

M=zeros(m,n);
Mx=M;
My=M;

Mx(2:m-1,:)=px(2:m-1,:)-px(1:m-2,:);
Mx(1,:)=px(1,:);
Mx(m,:)=-px(m-1,:);

My(:,2:n-1)=py(:,2:n-1)-py(:,1:n-2);
My(:,1)=py(:,1);
My(:,n)=-py(:,n-1);

M=Mx+My;


