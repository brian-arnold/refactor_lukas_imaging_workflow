%code by: Martin.Strauch@lfb.rwth-aachen.de
%
%smooth convex cone analysis (scca): convex cone analysis with a smoothness constraint
%
function [C, X, C_indices, reconstruction_accuracy, norm_error] = scca3d(A, width, height, depth, k, threshold)

original_A = A;

m = size(A,1);
n = size(A,2);

C = zeros(m, k);
X = zeros(k, n);

C_indices  = zeros(1,k);

%initialise with max norm vector
maxi = -999999999999;
max_pos = -1;
for (i=1:n)
   value = norm(A(:,i));
   if(value>maxi)
   maxi=value;
   max_pos=i;
   end
end
p = max_pos;


for (i = 1:k)
  C_indices(i) = p;  
  t = A(:,p);
  t=t/norm(t);
  s = A'*t;
   
  %non-negativity constraint
  s(find(s<0))=0;
  
  %non-negativity & spatial smoothness constraint:
  %s = local_smoothness_filter_3d(s, width, height, depth, p, threshold); 
  
  C(:,i) = t;
  X(i,:) = s;
   
  A = A - (t*s');
  
  %max column norm in A:
  colnorms=zeros(1,n);
  for (j = 1:n)
       colnorms(j) = norm(A(:,j));
  end
  
  [Y, I] = max(colnorms);
  p = I;

end

disp(C_indices);

%for NNCX, uncomment:
%C = original_A(:,C_indices);
%X = pinv(C) * original_A; 

% convex quadratic programming to find optimal non-negative X:
%  
%   cvx_begin quiet
%     
%        variable X(c,n);
%        X>=0;
%        minimize( norm(original_A-(C*X) ,'fro')  );
%        
%   cvx_end


fr_A                    = norm(original_A, 'fro');
reconstruction_error    = norm(original_A-(C*X), 'fro');

reconstruction_accuracy = 100 - ( (reconstruction_error * reconstruction_error) / (fr_A * fr_A) * 100);

norm_error = (reconstruction_error * reconstruction_error);




