%code by: Martin.Strauch@lfb.rwth-aachen.de
%
%F = load_animal(<mat_file>)
%
function A = reconstruct_movie2(F, channel)

U = F.U{channel};
V = F.V{channel};
S = F.S{channel};

%reconstruct the data from SVD:
A = U*S*V';

