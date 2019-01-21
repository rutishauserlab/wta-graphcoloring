%
% returns the hermitian part of a matrix
%
% ' is the complex conjugate in matlab
%
%
function B = getSymPart( A )

B = 1/2 .* ( A + A' );