%
%returns the skew symmetric part of a matrix
%
%urut/CC2014
function B = getSkewPart( A )

B = 1/2 .* ( A - A' );