%rescale range
%http://math.stackexchange.com/questions/43698/range-scaling-problem
%
%urut/aug15
function y=rescale_range(x,A,B,C,D)

%range x in [A-B]  
%to y in [C-D]


y=C+(D-C)*((x-A)/(B-A));