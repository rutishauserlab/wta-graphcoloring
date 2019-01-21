%returns the home directory
%
%from matlab central
%
%
function userdir=getuserdir

if ispc; userdir= getenv('USERPROFILE');
else; userdir= getenv('HOME');
end