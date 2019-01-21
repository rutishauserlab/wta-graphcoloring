%
% converts active units (DD) to set number, as a function of time
%
%matches only excitatory units
%
function DD_matchedSet = convert_DD_toSetNr(DD,allSets_sorted, indsExcit)

%%

toMatch = DD(indsExcit,:)';

[out1,out2]=ismember( toMatch, allSets_sorted ,'rows' );

DD_matchedSet=out2;

% %%==== below, same but not optimized
% %% --- next: plot set # as funct of time to see if it follows the gradient
% DD_matchedSet=zeros(1,size(DD,2));
% parfor i=1:size(DD,2)
%     
%     activeSet = (DD(:,i));
%     
%     % find which subset this is
%     matchedSet=0;
%     
%     
%     
%     [out] = ismember( allSets_sorted,activeSet(indsExcit)', 'rows');
%     
%     %[out1,out2] = ismember([1 2 3;2 3 4;5 6 7; 9 10 11],[ 9 10 11;1 2 3 ], 'rows')
% 
%     
%     matchedSet = find(out);
%    
% %     
% %     for k=1:size(allSets_sorted,1)
% %     
% %         out = ismember(v,v(2,:),'rows');
% %         
% %         
% %         if isequal(allSets_sorted(k,:), activeSet(indsExcit)' )
% %             matchedSet=k;
% %             break;
% %         end
% %     end
%     
%     
%     DD_matchedSet(i) = matchedSet;
% % end