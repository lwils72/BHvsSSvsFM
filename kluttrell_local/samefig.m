function samefig(figchild,figparent)

% samefig(figchild,figparent)

for i=1:numel(figchild)
figure(figchild(i))
end
set(figchild,'position',get(figparent,'position'))
