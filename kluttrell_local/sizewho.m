%function sizewho

ans=whos;
%GB=sum([ans.bytes])/1e9;
[num2str(sum([ans.bytes])/1e9),' GB']
