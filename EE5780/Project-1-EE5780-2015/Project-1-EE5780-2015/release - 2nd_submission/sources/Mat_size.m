function [matrix_dim,node_num]=Mat_size(R,C,L,V,I,MOS)
[L_I,W_I]=size(I);
[L_V,W_V]=size(V);
[L_L,W_L]=size(L);
[L_C,W_C]=size(C);
[L_R,W_R]=size(R);
[L_MOS,W_MOS]=size(MOS);

anyVs= V{2,2} ~=0 || V{2,3} ~=0;
anyIs= I{2,2} ~=0 || I{2,3} ~=0;
anyL= L{2,2} ~=0 || L{2,3} ~=0;
anyC= C{2,2} ~=0 || C{2,3} ~=0;
anyR= R{2,2} ~=0 || R{2,3} ~=0;
anyMOS= MOS{2,2} ~=0 || MOS{2,3} ~=0|| MOS{2,4} ~=0;

num_Vs= anyVs*(L_V-1);
num_Is= anyIs*(L_I-1);
num_L= anyL*(L_L-1);
num_C= anyC*(L_C-1);
num_R= anyR*(L_R-1);
num_MOS= anyMOS*(L_MOS-1);

% Creat Inietial Matrix
% note: for the case grounded node, the ckt is index-0, but the array is
% start from 1.


node=unique([reshape(cell2mat(MOS(2:L_MOS,2:4)),1,[]),...
    reshape(cell2mat(V(2:L_V,2:3)),1,[]),...
    reshape(cell2mat(R(2:L_R,2:3)),1,[]),...
    reshape(cell2mat(C(2:L_C,2:3)),1,[]),...
    reshape(cell2mat(I(2:L_I,2:3)),1,[]),...
    reshape(cell2mat(L(2:L_L,2:3)),1,[]),...
    ]);


node_num=length(node);

ground_v_num=1;

matrix_dim=node_num+num_Vs+num_L+ground_v_num;

end