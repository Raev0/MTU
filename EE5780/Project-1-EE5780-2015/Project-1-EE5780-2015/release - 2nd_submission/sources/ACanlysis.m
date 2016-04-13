function[Vo,V_index,Io,I_index,F]=ACanlysis(R,C,L,V,I,MOS,MOS_MODEL,AC)
% Solves a set of complex linear eqns
% Modification
% 1.1: 18th Mar; ZCT; using 2nd order model AWE to cal the AC responce(AWE=true)
% 1.0: 18th Mar; ZCT; using normal AC model to calculate the AC responce(AWE=false)
% just for 'rcsimple.ckt', which DC operating point is zero.
global AWE;

r2d=180/pi;

[L_I,W_I]=size(I);
[L_V,W_V]=size(V);
[L_L,W_L]=size(L);
[L_C,W_C]=size(C);
[L_R,W_R]=size(R);
[L_MODEL,W_MODEL]=size(MOS_MODEL);
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
% for 'rcsimple.ckt', DC bias is 0
node=unique([reshape(cell2mat(MOS(2:L_MOS,2:4)),1,[]),...
    reshape(cell2mat(V(2:L_V,2:3)),1,[]),...
    reshape(cell2mat(R(2:L_R,2:3)),1,[]),...
    reshape(cell2mat(C(2:L_C,2:3)),1,[]),...
    reshape(cell2mat(I(2:L_I,2:3)),1,[]),...
    reshape(cell2mat(L(2:L_L,2:3)),1,[]),...
    ]);

ground_node=1;
matrix_dim=length(node)-ground_node+num_L+num_Vs;

f_start=AC{3};
f_stop=AC{4};
f_num=AC{2};


f_step=(f_stop-f_start)/(f_num-1);
freq=log10(2*pi*f_start:2*pi*f_step:2*pi*f_stop);


%% Solve for dc operating point
% for 'rcsimple.ckt',  DC operating point is all zero,the input is only a step
%% small signal analysis
%small signal ac modal: M_G,B

%


M_C=spalloc(matrix_dim, matrix_dim, 6*matrix_dim);
M_G=spalloc(matrix_dim, matrix_dim, 6*matrix_dim);
M_B=spalloc(matrix_dim, 1, 6*matrix_dim);
Vt_s=spalloc(matrix_dim, f_num, 6*matrix_dim);

%%complex impedance equation stamping
% small signal voltage stamping

if anyVs && num_Vs==1
    Ve=str2num(cell2mat(V{2,4}));
    V_FROM=find(node==V{2,2})-ground_node;
    V_TO=find(node==V{2,3})-ground_node;
    if V_FROM~=0
        M_G(matrix_dim,V_FROM)=-1;
        M_G(V_FROM,matrix_dim)=1;
        M_B(matrix_dim)=-1;
    end
    if V_TO~=0
        M_G(matrix_dim,V_TO)=1;
        M_B(V_TO)=-1;
    end
end


% conductor stamping
if anyC
    for n=2:L_C
        C_value=C{n,4};
        C_FROM=find(node==C{n,2})-ground_node;
        C_TO=find(node==C{n,3})-ground_node;
        if C_FROM~=0 && C_TO~=0
            M_C(C_FROM,C_FROM)=M_C(C_FROM,C_FROM)+C_value;
            M_C(C_TO,C_TO)=M_C(C_TO,C_TO)+C_value;
            M_C(C_FROM,C_TO)=M_C(C_FROM,C_TO)-C_value;
            M_C(C_TO,C_FROM)=M_C(C_TO,C_FROM)-C_value;
        elseif C_FROM==0
            M_C(C_TO,C_TO)=M_C(C_TO,C_TO)+C_value;
        elseif C_TO==0
            M_C(C_FROM,C_FROM)=M_C(C_FROM,C_FROM)+C_value;
        end
        
    end
end
% inductor stamping
if anyL
    offset=length(node)-ground_node;
    for n=2:L_L
        
        L_value=L{n,4};
        L_FROM=find(node==L{n,2})-ground_node;
        L_TO=find(node==L{n,3})-ground_node;
        
        M_C(n+offset,n+offset)=L_value;
        
        if L_FROM~=0 && L_TO~=0
            M_G(L_FROM,n+offset)=1;
            M_G(L_TO,n+offset)=-1;
            M_G(n+offset,L_FROM)=-1;
            M_G(n+offset,L_TO)=1;
        elseif L_FROM==0
            M_G(L_TO,n+offset)=-1;
            M_G(n+offset,L_TO)=1;
        elseif L_To==0
            M_G(L_FROM,n+offset)=1;
            M_G(n+offset,L_FROM)=-1;
        end
        
    end
end
% Resistance stamping
if anyR
    
    for n=2:L_R
        
        R_value=R{n,4};
        R_FROM=find(node==R{n,2})-ground_node;
        R_TO=find(node==R{n,3})-ground_node;
        
        if R_FROM~=0 && R_TO~=0
            M_G(R_FROM,R_FROM)=M_G(R_FROM,R_FROM)+1/R_value;
            M_G(R_TO,R_TO)=M_G(R_TO,R_TO)+1/R_value;
            M_G(R_FROM,R_TO)=M_G(R_FROM,R_TO)-1/R_value;
            M_G(R_TO,R_FROM)=M_G(R_TO,R_FROM)-1/R_value;
        elseif L_FROM==0
            M_G(R_TO,R_TO)=M_G(R_TO,R_TO)+1/R_value;
        elseif L_To==0
            M_G(R_FROM,R_FROM)=M_G(R_FROM,R_FROM)+1/R_value;
        end
        
    end
end
%% solving the AC reponse
if AWE % AWE Metheod, reduce the model to a 2nd order model.
    warning('off');
    M_A=-inv(M_G)*M_C; %determine matrix A
    M_R=inv(M_G)*M_B; %determine matrix R
    M_m=zeros(2,2);
    for n=1:matrix_dim
        
        
        indicate_str=['\t(' num2str(n) '/' num2str(matrix_dim) ')Cal the AWE for all ' num2str(matrix_dim) ' nodes....\n'];
        fprintf(indicate_str);
        M_L=zeros(1,matrix_dim);
        M_L(n)=1;
        

        for count=1:f_num
            
            %M_moment=zeros(4,matrix_dim); % moments of each node
           
            m_0=M_L*M_R;% LR
            m_1=M_L*M_A*M_R;% LAR
            m_2=2*M_L*M_A*M_A*M_R; % 2LAAR
            m_3=6*M_L*M_A*M_A*M_A*M_R;% 6LAAAR
            
            M_m(1,1)=m_0;
            M_m(1,2)=m_1;
            M_m(2,1)=m_1;
            M_m(2,2)=m_2;
            
            bcef=-inv(M_m)*[m_2;m_3]; % bcef(2)=b1,bcef(1)=b2
            acef=[m_0;m_0*bcef(2)+m_1];% acef(1)=a0,acef(2)=a1
            
            % freq sweep
            f=f_start+(count-1)*f_step;
            w=2*pi*f;
            s=i*w;
            
            zero_line=acef(1)+acef(2)*s;
            pole_line=1+bcef(2)*s+bcef(1)*s^2;
            
            Vt_s(n,count)=zero_line/pole_line;
        end
    end
    
    warning('on');
else % noarmal normal AC model: solving Result matrix directly
    
    
    % s*M_C*Res=-M_G*Res+M_B*Ve
    for n=1:f_num
        f=f_start+(n-1)*f_step;
        indicate_str=['\t(' num2str(n) '/' num2str(f_num) ')Cal the AC response for ' num2str(f)  'Hz....\n'];
        fprintf(indicate_str);
        
        w=2*pi*f;
        s=i*w;
        M_temp=s.*M_C+M_G;
        V_temp=inv(M_temp);
        Vt_s(:,n)=V_temp*M_B*complex(Ve);
        %Vt_s(:,n)=inv(M_temp)*M_B*Ve;
        
    end
end



%subplot(2,1,1)
%plot(freq,20*log10(abs(Vt_s(2,:))));
%subplot(2,1,2)
%plot(freq,angle(Vt_s(2,:))*r2d);

%%bode: magnitude and phase

Vo=Vt_s;
V_index=node;
Io=[6,6,6,6,6];
I_index=[6,6,6,6,6];
F=freq;

end