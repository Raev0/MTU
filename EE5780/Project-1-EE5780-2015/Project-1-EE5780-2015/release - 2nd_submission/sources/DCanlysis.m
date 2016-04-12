function [Vo,V_index,Io,I_index]=DCanlysis(R,C,L,V,I,MOS,MOS_MODEL,Y_AC,J_AC)
global  NR_Maxnum;
global Gminstep;

global Conv_step;
global Gmin_conv_bound;
global Gmin_ini;
global Gminlimit;
global Linear_Zone;


Gmin_conv_bound=1e-9;
Conv_step=20;

% Solves a set of nonlinear algebraic eqns
%% Modification
% 3.3: 28th Mar; ZCT; Add boundry condition Linear_Zone.
% 3.2: 28th Mar; ZCT; Bugfix: change converge judgement condition.
% 3.1: 25th Mar; ZCT; Bugfix: update Y and J in each NR start.
% 3.0: 25th Mar; ZCT; Bugfix: gmin and NR should be 2 loop
% 2.2: 24th Mar; ZCT; remove 1.9 modification: Tevenin stamping in TR mode
% 2.1: 24th Mar; ZCT; add bounding JUMP out of NR LOOP
% 2.0: 24th Mar; ZCT; add gmin method to acc the iteration convergency
% 1.9: 23th Mar; ZCT; L no more as 0v voltage: Tevenin stamping in TR mode
% 1.8: 23th Mar; ZCT; bugfix. when tr cal, can't not treat L as 0v voltage
% 1.7: 16th Mar; ZCT: NR convergence problem
% 1.7: 16th Mar; ZCT: reaction to the singular Y of test_ind.ckt (pending)
% 1.6: 16th Mar; ZCT: connect v0 to voltage 0
% 1.5: 16th Mar; ZCT: And Branch current calculation
% 1.4: 16th Mar; ZCT: change the DC inductor, treated as 0V DC V
% 1.3: 15th Mar; ZCT; For singular Y, let v=0 (doubting)
% 1.2: 12th Mar; ZCT; sketch NR LOOP for certain times NR_Maxnum
% 1.1: 10th Mar; ZCT; add the MOSFET stamping of 1st step NR
% 1.0: 8th Mar; ZCT; add the linear stamping



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

V_t_k=zeros(matrix_dim,1);
V_t_kk=zeros(matrix_dim,1);
Y=zeros(matrix_dim,matrix_dim)+Y_AC;
% voltage node0
% voltage node1~n
% current: Vs 1~n
% current: L as Vs
% ground: as Vs

%% Branch Current: each element(R;L;MOS;V) have a branch
% I_index made of strings, Io made of branch current
I_index=char(horzcat(R{2:L_R,1},L{2:L_L,1},V{2:L_V,1},I{2:L_I,1},MOS{2:L_MOS,1}));
I_index(I_index=='N',:)=[]; % delete 'NA'

Io=zeros(length(I_index),1);

%% stamp the the independent source J
J=zeros(matrix_dim,1)+J_AC;
% Independent current source row (1~node_num
if anyIs
    for count=2:L_I
        
        DC_I=FindDC(I); % get the DC value of the I and V
        J(find(node==I{count,3}),1)=J(find(node==I{count,3}),1)+DC_I(count-1);
    end
end
% Independent Voltage source row node_num~node_num


if anyVs
    for count=2:L_V
        
        DC_V=FindDC(V); % get the DC value of the I and V
        J(node_num+count-1,1)=J(node_num+count-1,1)+DC_V(count-1);
        
        Y(node_num+count-1,find(node==V{count,2}))=1; %from
        Y(find(node==V{count,2}),node_num+count-1)=1; %from
        Y(node_num+count-1,find(node==V{count,3}))=-1; %from
        Y(find(node==V{count,3}),node_num+count-1)=-1; %from
        
        
    end
end
% v0 connect to voltage of 0V

Y(matrix_dim,1)=1; %from
Y(1,matrix_dim)=1; %from

%% stamping: inductance is treated as shorted
% inductance is treated as DC voltage of 0V


if anyL
    for count=2:L_L
        
        Y(node_num+num_Vs+count-1,find(node==L{count,2}))=1; %from
        Y(find(node==L{count,2}),node_num+num_Vs+count-1)=1; %from
        Y(node_num+num_Vs+count-1,find(node==L{count,3}))=-1; %to
        Y(find(node==L{count,3}),node_num+num_Vs+count-1)=-1; %to
        
        
    end
    
end

%% stamp the the linear Y matrix
% conducter concern
% inducter also concern, because in the DC analyse, Induscter is seen as
% short
if anyR
    for count=2:L_R
        
        Y(find(node==R{count,2}),find(node==R{count,2}))=...
            Y(find(node==R{count,2}),find(node==R{count,2}))+1/R{count,4};
        Y(find(node==R{count,3}),find(node==R{count,3}))=...
            Y(find(node==R{count,3}),find(node==R{count,3}))+1/R{count,4};
        
        Y(find(node==R{count,2}),find(node==R{count,3}))=...
            Y(find(node==R{count,2}),find(node==R{count,3}))-1/R{count,4};
        Y(find(node==R{count,3}),find(node==R{count,2}))=...
            Y(find(node==R{count,3}),find(node==R{count,2}))-1/R{count,4};
    end
end

%% stamp the the nonlinear componient (MOSFET)
% MNA Stamping with Newton Raphason Method k=1
% beta(MU Cox W_diffusion L_diffusion) TH lamda
if anyMOS
    MOS_NR_state=zeros(num_MOS,NR_Maxnum,5);
    % MOS_NR_state (,,1)................Vgs_k
    % MOS_NR_state (,,2)................Vds_k
    % MOS_NR_state (,,3)................ids_k
    % MOS_NR_state (,,4)................gm_k
    % MOS_NR_state (,,5)................Gds_k
    Ids=zeros(num_MOS,1);
    Type_mark=zeros(num_MOS,1);
    
    

    
    % calculate ids_ini
    for count=2:L_MOS
        

        
        

       
        
        MODE=find(cell2mat(MOS_MODEL(2:L_MODEL,1))==MOS{count,8})+1;
        
        MU=MOS_MODEL{MODE,3};
        COS=MOS_MODEL{MODE,4};
        LAMBDA=MOS_MODEL{MODE,5};
        TH=MOS_MODEL{MODE,2};
        W_diffusion=MOS{count,6};
        L_diffusion=MOS{count,7};
        PN_type=MOS{count,5};
        
        Beta=MU*COS*W_diffusion/L_diffusion;
        
        if strcmp(PN_type,'N') || strcmp(PN_type,'n')
            vgs_ini=3;
            vds_ini=1;
            Type_mark(count-1)=1;
            if vds_ini<vgs_ini-TH && vgs_ini>TH % Linear region
                %MOS_NR_state(count-1,1,3)=Beta*((NMOS_vgs_ini-TH)*NMOS_vds_ini-NMOS_vds_ini^2/2);
                ids_ini=Beta*((vgs_ini-TH)*vds_ini-vds_ini^2/2)*(1+LAMBDA*vds_ini);
            elseif vds_ini>=vgs_ini-TH && vgs_ini>TH % saturation region
                ids_ini=Beta*(vgs_ini-TH)^2*(1+LAMBDA*vds_ini)/2;
            elseif vgs_ini<=TH
                ids_ini=0;
            end
        elseif strcmp(PN_type,'P') || strcmp(PN_type,'p')
            Type_mark(count-1)=-1;
            vgs_ini=-2;
            vds_ini=3;
            if vds_ini>vgs_ini-TH && vgs_ini<TH % Linear region
                %MOS_NR_state(count-1,1,3)=-Beta*((PMOS_vgs_ini-TH)*PMOS_vds_ini-PMOS_vds_ini^2/2);
                ids_ini=-Beta*((vgs_ini-TH)*vds_ini-vds_ini^2/2)*(1-LAMBDA*vds_ini);
            elseif vds_ini<=vgs_ini-TH && vgs_ini<TH % saturation region
                ids_ini=-Beta*(vgs_ini-TH)^2*(1-LAMBDA*vds_ini)/2;
            elseif vgs_ini>=TH
                ids_ini=0;
            end
            
        end
        
        MOS_NR_state(count-1,1,1)=vgs_ini;% ini vgs
        MOS_NR_state(count-1,1,2)=vds_ini;% ini vds
        MOS_NR_state(count-1,1,3)=ids_ini;
    end
    
      
    
    %%  MNA Stamping with Newton Raphason Method k=2
    
    n=2;
    B_conv=false;
    delta_V_t=zeros(matrix_dim,1);
    
   
    % attach the gmin to each drain node
    pos_matrix=zeros(matrix_dim,matrix_dim);
    
    pos_matrix=[ diag(ones(1,node_num)) zeros(node_num,matrix_dim-node_num);...
        zeros(matrix_dim-node_num,matrix_dim)];
    pos_matrix(1,2:node_num)=-1;
    pos_matrix(2:node_num,1)=-1;
    pos_matrix(1,1)=node_num;
    
    Y_Gmin=Y;
    J_Gmin=J;
    
    VI_Gmin_k=zeros(matrix_dim,1);
    
    Gmin=Gmin_ini; 
 
    while Gmin>Gminlimit
    
    
 
    
    
    Gmin_matrix=Gmin*pos_matrix;
    Y_Gmin=Y+Gmin_matrix;
    
    while  ~B_conv
        

        
        Y_NR=Y_Gmin;
        J_NR=J;
        
        if n>NR_Maxnum

            
            break;
        end
        
        % Loop 1: Calculate the results V array
        for count=2:L_MOS
            % Loop: every MOSFET
            MODE=find(cell2mat(MOS_MODEL(2:L_MODEL,1))==MOS{count,8})+1;
            
            MU=MOS_MODEL{MODE,3};
            COS=MOS_MODEL{MODE,4};
            LAMBDA=MOS_MODEL{MODE,5};
            TH=MOS_MODEL{MODE,2};
            W_diffusion=MOS{count,6};
            L_diffusion=MOS{count,7};
            PN_type=MOS{count,5};
            
            Beta=MU*COS*W_diffusion/L_diffusion;
            
            Vgs_k=MOS_NR_state(count-1,n-1,1);
            Vds_k=MOS_NR_state(count-1,n-1,2);
            ids_k=MOS_NR_state(count-1,n-1,3);
            
            
            if strcmp(PN_type,'N') || strcmp(PN_type,'n')
                
                if Vds_k<Vgs_k-TH && Vgs_k>TH % Linear region
                    %if Vds_k<Vgs_k-TH  % Linear region
                    %gm_k=Beta*Vds_k;
                    %Gds_k=-Beta*Vds_k;
                    gm_k=Beta*Vds_k*(LAMBDA*Vds_k+1);
                    Gds_k=Beta*(LAMBDA*Vds_k+1)*(Vgs_k-Vds_k-TH)+Beta*LAMBDA*(Vds_k*(Vgs_k-TH)-Vds_k^2/2);
                    ieq_k=ids_k-gm_k*Vgs_k-Gds_k*Vds_k;
                    %elseif Vds_k>=Vgs_k-TH % saturation region
                elseif Vds_k>=Vgs_k-TH && Vgs_k>TH % saturation region
                    gm_k=Beta*(1+LAMBDA*Vds_k)*(Vgs_k-TH);
                    Gds_k=Beta*LAMBDA*(Vgs_k-TH)^2/2;
                    ieq_k=ids_k-gm_k*Vgs_k-Gds_k*Vds_k;
                    
                elseif Vgs_k<=TH
                    
                    gm_k=0;
                    Gds_k=0;
                    ieq_k=0;
                end
            elseif strcmp(PN_type,'P') || strcmp(PN_type,'p')
                
                if Vds_k>Vgs_k-TH && Vgs_k<TH % Linear region
                    %gm_k=-Beta*Vds_k;
                    %Gds_k=Beta*Vds_k;
                    gm_k=-Beta*Vds_k*(1-LAMBDA*Vds_k);
                    Gds_k=-Beta*(1-LAMBDA*Vds_k)*(Vgs_k-Vds_k-TH)+Beta*LAMBDA*(Vds_k*(Vgs_k-TH)-Vds_k^2/2);
                    ieq_k=ids_k-gm_k*Vgs_k-Gds_k*Vds_k;
                elseif Vds_k<=Vgs_k-TH && Vgs_k<TH % saturation region
                    gm_k=-Beta*(1-LAMBDA*Vds_k)*(Vgs_k-TH);
                    Gds_k=Beta*LAMBDA*(Vgs_k-TH)^2/2;
                    ieq_k=ids_k-gm_k*Vgs_k-Gds_k*Vds_k;
                elseif Vgs_k>=TH
                    gm_k=0;
                    Gds_k=0;
                    ieq_k=0;
                end

            end
            MOS_NR_state(count-1,n-1,4)=gm_k;
            MOS_NR_state(count-1,n-1,5)=Gds_k;
            
            d_node=MOS{count,2};
            g_node=MOS{count,3};
            s_node=MOS{count,4};
            
            
            % stamping gm_k, Gds_k, ieq_k(count th MOSFET N/P TYPE)
            % stamping Gds_k as Conductor: drain to source
            Y_NR(find(node==d_node),find(node==d_node))=...
                Y_NR(find(node==d_node),find(node==d_node))+Gds_k;
            Y_NR(find(node==s_node),find(node==s_node))=...
                Y_NR(find(node==s_node),find(node==s_node))+Gds_k;
            Y_NR(find(node==d_node),find(node==s_node))=...
                Y_NR(find(node==d_node),find(node==s_node))-Gds_k;
            Y_NR(find(node==s_node),find(node==d_node))=...
                Y_NR(find(node==s_node),find(node==d_node))-Gds_k;
            % stamping ieq_k as independent current source: drain to source
            
            J_NR(find(node==d_node),1)=J_NR(find(node==d_node),1)-ieq_k;
            J_NR(find(node==s_node),1)=J_NR(find(node==s_node),1)+ieq_k;
            
            % stamping gm_k as VCCS: p=drain,k=gate; l=q=source
            Y_NR(find(node==d_node),find(node==g_node))=...
                Y_NR(find(node==d_node),find(node==g_node))+gm_k;%(p,k)
            Y_NR(find(node==s_node),find(node==s_node))=...
                Y_NR(find(node==s_node),find(node==s_node))+gm_k;%(q,l)
            
            Y_NR(find(node==d_node),find(node==s_node))=...
                Y_NR(find(node==d_node),find(node==s_node))-gm_k;%(p,l)
            Y_NR(find(node==s_node),find(node==g_node))=...
                Y_NR(find(node==s_node),find(node==g_node))-gm_k;%(q,k)
            
            % Calculate after stamping
        end
        
        V_t_kk=inv(Y_NR)*J_NR;
        
        
        V_t_kk(2:node_num)=V_t_kk(2:node_num)-(V_t_kk(2:node_num)<Linear_Zone).*V_t_kk(2:node_num);

        
        
        delta_V_t_kk=V_t_kk-V_t_k;
        
        
        
        
        
        
        V_conv_node=(abs(delta_V_t_kk(2:node_num))<Linear_Zone);
        B_conv=all(V_conv_node) &&...
            all(abs(delta_V_t_kk(2:node_num))<=abs(delta_V_t(2:node_num)));
        
        
        
        if B_conv
                
                break; 
           
        end
        
  %      out of bound: update
         delta_V_t=delta_V_t_kk;
         
        V_t_kk(2:node_num)=V_t_kk(2:node_num)-(V_t_kk(2:node_num)>3).*(V_t_kk(2:node_num)-3);
        V_t_kk(2:node_num)=V_t_kk(2:node_num)-(V_t_kk(2:node_num)<-3).*(V_t_kk(2:node_num)+3);
         V_t_k=V_t_kk; 
         
        
        
        for count=2:L_MOS
            % update the NR_state map by V_t_kk value
            
            Vgs_k=MOS_NR_state(count-1,n-1,1);
            Vds_k=MOS_NR_state(count-1,n-1,2);
            ids_k=MOS_NR_state(count-1,n-1,3);
            gm_k=MOS_NR_state(count-1,n-1,4);
            Gds_k=MOS_NR_state(count-1,n-1,5);
            
            
            Vgs_kk=V_t_kk(find(node==MOS{count,3}),1)-V_t_kk(find(node==MOS{count,4}),1);
            Vds_kk=V_t_kk(find(node==MOS{count,2}),1)-V_t_kk(find(node==MOS{count,4}),1);
            
            delta_Vgs=Vgs_kk-Vgs_k;
            delta_Vds=Vds_kk-Vds_k;
            
            %% update ids_kk
            MODE=find(cell2mat(MOS_MODEL(2:L_MODEL,1))==MOS{count,8})+1;
            
            MU=MOS_MODEL{MODE,3};
            COS=MOS_MODEL{MODE,4};
            LAMBDA=MOS_MODEL{MODE,5};
            TH=MOS_MODEL{MODE,2};
            W_diffusion=MOS{count,6};
            L_diffusion=MOS{count,7};
            PN_type=MOS{count,5};
            
            
            
            Beta=MU*COS*W_diffusion/L_diffusion;
            
            if strcmp(PN_type,'N') || strcmp(PN_type,'n')
                
                if Vds_kk<Vgs_kk-TH && Vgs_kk>TH % Linear region
                    ids_kk=Beta*((Vgs_kk-TH)*Vds_kk-Vds_kk^2/2)*(1+LAMBDA*Vds_kk);
                elseif Vds_kk>=Vgs_kk-TH && Vgs_kk>TH % saturation region
                    ids_kk=Beta*(Vgs_kk-TH)^2*(1+LAMBDA*Vds_kk)/2;
                elseif Vgs_kk<=TH
                    ids_kk=0;
                end
            elseif strcmp(PN_type,'P') || strcmp(PN_type,'p')
                
                if Vds_kk>Vgs_kk-TH && Vgs_kk<TH % Linear region
                    ids_kk=-Beta*((Vgs_kk-TH)*Vds_kk-Vds_kk^2/2)*(1-LAMBDA*Vds_kk);
                elseif Vds_kk<=Vgs_kk-TH && Vgs_kk<TH % saturation region
                    ids_kk=-Beta*(Vgs_kk-TH)^2*(1-LAMBDA*Vds_kk)/2;
                elseif Vgs_kk>=TH
                    ids_kk=0;
                end
                
            end
            
            

            
            
            % update for next NR at new linear zone
            MOS_NR_state(count-1,n,1)=Vgs_kk;
            MOS_NR_state(count-1,n,2)=Vds_kk;
            MOS_NR_state(count-1,n,3)=ids_kk;
            
            
        
            
           
            
            
        end
        
        
        n=n+1;
    end
      Gmin=Gmin/10;
  end   
    
    
    
    
    Ids=MOS_NR_state(:,n-1,3).*Type_mark;
else
    if det(Y)==0
    else
        V_t_kk=inv(Y)*J;
    end
    
end

% add the node which was remove due to conductor
% V_t is shrink ,index by node
% V_r is lefted one, index by node_remove
% V_o is the whole


Vo=V_t_kk;
V_index=node;
% update Io:
% (1~num_R) is Resister current;
% (num_R+1~num_R+num_L) is inductor current;
% (num_R+num_L+1~num_R+num_L+num_Vs) is Voltage current;
% (num_R+num_L+1~num_R+num_L+num_Vs+num_Is) is Current Source current;
% (num_R+num_L+num_MOS+1~num_R+num_L+num_MOS+num_Vs+num_Is) is MOS ds current;


for n=1:length(I_index)
    if anyR && n<=num_R
        R_value=R{n+1,4};
        FROM_node_v=V_t_kk(find(node==R{n+1,2}));
        To_node_v=V_t_kk(find(node==R{n+1,3}));
        Io(n)=(FROM_node_v-To_node_v)/R_value;
    elseif anyL && n<=num_R+num_L && n>num_R
        % inductance was treated as 0 V voltage
        index=node_num+n-num_R+num_Vs; % (n-num_R)+num_Vs+node_num
        Io(n)=V_t_kk(index);
        
        
    elseif anyVs && n<=num_R+num_L+num_Vs && n>num_R+num_L
        index=node_num+n-num_R-num_L;       % (n-num_R-num_L)+node_num
        Io(n)=V_t_kk(index);
    elseif anyIs && n<=num_R+num_L+num_Vs+num_Is && n>num_R+num_L+num_Vs
        index=(n-num_R-num_L-num_Vs)+1;   % (n-num_R-num_L-num_Vs)+1
        Io(n)=I{index,5};
    elseif anyMOS && n>num_R+num_L+num_Vs+num_Is
        index= n-num_R-num_L-num_Vs-num_Is; %n-num_R-num_L-num_Vs-num_Is
        
        Io(n)=Ids(index); % PMOS should be isd??????
        
        
    end
end


end
