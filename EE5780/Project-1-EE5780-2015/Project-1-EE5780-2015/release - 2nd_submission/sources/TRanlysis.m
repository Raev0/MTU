function [Vo,V_indexo,Io,I_indexo,T]=TRanlysis(R,C,L,V,I,MOS,MOS_MODEL,TR)
%Solves a set of nonlinear different eqns
% Modification
% 2.0: 25th Mar; ZCT; MOSFET C stamping, new norton equivelent(pending)
% 1.2: 24th Mar; ZCT; L no more as 0v voltage: L tevenin method
% 1.1: 16th Mar; ZCT; deal with the PWL voltage
% 1.0: 16th Mar; ZCT; complete the logic flow
global waiting_time;
t_end=TR{1,3};
t_delta=TR{1,2};
step=int32(t_end/t_delta);

t=0:t_delta:t_end;


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

[matrix_dim,node_num]=Mat_size(R,C,L,V,I,MOS);


%Norton equivelent C and L, thus dimension remains.


%establish a table of num_Vs x step; indicate the Vs instant value at
%each step
fprintf('\tScanning the input voltages....\n');
Vs_step=cell(num_Vs,step+1);
for n=2:L_V
    if strcmp(V{n,4},'DC') || strcmp(V{n,4},'dc')
        Vs_step(n-1,:)=V(n,5);
    elseif strcmp(V{n,4},'PWL') || strcmp(V{n,4},'pwl')
        Vs_step{n-1,1}=V{n,5};
        m=5;
        for count=2:step+1
            if double(count)==501
                xisix=1;
            end
            
            if isempty(V{n,m+1})
                Vs_k=V{n,m};
            elseif double(count)*t_delta==V{n,m+1}
                Vs_k=V{n,m+2};
                
            else
                if double(count)*t_delta>V{n,m+1}
                    m=m+2;
                end
                
                if double(count)*t_delta<V{n,m+1}
                    if m==5
                        last_t=0;
                    else
                        last_t=V{n,m-1};
                    end
                    last_v=V{n,m};
                    
                    next_v=V{n,m+2};
                    
                    
                    next_t=V{n,m+1};
                    
                    Vs_k=last_v+(double(count)*t_delta-last_t)*(next_v-last_v)/(next_t-last_t);
                    
                    
                    
                    
                end
                
                
            end
            Vs_step{n-1,count}=Vs_k;
            
        end
    end
end
pause(waiting_time);
fprintf('\t\tfinished!\n');

%% Circuit Calulation
%%IF DC analysis
%%  initial DC soltion for V_0 I_0 t_kk=t_k+delta_t
V_temp=V(1:L_V,1:5);
V_temp(2:L_V,5)=Vs_step(:,1);
Y_AC0=zeros(matrix_dim,matrix_dim);
J_AC0=zeros(matrix_dim,1);
[DC_V,V_index,DC_I,I_index]=DCanlysis(R,C,L,V_temp,I,MOS,MOS_MODEL,Y_AC0,J_AC0);




% energy elements states
num_EE=num_C+2*num_MOS+num_L; % energy elements number: C,C_MOS,L
EE_TR_state=zeros(num_EE,step+1,2);
% EE_TR_state (,,1)................i
% EE_TR_state (,,2)................v


%%  rocord v
AC_Voltage=zeros(length(V_index),step+1);
AC_Current=zeros(length(I_index),step+1);

AC_Voltage(:,1)=DC_V(1:length(V_index),1);
AC_Current(:,1)=DC_I;

V_STEP=DC_V;%ini
I_STEP=DC_I;%ini

fprintf('\tTR methods is pending.... please wait\n');
pause(waiting_time+1);
for count=1:step
    Y_AC=zeros(matrix_dim,matrix_dim);
    J_AC=zeros(matrix_dim,1);
    
    if count==1
        % initial the EE_TR_state(num_EE,1,:)
        EE_TR_state(1:(num_C+2*num_MOS),1,1)=0; % MOSFET C and C Current
        EE_TR_state((num_C+2*num_MOS):num_EE,1,2)=0; % L voltage
        % update energy elements states
        if anyC
            for n=1:num_C
                Vc0=DC_V(find(V_index==C{n+1,2}))-DC_V(find(V_index==C{n+1,3}));
                EE_TR_state(n,1,2)=Vc0;
            end
        end
        if anyL
            for n=1:num_L
                
                iL0=DC_I(Loc(I_index,cell2mat(L{n+1,1})));
                EE_TR_state(num_C+2*num_MOS+n,1,1)=iL0;
            end
        end
        if anyMOS
            %% MOSFET C equavelent
            % initial the C value of mosfet: Cgs_k, Cgd_k
            C_MOS=cell(2*num_MOS+1,5);%(mos,gs/ds,from,to,value)
            C_MOS(1,:)={'MOS','gs/gd','from','to','value'};
            for n=1:num_MOS
                
                MODE=find(cell2mat(MOS_MODEL(2:L_MODEL,1))==MOS{n+1,8})+1;
                
                Cox=MOS_MODEL{MODE,4};
                TH=MOS_MODEL{MODE,2};
                W_diffusion=MOS{n+1,6};
                L_diffusion=MOS{n+1,7};
                PN_type=MOS{n+1,5};
                %vgs, vds , TH
                %determin the MOS state
                %determine Cgs and Cgd (base not connect to componients)
                
                g_node=MOS{n+1,3};
                s_node=MOS{n+1,4};
                d_node=MOS{n+1,2};
                vds=DC_V(find(V_index==g_node))-DC_V(find(V_index==s_node));
                vgs=DC_V(find(V_index==d_node))-DC_V(find(V_index==s_node));
                
                if strcmp(PN_type,'P') || strcmp(PN_type,'p')
                    linear_state=(vds>vgs-TH && vgs<TH);
                    coff_state=(vgs>=TH);
                    satu_state=(vds<=vgs-TH && vgs<TH);
                elseif strcmp(PN_type,'N') || strcmp(PN_type,'n')
                    linear_state=(vds<vgs-TH && vgs>TH);
                    coff_state=(vgs<=TH);
                    satu_state=(vds>=vgs-TH && vgs>TH);
                end
                
                if linear_state
                    Cgs_kk=Cox*W_diffusion*L_diffusion/2;
                    Cgd_kk=Cox*W_diffusion*L_diffusion/2;
                elseif satu_state
                    Cgs_kk=Cox*W_diffusion*L_diffusion*2/3;
                    Cgd_kk=Cox*W_diffusion*L_diffusion/3;
                elseif coff_state
                    Cgs_kk=0;
                    Cgd_kk=0;
                    
                end
                %
                
                

                
                
                C_MOS(2*n,:)={n,1,g_node,s_node,Cgs_kk};
                C_MOS(2*n+1,:)={n,2,g_node,d_node,Cgd_kk};
                
                EE_TR_state(num_C+2*n-1,1,2)=vgs;
                EE_TR_state(num_C+2*n,1,2)=vds;
                
            end
        end
    end
    %%  for delta T, calculate the campanion Model for C and L (TR)
    % linear C
    % 1. stamping Geq 2. stamping Ieq
    if anyC
        for n=2:L_C
            
            Vc_k=EE_TR_state(n-1,count,2);
            ic_k=EE_TR_state(n-1,count,1);
            Ieq=2*C{n,4}*Vc_k/t_delta+ic_k;
            Geq=2*C{n,4}/t_delta;
            %1. stamping Geq
            Y_AC(find(V_index==C{n,2}),find(V_index==C{n,2}))=...
                Y_AC(find(V_index==C{n,2}),find(V_index==C{n,2}))+Geq;
            Y_AC(find(V_index==C{n,3}),find(V_index==C{n,3}))=...
                Y_AC(find(V_index==C{n,3}),find(V_index==C{n,3}))+Geq;
            
            Y_AC(find(V_index==C{n,2}),find(V_index==C{n,3}))=...
                Y_AC(find(V_index==C{n,2}),find(V_index==C{n,3}))-Geq;
            Y_AC(find(V_index==C{n,3}),find(V_index==C{n,2}))=...
                Y_AC(find(V_index==C{n,3}),find(V_index==C{n,2}))-Geq;
            
            %2. stamping Ieq
            J_AC(find(V_index==C{n,2}),1)=J_AC(find(V_index==C{n,2}),1)+Ieq;
            J_AC(find(V_index==C{n,3}),1)=J_AC(find(V_index==C{n,3}),1)-Ieq;
        end
    end
    % nonlinear
    if anyMOS
        for n=1:num_MOS
            
            Vgs_k=EE_TR_state(2*n-1+num_C,count,2);
            Vgd_k=EE_TR_state(2*n+num_C,count,2);
            
            igs_k=EE_TR_state(2*n-1+num_C,count,1);
            igd_k=EE_TR_state(2*n+num_C,count,1);
            
            Ieq_gs=2*C_MOS{2*n,5}*Vgs_k/t_delta+igs_k;
            Geq_gs=2*C_MOS{2*n,5}/t_delta;
            
            
            Ieq_gd=2*C_MOS{2*n+1,5}*Vgd_k/t_delta+igd_k;
            Geq_gd=2*C_MOS{2*n+1,5}/t_delta;
            
            
                g_node=MOS{n+1,3};
                s_node=MOS{n+1,4};
                d_node=MOS{n+1,2};
                %1. stamping Geq of g s
                Y_AC(find(V_index==g_node),find(V_index==g_node))=...
                    Y_AC(find(V_index==g_node),find(V_index==g_node))+Geq_gs;
                Y_AC(find(V_index==s_node),find(V_index==s_node))=...
                    Y_AC(find(V_index==s_node),find(V_index==s_node))+Geq_gs;
                
                Y_AC(find(V_index==g_node),find(V_index==s_node))=...
                    Y_AC(find(V_index==g_node),find(V_index==s_node))-Geq_gs;
                Y_AC(find(V_index==s_node),find(V_index==g_node))=...
                    Y_AC(find(V_index==s_node),find(V_index==g_node))-Geq_gs;
                %2. stamping Geq of g d
                Y_AC(find(V_index==g_node),find(V_index==g_node))=...
                    Y_AC(find(V_index==g_node),find(V_index==g_node))+Geq_gd;
                Y_AC(find(V_index==d_node),find(V_index==d_node))=...
                    Y_AC(find(V_index==d_node),find(V_index==d_node))+Geq_gd;
                
                Y_AC(find(V_index==g_node),find(V_index==d_node))=...
                    Y_AC(find(V_index==g_node),find(V_index==d_node))-Geq_gd;
                Y_AC(find(V_index==d_node),find(V_index==g_node))=...
                    Y_AC(find(V_index==d_node),find(V_index==g_node))-Geq_gd;
                %3. stamping Ieq of g s 
                J_AC(find(V_index==g_node),1)=J_AC(find(V_index==g_node),1)+Ieq_gs;
                J_AC(find(V_index==s_node),1)=J_AC(find(V_index==s_node),1)-Ieq_gs;
                %4. stamping Ieq of g d 
                J_AC(find(V_index==g_node),1)=J_AC(find(V_index==g_node),1)+Ieq_gd;
                J_AC(find(V_index==d_node),1)=J_AC(find(V_index==d_node),1)-Ieq_gd;
        end
    end
    % inductance
    if anyL
        for n=2:L_L
            Req_k=2*L{n,4}/t_delta;
            VL_k=EE_TR_state(n-1+num_C+2*num_MOS,count,2);
            iL_k=EE_TR_state(n-1+num_C+2*num_MOS,count,1);
            Veq_k=2*L{n,4}*iL_k/t_delta+VL_k;
            %1. stamping Req_k
            Y_AC(node_num+num_Vs+n-1,node_num+num_Vs+n-1)=-Req_k;
            
            %2. stamping Veq
            
            J_AC(node_num+num_Vs+n-1,1)=-Veq_k;
        end
    end
    %% update the voltage
    if anyVs
        V_temp(2:L_V,5)=Vs_step(:,count+1);
    end
    %%  DC_solution for t_kk
    indicate_str=['\t(' num2str(count) '/' num2str(step) ')Calculate the DC response for ' num2str(double(count)*t_delta)  's....\n'];
    fprintf(indicate_str);
    [V_STEP,V_index,I_STEP,I_index]=DCanlysis(R,C,L,V_temp,I,MOS,MOS_MODEL,Y_AC,J_AC);
    % update energy elements states
    if anyC
        for n=1:num_C
            Vc_k=EE_TR_state(n,count,2);
            ic_k=EE_TR_state(n,count,1);
            Ieq_k=2*C{n+1,4}*Vc_k/t_delta+ic_k;
            Geq_k=2*C{n+1,4}/t_delta;
            
            % update
            Vc_kk=V_STEP(find(V_index==C{n+1,2}))-V_STEP(find(V_index==C{n+1,3}));
            ic_kk=Vc_kk*Geq_k-Ieq_k;
            
            EE_TR_state(n,count+1,2)=Vc_kk;
            EE_TR_state(n,count+1,1)=ic_kk;
        end
        
    end
    if anyL
        for n=1:num_L
            
            VL_k=EE_TR_state(n+num_C+2*num_MOS,count,2);
            iL_k=EE_TR_state(n+num_C+2*num_MOS,count,1);
            Req_k=2*L{n+1,4}/t_delta;
            Veq_k=2*L{n+1,4}*iL_k/t_delta+VL_k;
            
            iL_kk=I_STEP(Loc(I_index,cell2mat(L{n+1,1})));
            vL_kk=iL_kk*Req_k-Veq_k;
            
            % vl update
            EE_TR_state(num_C+2*num_MOS+n,count+1,1)=iL_kk;
            EE_TR_state(num_C+2*num_MOS+n,count+1,2)=vL_kk;
        end
    end
    if anyMOS
        for n=1:num_MOS
            %% update the C value: Cgs_k, Cgd_k , Cgs_kk , Cgd_kk
            
            MODE=find(cell2mat(MOS_MODEL(2:L_MODEL,1))==MOS{n+1,8})+1;
            
            Cox=MOS_MODEL{MODE,4};
            TH=MOS_MODEL{MODE,2};
            W_diffusion=MOS{n+1,6};
            L_diffusion=MOS{n+1,7};
            PN_type=MOS{n+1,5};
            
            g_node=MOS{n+1,3};
            s_node=MOS{n+1,4};
            d_node=MOS{n+1,2};
            vgs=V_STEP(find(V_index==g_node))-V_STEP(find(V_index==s_node));
            vds=V_STEP(find(V_index==d_node))-V_STEP(find(V_index==s_node));
      
            %vgs, vds , TH
            %determin the MOS state
            %determine Cgs and Cgd (base not connect to componients)
            if strcmp(PN_type,'P') || strcmp(PN_type,'p')
                linear_state=(vds>vgs-TH && vgs<TH);
                coff_state=(vgs>=TH);
                satu_state=(vds<=vgs-TH && vgs<TH);
            elseif strcmp(PN_type,'N') || strcmp(PN_type,'n')
                linear_state=(vds<vgs-TH && vgs>TH);
                coff_state=(vgs<=TH);
                satu_state=(vds>=vgs-TH && vgs>TH);
            end
            
            if  linear_state
                Cgs_kk=Cox*W_diffusion*L_diffusion/2;
                Cgd_kk=Cox*W_diffusion*L_diffusion/2;
            elseif satu_state
                Cgs_kk=Cox*W_diffusion*L_diffusion*2/3;
                Cgd_kk=Cox*W_diffusion*L_diffusion/3;
            elseif coff_state
                Cgs_kk=0;
                Cgd_kk=0;
            end
            %
            C_MOS(2*n,:)={n,1,g_node,s_node,Cgs_kk};
            C_MOS(2*n+1,:)={n,2,g_node,d_node,Cgd_kk};
            
            %% calculate Vc_mos Ic_mos
            
            Vgs_k=EE_TR_state(2*n-1+num_C,count,2);
            Vgd_k=EE_TR_state(2*n+num_C,count,2);
            
            igs_k=EE_TR_state(2*n-1+num_C,count,1);
            igd_k=EE_TR_state(2*n+num_C,count,1);
            
            Ieq_gs_k=2*C_MOS{2*n,5}*Vgs_k/t_delta+igs_k;
            Geq_gs_k=2*C_MOS{2*n,5}/t_delta;
            
            Ieq_gd_k=2*C_MOS{2*n+1,5}*Vgd_k/t_delta+igd_k;
            Geq_gd_k=2*C_MOS{2*n+1,5}/t_delta;
            
            Vgs_kk=V_STEP(find(V_index==g_node))-V_STEP(find(V_index==s_node));
            Vgd_kk=V_STEP(find(V_index==g_node))-V_STEP(find(V_index==d_node));
            
            
            igs_kk=Vgs_kk*Geq_gs_k-Ieq_gs_k;
            igd_kk=Vgd_kk*Geq_gd_k-Ieq_gd_k;
            % ic update
            
            if coff_state
                igs_kk=0;
                igd_kk=0;
            end
            
            
            
            EE_TR_state(2*n-1+num_C,count+1,2)=Vgs_kk;
            EE_TR_state(2*n+num_C,count+1,2)=Vgd_kk;
            EE_TR_state(2*n-1+num_C,count+1,1)=igs_kk;
            EE_TR_state(2*n+num_C,count+1,1)=igd_kk;
            
            
            
        end
        
    end
    
    %%  if t_kk=t_final end else t_k=t_kk, repeat
    
    AC_Voltage(:,count+1)=V_STEP(1:length(V_index),1);
    AC_Current(:,count+1)=I_STEP;
    
    
end
Vo=AC_Voltage;
Io=AC_Current;
T=t;
V_indexo=V_index;
I_indexo=I_index;
end