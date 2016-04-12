clear all;clc;

global waiting_time;
global AWE;
global NR_Maxnum;
global ckt_name;

global Gmin_ini;
global Gminlimit;

global Linear_Zone;

%warning('off');
%% AC solver setting
%AWE=false;
AWE=false;
%% NR Convergency setting
NR_Maxnum=50;
Linear_Zone=1e-4;
%% Gmin setting
Gmin_ini=1e-10; % gmin initial conductance
Gminlimit=1e-12;% gmin contraint 1




%% set NR loop number
waiting_time=1;
%% load the file
[FileName,PathName] = uigetfile('*.ckt','Select the curcuit-file');

ckt_name=FileName;
if isequal(FileName,0)
    disp('User selected Cancel')
else
    file_path=fullfile(PathName,FileName);
    
    
    %% ckt file scan
    % Low-Level File I/O
    fprintf('Run Programme:scanning the file....\n');
    [R,C,L,V,I,MOS,MOS_MODEL,DC,TR,AC,OUT_REQ] = ScanCKT(file_path);
    [matrix_dim,node_num]=Mat_size(R,C,L,V,I,MOS);
    
    pause(waiting_time);
    fprintf('\t\tfinished!\n');
    
    anyDC=DC;
    anyTR=~isempty(TR{1,1});
    anyAC=~isempty(AC{1,1});
    
    %% Circuit Calulation
    if anyDC
        
        %%IF DC analysis
        Null_1=zeros(matrix_dim,matrix_dim);
        Null_2=zeros(matrix_dim,1);
        fprintf('Run Programme:Caculating the DC steady response....\n');
        
        [V_DC,V_index,I_DC,I_index]=DCanlysis(R,C,L,V,I,MOS,MOS_MODEL,Null_1,Null_2);
        if ~anyTR
            Outstream(V_DC,V_index,I_DC,I_index,Null_1,Null_1,'DC');
        end
    end
    
    %% TR analysis
    if anyTR
        fprintf('Run Programme:Caculating the Transient response....\n');
        [V_TR,V_index,I_TR,I_index,T]=TRanlysis(R,C,L,V,I,MOS,MOS_MODEL,TR);
        pause(waiting_time);
        fprintf('\t\tfinished!\n');
        % ploting according to the request
        Outstream(V_TR,V_index,I_TR,I_index,T,OUT_REQ,'TR');
        if anyDC
            
            Outstream(V_DC,V_index,I_DC,I_index,Null_1,Null_1,'DC');
        end
        
    end
    
    %% Rf anlysis
    if anyAC
        fprintf('Run Programme:Caculating the AC response....\n');
        [V_AC,V_index,I_AC,I_index,F]=ACanlysis(R,C,L,V,I,MOS,MOS_MODEL,AC);
        pause(waiting_time);
        fprintf('\t\tfinished!\n');
        % ploting according to the request
        Outstream(V_AC,V_index,I_AC,I_index,F,OUT_REQ,'AC');
    end
warning('on');    
end
