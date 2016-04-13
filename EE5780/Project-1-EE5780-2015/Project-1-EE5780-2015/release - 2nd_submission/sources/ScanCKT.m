function [R,C,L,V,I,MOS,MOS_MODEL,DC,TR,AC,OUT_REQ] = ScanCKT(file_path)

R={'Name','From','To','Value';...
    'NA',0,0,'NA',};
C={'Name','From','To','Value';...
    'NA',0,0,'NA',};
L={'Name','From','To','Value';...
    'NA',0,0,'NA',};
V={'Name','From','To','TYPE','Value';...
    'NA',0,0,'NA','NA'};
I={'Name','From','To','TYPE','Value';...
    'NA',0,0,'NA','NA'};
MOS={'Name','Drain','Gate','Source','TYPE','Width','Length','charicteristic';...
    'NA',0,0,0,'NA','NA','NA','NA'};
MOS_MODEL={'charicteristic','TH','MU','COX','LAMBDA';...
    'NA',0,0,0,0};
OUT_REQ={'Type','node','branch';...
    0,'NA','NA'};

DC=false;
TR=cell(1,3);
AC=cell(1,4);

cktfile=fopen(file_path);                               % open the ckt file           
% initiate the cells and variables          
R_i=2;
C_i=2;
L_i=2;
V_i=2;
I_i=2;
MOS_i=2;
MODEL_i=2;
OUT_i=2;
while ~feof(cktfile)                                      % until the end of file
    
    tline=fgetl(cktfile);                                 % Read line by line
    
    if length(tline)==0
        
    else
        if strcmp(tline(1),'R')||strcmp(tline(1),'r')        % Read the first letter
            R(R_i,:)=textscan(tline, '%s%d16%d16%f');           % scan this line
            R_i=R_i+1;
            continue                                         % continue scan
        elseif strcmp(tline(1),'C')||strcmp(tline(1),'c')       % Read the first letter
            C(C_i,:)=textscan(tline, '%s%d16%d16%f');           % scan this line
            C_i=C_i+1;
            continue                                         % continue scan
        elseif strcmp(tline(1),'L')||strcmp(tline(1),'l')        % Read the first letter
            L(L_i,:)=textscan(tline, '%s%d16%d16%f');           % scan this line
            L_i=L_i+1;
            continue                                         % continue scan
        elseif strcmp(tline(1),'V')||strcmp(tline(1),'v')      % Read the first letter
            temp=textscan(tline, '%s%d16%d16%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f');% scan this line
            V(V_i,1:length(temp))=temp;
            V_i=V_i+1;
            continue                                         % continue scan
        elseif strcmp(tline(1),'I')||strcmp(tline(1),'i')       % Read the first letter
            temp=textscan(tline, '%s%d16%d16%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f');% scan this line
            I(I_i,1:length(temp))=temp;
            I_i=I_i+1;
            continue                                         % continue scan
        elseif strcmp(tline(1),'M')||strcmp(tline(1),'m')       % Read the first letter
            MOS(MOS_i,:)=textscan(tline, '%s%d16%d16%d16%s%f%f%d16');           % scan this line
            MOS_i=MOS_i+1;
            continue                                         % continue scan
        elseif strcmp(tline(1),'.')
            
            if strcmp(tline(2:3),'DC')||strcmp(tline(2:3),'dc')
                DC=true;
                
            end
            if strcmp(tline(2:3),'TR')||strcmp(tline(2:3),'tr')       % Read the first letter
                
                TR(1,2:3)=textscan(tline, '%*s%f%f');           % scan this line
                TR{1,1}=true;
            end
            if strcmp(tline(2:3),'AC')||strcmp(tline(2:3),'dc')       % Read the first letter
                AC(1,2:4)=textscan(tline, '%*s%f%f%f');           % scan this line
                AC{1,1}=true;
                V(1,1:5)={'Name','From','To','Value(s domain)',[]};
            end
            
            if(length(tline)>3)
                if strcmp(tline(2:6),'MODEL')||strcmp(tline(2:6),'model')       % Read the first letter
                    MOS_MODEL(MODEL_i,:)=textscan(tline, '%*s%d16%*s%f%*s%f%*s%f%*s%f');           % scan this line
                    MODEL_i=MODEL_i+1;
                end
                
                if strcmp(tline(2:8),'PRINTNV')||strcmp(tline(2:6),'printnv')       % Read the first letter
                    OUT_REQ{OUT_i,1}=1;
                    OUT_REQ(OUT_i,2)=textscan(tline, '%*s%d16');           % scan this line
                    OUT_REQ{OUT_i,3}='NA';   
                    OUT_i=OUT_i+1;
                end
                
                if strcmp(tline(2:8),'PRINTBI')||strcmp(tline(2:6),'printbi')       % Read the first letter
                    OUT_REQ{OUT_i,1}=2;
                    OUT_REQ(OUT_i,3)=textscan(tline, '%*s%s');           % scan this line
                    OUT_REQ{OUT_i,2}='NA';   
                    OUT_i=OUT_i+1;
                end
                
                if strcmp(tline(2:7),'PLOTNV')||strcmp(tline(2:6),'plotnv')       % Read the first letter
                    OUT_REQ{OUT_i,1}=3;
                    OUT_REQ(OUT_i,2)=textscan(tline, '%*s%d16');           % scan this line
                    OUT_REQ{OUT_i,3}='NA';   
                    OUT_i=OUT_i+1;
                end
                
                if strcmp(tline(2:7),'PLOTBI')||strcmp(tline(2:6),'plotbi')       % Read the first letter
                    OUT_REQ{OUT_i,1}=4;
                    OUT_REQ(OUT_i,3)=textscan(tline, '%*s%s');           % scan this line
                    OUT_REQ{OUT_i,2}='NA';   
                    OUT_i=OUT_i+1;
                end
                
                continue % continue scan
            end
            
        end
    end
    
end
end

