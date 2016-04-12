function Outstream(V_array,V_index,I_array,I_index,Xaxle,OUT_REQ,Solvermode)
global ckt_name;
global waiting_time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% type 1. PRINTNV print the node voltage
% type 2. PRINTBI print the element current
% type 3. PLOTNV plot the node voltage
% type 4. PLOTBI plot the element current

plot_count=1;
r2d=180/pi;
fprintf('generating result\n');
pause(waiting_time);


if strcmp(Solvermode,'DC') || strcmp(Solvermode,'dc')
    fprintf('DC response result\n');
     pause(waiting_time);
    for n=2:length(V_index)
         titlestr=['Steady voltage of node '  num2str(V_index(n)) ' in ckt <' ckt_name '> is '];
         fprintf(titlestr);
         fprintf(num2str(V_array(n)));
         fprintf(' V\n');
    end
    for n=1:length(I_index)
         element=I_index(n,:);
         titlestr=['Steady current of element '  num2str(element) ' in ckt <' ckt_name '> is '];
         fprintf(titlestr);
         fprintf(num2str(I_array(n)));
         fprintf(' A\n');
    end
elseif strcmp(Solvermode,'TR') || strcmp(Solvermode,'tr')
     fprintf('Transient response result\n');
     pause(waiting_time);
    for n=2:length(OUT_REQ)
        
        
        switch OUT_REQ{n,1}
            
            case 1
                node_index=find(V_index==OUT_REQ{n,2});
                last_index=length(Xaxle);
                titlestr=['Steady voltage of node '  num2str(OUT_REQ{n,2}) ' in ckt <' ckt_name '> is '];
                fprintf(titlestr);
                fprintf(num2str(V_array(node_index,last_index)));
                fprintf(' V\n');

            case 2
                ele_index=Loc(I_index,cell2mat(OUT_REQ{n,3}));
                last_index=length(Xaxle);
                element=cell2mat(OUT_REQ{n,3});
                titlestr=['Steady current of node '  num2str(element) ' in ckt <' ckt_name '> is '];
                fprintf(titlestr);
                fprintf(num2str(I_array(ele_index,last_index)));
                fprintf(' A\n');
            case 3
                index=find(V_index==OUT_REQ{n,2});
                h=figure(plot_count);
                plot(Xaxle,V_array(index,:));
                xlabel('time(s)');
                ylabel('node voltage(V)');
                titlestr=['Voltage of node '  num2str(OUT_REQ{n,2}) ' in ckt <' ckt_name '>'];
                set(h,'name',titlestr,'NumberTitle','off');
                plot_count=plot_count+1;
            case 4
                index=Loc(I_index,cell2mat(OUT_REQ{n,3}));
                h=figure(plot_count);
                plot(Xaxle,I_array(index,:));
                xlabel('time(s)');
                ylabel('node current(A)');
                
                element=cell2mat(OUT_REQ{n,3});
                titlestr=['Current of element '  num2str(element) ' in ckt <' ckt_name '>'];
                set(h,'name',titlestr,'NumberTitle','off');
                plot_count=plot_count+1;
            otherwise
        end
    end
    
    
    
elseif strcmp(Solvermode,'AC') || strcmp(Solvermode,'ac')
     fprintf('AC response result\n');
     pause(waiting_time);
    for n=2:length(OUT_REQ)
        
        
        switch OUT_REQ{n,1}
            case 3
                index=find(V_index==OUT_REQ{n,2});
                h=figure(plot_count);
                subplot(2,1,1)
                plot(Xaxle,20*log10(abs(V_array(index,:))));
                xlabel('angular frequent (rad/s) in logarithmic scale');
                ylabel('Magnitude(dB)');
                titlestr=['Bode diagram of node '  num2str(OUT_REQ{n,2}) ' in ckt <' ckt_name '>'];
                set(h,'name',titlestr,'NumberTitle','off');
                subplot(2,1,2)
                plot(Xaxle,angle(V_array(index,:))*r2d);
                xlabel('angular frequent (rad/s) in logarithmic scale');
                ylabel('Phase(deg)');

                %%bode: magnitude and phase 
                plot_count=plot_count+1;
        end
    end
    
end
    pause(waiting_time);
    fprintf('\t\tfinished!\n');
end
