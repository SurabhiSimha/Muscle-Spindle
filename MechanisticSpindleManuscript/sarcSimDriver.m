function [hsB,dataB,hsC,dataC] = sarcSimDriver(t,delta_f_activated,delta_cdl)
time_step = t(2)-t(1);
CHECK_SLACK = 1;  

% Make a half-sarcomere
hsB = halfSarcBag();
hsC = halfSarcChain();

% Loop through the time-steps
for a = 1
    for i=1:numel(t)
        
        %This applies an a priori condition in which slack is likely
        %(Speeds up simulations by factor of ~3)
        if hsB.hs_force <= 1e3
            CHECK_SLACK = 1;
        elseif hsC.hs_force <= 1e3
            CHECK_SLACK = 1;
        else
            CHECK_SLACK = 0;
        end
        
        if CHECK_SLACK
%             slack_mode = 0;
            %CHECK IF SLACK
            %ADJUST HS LENGTH SO FORCE IS BACK TO 0
%             hs.hs_force = 0;
            xB = find_hsl_from_force(hsB,0.0);
            xC = find_hsl_from_force(hsC,0.0);
          
            
            if hsB.cmd_length<xB
                hsB.hs_force = 0;
                slack_modeB = 1;
%                 new_length = max(x,hs.cmd_length);
                adj_lengthB = xB - hsB.hs_length;
                hsB.forwardStep(0.0,adj_lengthB,0,0,0,1);
            else
                slack_modeB = 0;
            end
            
            if hsC.cmd_length<xC
                hsC.hs_force = 0;
                slack_modeC = 1;
                adj_lengthS = xC - hsC.hs_length;
                hsC.forwardStep(0.0,adj_lengthS,0,0,0,1);
            else
                slack_modeC = 0;
            end
            
            if slack_modeB
               %CROSS BRIDGE EVOLUTION TAKES UP SLACK  
               % Any cb cycling here must be applied to shortening
               % against zero load.
                
               % First, we evolve the distribution
                hsB.forwardStep(time_step,0.0,0.0,delta_f_activated(a,i),1,0)
                
               % Next, we iteratively search for the sarcomere length
               % that would give us zero load
                xB = find_hsl_from_force(hsB,0);
                
               % We then compute the new hs length applied to the
               % sarcomere based on whether the command length is now
               % greater than the sarcomere length (back into
               % length-control) or not (still in isotonic mode). The
               % length adjustment is the the calculated new length - the 
               % current measurement of hs length. 
               
                new_lengthD = max(xB,hsB.cmd_length);
                adj_lengthB = new_lengthD - hsB.hs_length;
                
               % Finally, we shift the distribution by the adjusted length
                hsB.forwardStep(time_step,adj_lengthB,delta_cdl(a,i),0,0,1);
                
            else %length control
                delta_hslB = delta_cdl(a,i);
                hsB.forwardStep(time_step,delta_hslB,delta_cdl(a,i),delta_f_activated(a,i),1,1);
            end
            
            if slack_modeC
                %CROSS BRIDGE EVOLUTION TAKES UP SLACK
                % Any cb cycling here must be applied to shortening
                % against zero load.
                
                % First, we evolve the distribution
                hsC.forwardStep(time_step,0.0,0.0,delta_f_activated(a,i),1,0)
                
                % Next, we iteratively search for the sarcomere length
                % that would give us zero load
                xC = find_hsl_from_force(hsC,0);
                
                % We then compute the new hs length applied to the
                % sarcomere based on whether the command length is now
                % greater than the sarcomere length (back into
                % length-control) or not (still in isotonic mode). The
                % length adjustment is the the calculated new length - the
                % current measurement of hs length.
                
                new_lengthS = max(xC,hsC.cmd_length);
                adj_lengthS = new_lengthS - hsC.hs_length;
                
                % Finally, we shift the distribution by the adjusted length
                hsC.forwardStep(time_step,adj_lengthS,delta_cdl(a,i),0,0,1);
                
            else %length control
                delta_hslS = delta_cdl(a,i);
                hsC.forwardStep(time_step,delta_hslS,delta_cdl(a,i),delta_f_activated(a,i),1,1);
            end
            
        else %Assume length control
            delta_hslB = delta_cdl(a,i);
            hsB.forwardStep(time_step,delta_hslB,delta_cdl(a,i),delta_f_activated(a,i),1,1);
            delta_hslC = delta_cdl(a,i);
            hsC.forwardStep(time_step,delta_hslC,delta_cdl(a,i),delta_f_activated(a,i),1,1);
        end
        
        % Store data
        
        dataB(a).f_activated(i) = hsB.f_activated;
        dataB(a).f_bound(i) = hsB.f_bound;
        dataB(a).f_overlap(i) = hsB.f_overlap;
        dataB(a).cb_force(i) = hsB.cb_force;
        
        
        dataB(a).passive_force(i) = hsB.passive_force;
        dataB(a).hs_force(i) = hsB.hs_force;
        dataB(a).hs_length(i) = hsB.hs_length;
        dataB(a).cmd_length(i) = hsB.cmd_length;
        dataB(a).bin_pops(:,i) = hsB.bin_pops;
        dataB(a).no_detached(i) = hsB.no_detached;
        
        dataC(a).f_activated(i) = hsC.f_activated;
        dataC(a).f_bound(i) = hsC.f_bound;
        dataC(a).f_overlap(i) = hsC.f_overlap;
        dataC(a).cb_force(i) = hsC.cb_force;
        
        
        dataC(a).passive_force(i) = hsC.passive_force;
        dataC(a).hs_force(i) = hsC.hs_force;
        dataC(a).hs_length(i) = hsC.hs_length;
        dataC(a).cmd_length(i) = hsC.cmd_length;
        dataC(a).bin_pops(:,i) = hsC.bin_pops;
        dataC(a).no_detached(i) = hsC.no_detached;
        
    end
    dataB(a).t = t;
    dataC(a).t = t;
end

