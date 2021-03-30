function [lickFlag relLickTime] = detectLickMarch21 (lickSenseDuration, niCardSession)


%This function detects the lick on the  spout through the output of the Janelia lick
%sensor board that is connected to a "Static"** ni-card digital inputs. 
%INPUTS. lickSenseDuration: period during which the lick sensor is checked.
%niCardSession: the NI-card session that is created bofore the function,
%and it includes three ports (right lick sensor, left lick sensor and the
%lever sensor) but we just monitor the port 1 which is connected to the
%right lick sensor
%call and indicate the NI-card channel that the sensor is connected to.

%OUTPUTS. lickFlag: 1 if the lick is detected on the spout during lickSenseDuration
%**Static digital input chnnel means that the channel can't be monitored continously
%through the background and forground commands and should be monitored with
%SingleScan commands.


startFlag = 0;
lickFlag = 0;

while (1)

    digitalInput = inputSingleScan(niCardSession);
    inputTime = GetSecs();
    port1 = digitalInput(1);
%     port2 = digitalInput(2);
    

    if ~startFlag
        startTime = inputTime;
        startFlag = 1;
    end

    if (inputTime>(startTime+lickSenseDuration))
        break;
    end



    if port1
        lickFlag = 1;

        relLickTime = inputTime - startTime;
        return

    end


end
relLickTime = inputTime - startTime;

end