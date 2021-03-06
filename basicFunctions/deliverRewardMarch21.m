function [totalSteps,exactRewardValue] = deliverReward (volume_uL, syringeSize_mL, niCardSession)

% We will use this function to deliver the reward by controling the step motor via the NI-card
% volume_uL: amount of reward to be delivered in micro liters,
% syringeSize_mL: size of the syringe in mili liters, niCardSession: the
% handle to the ni-card session with one digital output that is created
% before calling this function (creating session each time in the function
% would take a variable time which is not desirable). 
% totalSteps: number of steps to deliver the volume_uL that is depend on
% the syringeSize_mL. exactRewardValue: the exact theoretical value of
% reward based on the total number of steps that should be very close to
% the volume_uL
% This function was rewrited based on the written function by Danniella for Arduino! the
% timing parts of that function are ignored (the timing would depend on the
% velocity of the motor which we don't have any reliable way to measure and since it happens at the end of each trial, we decided that we don't need the exact end time of reward)
% 

% enableSession.outputSingleScan(0);
% 
% delay = 0.1;
% 
% startTime = GetSecs();
% while(GetSecs - startTime) < delay
%     ;
% end


if (syringeSize_mL == 5) 
    diameter_mm = 12.06; %in mm
elseif (syringeSize_mL == 10)
    diameter_mm = 14.5; %in mm
else
    print("didn't recognize a valid syringe size. available sizes '5' or '10.'"); 
    return; 
end


% test changes in github




% // determine vol per revolution, area of small cylinder with h=0.8mm
%   // 0.8mm length per thread. 1thread=1cycle. 1 like=1prayer.

volPerRevolution_uL = 0.8 * ( diameter_mm/2 )*( diameter_mm/2 ) * pi ; 

% // determine how many revolutions needed for the desired volume
correctionFactor = 1;%MS1=MS3=0 and MS2=1 %2.25; %MS1, MS2 and MS3 are off for the new soundless board tested on box1 3/10/2020 
howManyRevolutions = correctionFactor*volume_uL / volPerRevolution_uL ;

%   // determine total steps needed to reach desired revolutions, @200
%   steps/revolution (step motor specification)
%   // use *4 as a multiplier because it's operating at 1/4 microstep mode (MS1=0,MS2=1,MS3=0).
%   // round to nearest int because totalSteps is unsigned long
totalSteps = round(200 * howManyRevolutions * 4);

exactRewardValue = totalSteps*volPerRevolution_uL/800;

%   // determine shortest delivery duration, total steps * 2 ms per step.
%   (where does 1 ms come from?)
%   // minimum 1 ms in high, 1 ms in low for the shortest possible step function.
% minimumDeliveryDuration_ms = totalSteps*2; 

%   // make sure delivery duration the user wants is long enough
% if (local_deliveryDuration_ms < minimumDeliveryDuration_ms)
%     print("duration too low. duration needs to be >");
%     print(minimumDeliveryDuration_ms); 
%     print("with that diameter and reward volume.");
%     return;
% end

%   // determine duration of each step for the timer oscillate function
% stepDuration_ms = local_deliveryDuration_ms / totalSteps;

%to enable the pump, moved to the begining and the end of each session
% niCardSession.outputSingleScan(0)

% we can change it based on how fast the reward seems to be delivered and
% how load is the pump
waitBetweenPulses = 0.001; % 2.5 ms off and 2.5 ms on periods, based on what seems to be used with Soyoun to deliver 4 microL in 200 ms

for loop=1:totalSteps
    
    startTime = GetSecs;
    while (GetSecs - startTime) < waitBetweenPulses
        ;
    end
    
    niCardSession.outputSingleScan(0)
    
    startTime = GetSecs;
    while (GetSecs - startTime) < waitBetweenPulses
        ;
    end
    
    niCardSession.outputSingleScan(1)
%     niCardSession.outputSingleScan(1)
%     niCardSession.outputSingleScan(1)
%     niCardSession.outputSingleScan(1)
%     niCardSession.outputSingleScan(1)
%     niCardSession.outputSingleScan(0)
%     niCardSession.outputSingleScan(0)
%     niCardSession.outputSingleScan(0)
%     niCardSession.outputSingleScan(0)
%     niCardSession.outputSingleScan(0)
    
end

%to disable the pump, moved to the begining and the end of each session
% enableSession.outputSingleScan(1);

end
