clear all

niDevName = 'Dev1';
niPortLine1 = 'port0/line0';


rewardStepMotorCtl1 = daq.createSession('ni');
rewardStepMotorCtl1.addDigitalChannel(niDevName,niPortLine1,'OutputOnly');

rewardStepMotorEnable = daq.createSession('ni');
motorEnablePortLine = 'port0/line5';
rewardStepMotorEnable.addDigitalChannel(niDevName,motorEnablePortLine,'OutputOnly');
% rewardStepMotorEnable.outputSingleScan(1);
rewardStepMotorEnable.outputSingleScan(0);

manualReward = 0;
rewardVol = 4; 

deliverRewardMarch21(rewardVol,5,rewardStepMotorCtl1);

manualReward = manualReward + rewardVol;
disp(['Total Reward: ', num2str(manualReward)]);

rewardStepMotorEnable.outputSingleScan(1);

rewardStepMotorCtl1.release();