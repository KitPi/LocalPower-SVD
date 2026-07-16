

PowerIteration = [8,4,2,8,4,2,8,4,2];

NumThread = [8,8,8,4,4,4,2,2,2];

TotalTime = [9.06261133333333,8.605155,8.14293933333333,11.2341593333333,10.6867413333333,9.69740633333333,14.9222093333333,14.0515413333333,12.219839];
%ProcessTime = [1.70064329166667,1.45637983333333,0.979620791666667,3.43791783333333,2.92534958333333,1.93053083333333,6.8825485,5.90270366666667,3.92319266666667];
ProcessTime = [57.0283173333333,48.9525693333333,30.9574983333333,28.369276,24.5985803333333,16.342885,15.0615596666667,13.018988,8.49762333333333]

Error = [4.233081,3.697486,3.62206366666667,3.61976,3.47782166666667,4.43357633333333,4.39731933333333,3.74245,3.790614];



[X, Y] = meshgrid(unique(PowerIteration), unique(NumThread));

t_time = griddata(PowerIteration, NumThread, TotalTime, X, Y, 'linear');
p_time = griddata(PowerIteration, NumThread, ProcessTime./TotalTime, X, Y, 'linear');
error = griddata(PowerIteration, NumThread, Error, X, Y, 'linear');


surf(X,Y,t_time);
xlabel('Power Iterations');
ylabel('Num Processes');
zlabel('Total Time');
title('3D Surface Plot of Total Time vs Power Iterations and Num Processes'); grid on;

figure;
surf(X,Y,p_time);
xlabel('Power Iterations');
ylabel('Num Processes');
zlabel('Useful Time');
title('3D Surface Plot of Useful Time vs Power Iterations and Num Processes'); grid on;

figure;
surf(X,Y,error);
xlabel('Power Iterations');
ylabel('Num Processes');
zlabel('Error');
title('3D Surface Plot of Error vs Power Iterations and Num Processes'); grid on;


TotalTime = [12.496779,12.558154,12.509493,12.1502836666667,12.44515,12.639199,12.6273286666667,12.5934426666667,12.6201093333333];
%ThreadTime = [7.12853966666667,6.11907116666667,3.86968729166667,7.092319,6.14964508333333,4.08572125,7.53077983333333,6.509494,4.24881166666667];
ThreadTime = [13.6051463333333,11.6510386666667,7.83696633333333,13.7516713333333,11.7013983333333,7.72212333333333,13.765097,11.8054073333333,7.84638533333333]

Error = [3.580631,3.17998366666667,3.81608433333333,4.53182166666667,4.66895966666667,4.28700566666667,4.201654,4.30027966666667,4.42504833333333];

t_time = griddata(PowerIteration, NumThread, TotalTime, X, Y, 'linear');
p_time = griddata(PowerIteration, NumThread, ThreadTime./TotalTime, X, Y, 'linear');
error = griddata(PowerIteration, NumThread, Error, X, Y, 'linear');

%figure;
%surf(X,Y,t_time);
%xlabel('Power Iterations');
%ylabel('Num Threads');
%zlabel('Total Time');
%title('3D Surface Plot of Total Time vs Power Iterations and Num Threads'); grid on;

figure;
surf(X,Y,p_time);
xlabel('Power Iterations');
ylabel('Num Threads');
zlabel('Useful Time');
title('3D Surface Plot of Useful Time vs Power Iterations and Num Threads'); grid on;

%figure;
%surf(X,Y,error);
%xlabel('Power Iterations');
%ylabel('Num Threads');
%zlabel('Error');
%title('3D Surface Plot of Error vs Power Iterations and Num Threads'); grid on;
