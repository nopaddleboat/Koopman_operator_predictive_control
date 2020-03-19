# Koopman_operator_predictive_control
#This project contains simulation codes for the paper "Discrete System Linearization using Koopman Operators for Predictive Control and Its Applications"

data:
net12.mat    :   RNN model for the PEA system
xstar_np.mat   : LME model mentioned in [1]

codes:
prediction_accuracy  : compare the accuracy of linearization with Taylor series and Koopman operators
predictive_control_koopman1  : predictive control, liearization using Koopman operators, one linear model 
predictive_control_koopman2  : predictive control, liearization using Koopman operators, two linear models 
predictive_control_taylor    : predictive control, liearization using Taylor series 




[1] Xie, Shengwen, and Juan Ren. "Recurrent-neural-network-based Predictive Control of Piezo Actuators for Precision Trajectory Tracking." 2019 American Control Conference (ACC). IEEE, 2019.
