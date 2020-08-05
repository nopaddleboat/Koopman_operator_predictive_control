 Koopman_operator_predictive_control
* This project contains simulation codes for the paper "Discrete System Linearization using Koopman Operators for Predictive Control and Its Applications"[2].
* All the codes are Not commented and NOT optimized, but they should be easy to follow with the paper aside.
* Feel free to contact me through swxie@outlook.com if you have any questions regarding the code and paper.

Data:\
* net12.mat    :   RNN model for the PEA system
* xstar_np.mat   : LME model mentioned in [1]
* These two models are generated based on the real system input-output data with the method in [1].

Codes:\
* prediction_accuracy  : compare the accuracy of linearization with Taylor series and Koopman operators
* predictive_control_koopman1  : predictive control, linearization using Koopman operators, one linear model 
* predictive_control_koopman2  : predictive control, linearization using Koopman operators, two linear models 
* predictive_control_taylor    : predictive control, linearization using Taylor series 

#Note that the results generated from these codes may be different from the results reported in the paper, this is due to the characteristics of randomness in the method.

[1] Xie, Shengwen, and Juan Ren. "Recurrent-neural-network-based Predictive Control of Piezo Actuators for Trajectory Tracking." IEEE/ASME Transactions on Mechatronics (2019).\
[2] Xie, Shengwen, and Juan Ren. "Linearization of Recurrent-neural-network-based models for Predictive Control of Nano-positioning Systems using Data-driven Koopman Operators" IEEE Access (2020). DOI:10.1109/ACCESS.2020.3013935.
