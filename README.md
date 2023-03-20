# HME

to run HME:  
Input .root file must be provided and specified in hme.cpp. <br>
cd hme/ <br>
./run.sh 

# DNN

- generate_code_template.cpp - generates .cpp file with main function and names of all branches declared, everything else should be written by you <br>
- fill_data.cpp - selects and calculates features required for DNN training and writes them to a .csv file <br>
- ML_test_neural_net.ipynb - notebook with the actual model

.cpp files require a .root file input (add location)

to run:  
cd DNN/ <br>
./run.sh fill
