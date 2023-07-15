%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% Veerendra P Gadekar 
% Email: ic36871@imail.iitm.ac.in; gpveerendra09@gmail.com

% script to run gapolyfitn iteratively 
% with different maxterm and maxpower parameters

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% please provide the path to gapolyfitn packages below
addpath /path/to/gapolyfitn/gapolyfitndir/
savepath
addpath /path/to/gapolyfitn/gfit2/
savepath
addpath /path/to/gapolyfitn/randMat/
savepath
addpath /path/to/gapolyfitn/Polyfitn/
savepath
addpath /path/to/gapolyfitn/GAToolbox/UoS-CODeM-GA-Toolbox-8ca5c5c/
savepath
addpath /path/to/gapolyfitn/Polyfitn/PolyfitnTools/
savepath
addpath /path/to/gapolyfitn/multicore/
savepath

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% log
% term 2 power 2

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_T2_P2_itr',  v ,'.txt');
    eval(file)
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_bpd_T2_P2_itr',  v ,'.txt');
    eval(file)
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_bpd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ofd_T2_P2_itr',  v ,'.txt');
    eval(file)
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ofd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_hp_T2_P2_itr',  v ,'.txt');
    eval(file)
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_hp.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ap_T2_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ap.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_fl_T2_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_fl.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% term 3 power 2
for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_T3_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_bpd_T3_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_bpd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ofd_T3_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ofd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_hp_T3_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_hp.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ap_T3_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ap.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_fl_T3_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_fl.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% term 4 power 2
for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_T4_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_bpd_T4_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_bpd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end


for v = 1:10
    eval(file)    
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ofd_T4_P2_itr',  v ,'.txt');
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ofd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_hp_T4_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_hp.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end


for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ap_T4_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ap.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_fl_T4_P2_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_fl.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 2
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


% term 2 power 3
for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_T2_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_bpd_T2_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_bpd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end


for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ofd_T2_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ofd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_hp_T2_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_hp.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ap_T2_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ap.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_fl_T2_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_fl.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 2
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% term 3 power 3
for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_T3_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end


for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_bpd_T3_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_bpd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ofd_T3_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ofd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_hp_T3_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_hp.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:25
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt_a/classic_log_sqrt_ap_T3_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ap.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_fl_T3_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_fl.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 3
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% term 4 power 3
for v = 4:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_T4_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_bpd_T4_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_bpd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ofd_T4_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ofd.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_hp_T4_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_hp.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_ap_T4_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_ap.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

for v = 1:10
    file = sprintf('%s%d','diary ./results/output_gapolyfitn_matlab_classic_log_sqrt/classic_log_sqrt_fl_T4_P3_itr',  v ,'.txt');
    eval(file)    
    objtrain_log   = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_log_ga.xlsx")
    objtrain_log(:, 1) = [];
    paramtrain_log = readtable("./data/input_for_gapolyfitn_matlab/dbscan_train_data_classic_log_sqrt_fl.xlsx")
    paramtrain_log(:, 1) = [];

    indepvar = paramtrain_log.Variables
    depvar = objtrain_log.Variables
    maxTerms = 4
    maxPower = 3
    options.VERBOSE = 1;
    options.MAXGEN = 1000;
    options.DOSAVE = 0;
    [polymodel, Best, IndAll] = gapolyfitn(indepvar, depvar, maxTerms, maxPower, options)
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

