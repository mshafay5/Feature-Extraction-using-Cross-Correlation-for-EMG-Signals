clear all
clc
close all
warning off
% Data Import and Plots 

T = readtable('emg_healthy.txt');
h_time = table2array(T(:,1));
h_amplitude = (table2array(T(:,2)));
subplot 311
plot(h_time(1:1000),h_amplitude(1:1000))
title("Healthy EMG")

T = readtable('emg_myopathy.txt');
m_time = table2array(T(:,1));
m_amplitude = (table2array(T(:,2)));
subplot 312
plot(m_time(1:1000),m_amplitude(1:1000))
title("Myopathy EMG")

T = readtable('emg_neuropathy.txt');
n_time = table2array(T(:,1));
n_amplitude = (table2array(T(:,2)));
subplot 313
plot(n_time(1:1000),n_amplitude(1:1000))
title("Neuropathy EMG")

% Data Scaling

hyScaled = (h_amplitude - nanmean(h_amplitude)) / nanstd(h_amplitude);
hyScaled(isnan(h_amplitude)) = 0;
myScaled = (m_amplitude - nanmean(m_amplitude)) / nanstd(m_amplitude);
myScaled(isnan(m_amplitude)) = 0;
nyScaled = (n_amplitude - nanmean(n_amplitude)) / nanstd(n_amplitude);
nyScaled(isnan(n_amplitude)) = 0;


%AutoCorelation of Healthy with Healthy and cross corelation of Healthy
%with Neuropathic and Healthy with Myopathic

[H_H_Corr_y,H_H_Corr_t] = xcorr(hyScaled((1:1000),1),hyScaled((1:1000),1));
[H_M_Corr_y,H_M_Corr_t] = xcorr(hyScaled((1:1000),1),myScaled((1:1000),1));
[H_N_Corr_y,H_N_Corr_t] = xcorr(hyScaled((1:1000),1),nyScaled((1:1000),1));

% Statistical Parameters

%Initializing Vectors and Arrays

H_H_Corr_y = zeros(1999,50);
H_M_Corr_y = zeros(1999,50);
H_N_Corr_y = zeros(1999,50);

H_H_Mean = zeros(1,50);
H_M_Mean = zeros(1,50);
H_N_Mean = zeros(1,50);

H_H_StdDev = zeros(1,50);
H_M_StdDev = zeros(1,50);
H_N_StdDev = zeros(1,50);

H_H_Pow = zeros(1,50);
H_M_Pow = zeros(1,50);
H_N_Pow = zeros(1,50);

H_H_Activity = zeros(1,50);
H_M_Activity = zeros(1,50);
H_N_Activity = zeros(1,50);

H_H_Mobility = zeros(1,50);
H_M_Mobility = zeros(1,50);
H_N_Mobility = zeros(1,50);

H_H_Complexity = zeros(1,50);
H_M_Complexity = zeros(1,50);
H_N_Complexity = zeros(1,50);

H_H_Kurtosis = zeros(1,50);
H_M_Kurtosis = zeros(1,50);
H_N_Kurtosis = zeros(1,50);

H_H_Skewness = zeros(1,50);
H_M_Skewness = zeros(1,50);
H_N_Skewness = zeros(1,50);

%Statistics of 50 sets are computed for each parameter

for i = 1:50
    H_H_Corr_y(:,i) = xcorr(hyScaled(((((i-1)*1000)+1):1000*i),1),hyScaled((1:1000),1));
    H_M_Corr_y(:,i) = xcorr(hyScaled(((((i-1)*1000)+1):1000*i),1),myScaled(((((i-1)*1000)+1):1000*i),1));
    H_N_Corr_y(:,i) = xcorr(hyScaled(((((i-1)*1000)+1):1000*i),1),nyScaled(((((i-1)*1000)+1):1000*i),1));
    
    H_H_Mean(1,i) = mean(H_H_Corr_y(:,i));
    H_M_Mean(1,i) = mean(H_M_Corr_y(:,i));
    H_N_Mean(1,i) = mean(H_N_Corr_y(:,i));
    
    H_H_StdDev(1,i) = std(H_H_Corr_y(:,i));
    H_M_StdDev(1,i) = std(H_M_Corr_y(:,i));
    H_N_StdDev(1,i) = std(H_N_Corr_y(:,i));
    
    H_H_Pow(1,i) = (1/numel(H_H_Corr_y(:,i)))*sum(abs(H_H_Corr_y(:,i)).^2);
    H_M_Pow(1,i) = (1/numel(H_M_Corr_y(:,i)))*sum(abs(H_M_Corr_y(:,i)).^2);
    H_N_Pow(1,i) = (1/numel(H_N_Corr_y(:,i)))*sum(abs(H_N_Corr_y(:,i)).^2);
    
    H_H_Activity(1,i) = var(H_H_Corr_y(:,i));
    H_M_Activity(1,i) = var(H_M_Corr_y(:,i));
    H_N_Activity(1,i) = var(H_N_Corr_y(:,i));
    
    H_H_Mobility(1,i) = var(diff(H_H_Corr_y(:,i))/H_H_Activity(:,i));
    H_M_Mobility(1,i) = var(diff(H_M_Corr_y(:,i))/H_M_Activity(:,i));
    H_N_Mobility(1,i) = var(diff(H_N_Corr_y(:,i))/H_N_Activity(:,i));
    
    H_H_Complexity(1,i) = (var(diff(diff(H_H_Corr_y(:,i)))/var(diff(H_H_Corr_y(:,i)))))/H_H_Mobility(:,i);
    H_M_Complexity(1,i) = (var(diff(diff(H_M_Corr_y(:,i)))/var(diff(H_M_Corr_y(:,i)))))/H_M_Mobility(:,i);
    H_N_Complexity(1,i) = (var(diff(diff(H_N_Corr_y(:,i)))/var(diff(H_N_Corr_y(:,i)))))/H_N_Mobility(:,i);
    
    H_H_Kurtosis(1,i) = kurtosis(H_H_Corr_y(:,i));
    H_M_Kurtosis(1,i) = kurtosis(H_M_Corr_y(:,i));
    H_N_Kurtosis(1,i) = kurtosis(H_N_Corr_y(:,i));
    
    H_H_Skewness(1,i) = skewness(H_H_Corr_y(:,i));
    H_M_Skewness(1,i) = skewness(H_M_Corr_y(:,i));
    H_N_Skewness(1,i) = skewness(H_N_Corr_y(:,i));


end

%Table is populated using known correlations of healthy-healthy, healthy -
%myopathic, and healthy - neuropathic
for i = 1:50
 statTable(i,:) = [H_H_Mean(1,i),H_H_StdDev(1,i),H_H_Pow(1,i),H_H_Activity(1,i),H_H_Mobility(1,i),H_H_Complexity(1,i),H_H_Kurtosis(1,i),H_H_Skewness(1,i),1];
end
for i = 1:50
 statTable(i+50,:) = [H_M_Mean(1,i),H_M_StdDev(1,i),H_M_Pow(1,i),H_M_Activity(1,i),H_M_Mobility(1,i),H_M_Complexity(1,i),H_M_Kurtosis(1,i),H_M_Skewness(1,i),2];
end
for i = 1:50
 statTable(i+100,:) = [H_N_Mean(1,i),H_N_StdDev(1,i),H_N_Pow(1,i),H_N_Activity(1,i),H_N_Mobility(1,i),H_N_Complexity(1,i),H_N_Kurtosis(1,i),H_N_Skewness(1,i),2];
end

% Classifying as healthy / diseased
%there is no label in testStat, But to include label in TrainStat, we used
%the following commands
healthyDiseaseTrainingStats = statTable([(1:37),(51:87),(101:137)],:);
healthyDiseaseTestStats = statTable([(38:50),(88:100),(138:150)],(1:8));

[healthyDiseaseClassifier,healthyDiseaseClassifierAccuracy] = HealthytrainClassifier(healthyDiseaseTrainingStats);
healthyDiseaseResults = healthyDiseaseClassifier.predictFcn(healthyDiseaseTestStats); % First 13 should be healthy, rest are either myopathic/neuropathic

%Plotting the confusion matrix for Healthy/Diseased Classifier
figure()
classifiedHealthyTest = healthyDiseaseResults;
HealthyDiseasedLabels(1:13) = 1; 
HealthyDiseasedLabels(14:39) = 2;
HealthyDiseasedLabels = transpose(HealthyDiseasedLabels);
confusionDiseased = confusionmat(HealthyDiseasedLabels,classifiedHealthyTest);
confusionchart(confusionDiseased)


for i = 1:39
    if healthyDiseaseResults(i,:) == 1
        healthyDiseaseResultsString(i,:) = "Healthy Patient";
    else
        healthyDiseaseResultsString(i,:) = "Diagnosed Patient";
    end
end

% Further classifying as myopathic / neuropathic

diseaseStatTable = statTable((51:150),:);
diseaseStatTable(1:50,9) = 1;

diseaseTrainingStats = diseaseStatTable([(1:37),(51:87)],:); 
diseaseTestStats = diseaseStatTable([(38:50),(88:100)],(1:8));

% Neuropathy & Myopathy Classifier
[diseaseClassifier,diseaseClassifierTrainingAccuracy] = DiseaseClassifierSVM(diseaseTrainingStats);
diseaseResults = diseaseClassifier.predictFcn(diseaseTestStats);

%Plotting confusion matrix for Neuropathy & Myopathy Classifier
figure()
NeuroMayoLabels(1:13) =1; 
NeuroMayoLabels(14:26) = 2;
NeuroMayoLabels = transpose(NeuroMayoLabels);
confusionDiseased = confusionmat(NeuroMayoLabels,diseaseResults);
confusionchart(confusionDiseased)

for i = 1:26
    if diseaseResults(i,:) == 1
        diseaseResultsString(i,:) = "Myopathic";
    else
        diseaseResultsString(i,:) = "Neuropathic";
    end
end

% Results

healthyDiseaseResultsString
diseaseResultsString