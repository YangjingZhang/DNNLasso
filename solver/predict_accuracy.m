%% precision = TP/(TP+FP)
%% recall = TP/(TP+FN)
%% fscore = 2*precision*recall/(precision+recall)
%% input symmetric n*n matrix true and estimate
function [precision,recall,fscore,mcc,acc] = predict_accuracy(true,estimate)
n = size(true,1);
tmp = ones(n,n);
tmp = triu(tmp,1);
ind = (tmp(:) > 0);
true = true(ind);
true = (true ~= 0);
estimate = estimate(ind);
estimate = (estimate ~= 0);
tp = sum((estimate ~= 0).*(true ~= 0));
fp = sum((estimate ~= 0).*(true == 0));
fn = sum((estimate == 0).*(true ~= 0));
tn = sum((estimate == 0).*(true == 0));
if tp + fp == 0 || tp + fn == 0
    precision = 0;
    recall = 0;
    fscore = 0;
    mcc = 0;
else
    precision = tp/(tp + fp);
    recall = tp/(tp + fn);
    fscore = 2*precision*recall/(precision + recall);
    mcc = (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn));
end
acc = (tp + tn)/(tp + fp + tn + fn);

end
