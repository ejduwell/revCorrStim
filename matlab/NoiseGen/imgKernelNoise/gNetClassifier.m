function [labelOut,labelConf,scores] = gNetClassifier(imgIn,netStr)

netCmd=strcat("net = ",netStr);
eval(netCmd);

inputSize = net.Layers(1).InputSize;

I = imresize(imgIn,inputSize(1:2));

[labelOut,scores] = classify(net,I);

labelConf = 100*scores(classNames == labelOut);

end