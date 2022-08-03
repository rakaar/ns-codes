close all
spike_train = squeeze(spikes(1,1,1,250:400));
tau_syn = 10;
kernel_kt = [0 exp(-[0:length(spike_train)-1])./tau_syn];
conv_vector = conv(spike_train, kernel_kt);
my_conv = zeros(1, length(spike_train));

% nearest_spike_time = 0;
% for i=1:length(spike_train)
%     if spike_train(i) == 1
%         nearest_spike_time = i;
%     end
%     if nearest_spike_time ~= 0 && i <= length(spike_train)-1
%         my_conv(i+1)  = 0.1*exp(-(i-nearest_spike_time));
%     end
% end

for i=1:length(spike_train)
    if i >=  11
        x = conv(spike_train(i-10:i),kernel_kt);
        my_conv(1,i) = x(11);
    else
         x = conv(spike_train(1:i),kernel_kt);
         my_conv(1,i) = x(11);
    end
end


figure
hold on
stem(spike_train/4)
plot(conv_vector)
plot(my_conv)

legend('spike','matlab conv', ',myconv')
hold off
grid

figure
    c = conv_vector(1,1:length(spike_train));
    plot(c-my_conv)
    title('diff')
 grid

figure
    plot(my_conv)
    title('myconv')
grid