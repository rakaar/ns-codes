%%  Application of Loebel's model for SSA
% The following code has been adapted from that of Alex Loebel
% Scalars are denoted by letters as in the article, or else by lower-case words.
% Vectros and matrices are denoted by letters as in the article, or else by capitalized words.

% This file generates and brings to equilibrium a network that can be used
% to run Nevo's protocols

%close all
%clear all

all_plots   = 0; % Controls which plots to display
plot_mixing = 1; % Plots the best frequencies for all neurons
plot_TC     = 1; % Plots the shape of a typical tuning curve

find_FRA = 0;
id_crit = 0;

Sel_Columns = [8 9 10]; % Selects columns whose mean activity will be plotted
n_sing      = 1; % Number of single neurons that will be tracked in each selected column
inc_non_act = 1; % Includes non-active neurons in the single neurons tracked

graph_col = 'bgrcmyk'; % Colors for the graphs of the different cases
n_sel     = length(Sel_Columns);

Conds{1} = 'Low';
Conds{2} = 'High';
Conds{3} = 'Equal';
Conds{4} = 'Diverse Broad';
Conds{5} = 'Diverse Narrow';
Conds{6} = 'Deviant Alone F1';
Conds{7} = 'Deviant Alone F2';
Conds{8} = 'Fixed Low';
Conds{9} = 'Fixed High';

nev_cond_code{1} = 'L';
nev_cond_code{2} = 'H';
nev_cond_code{3} = 'E';
nev_cond_code{4} = 'DB';
nev_cond_code{5} = 'DN';
nev_cond_code{6} = 'DA1';
nev_cond_code{7} = 'DA2';
nev_cond_code{8} = 'FL';
nev_cond_code{9} = 'FH';

%% Parameters of the Model:
% Number of columns:
P = 21;

% Number of frequencies to which the network is sensitive
freq_stretch = 1; % Choose an integer >=1; this multiplies P to give M.
M = freq_stretch*P; % Thus, column best frequencies will be certain presentable frequencies.

ring_net = 0;

% Numbers of cells in each column (excitatory and inhibitory):
NE = 100;
NI = 100; 

% Controlling connection strengths:
factor   = 1; % multiplies all intra-column connections (controls spontaneous population spikes)
factor_1 = 1; % multiplies nearest-neighbor inter-column connections (aids spread of population spikes)
factor_2 = 1; % multiplies 2nd-nearest-neighbor inter-column connections (aids spread of population spikes)

% Background input drawn from a uniform distribution (neurons are indexed by input strength):
bg_low_E  = -9.9; % lowest background input (in Hz) to excitatory population
bg_high_E = 9.9; % highest background input (in Hz) to excitatory population
bg_low_I  = bg_low_E; % lowest background input (in Hz) to inhibitory population
bg_high_I = bg_high_E; % highest background input (in Hz) to inhibitory population
% *) In Loebel & Tsodyks 2002, the above values were reached by requiring a
% spontaneous activity of a few Hz. 

% Time constants:
tau_E     = 0.001; % excitatory neurons' time constant (in seconds)
tau_I     = 0.001; % inhibitory neurons' time constant (in seconds) 
tau_ref_E = 0.003; % tau refractory of excitatory neurons (in seconds)
tau_ref_I = 0.003; % tau refractory of inhibitory neurons (in seconds)
tau_rec   = 0.800; % recovery time constant of intracortical synapses (in seconds)
tau_rec_s = 0.300; % recovery time constant of sensory input synapses (in seconds)
% *) This value was chosen since Eli says that LFP recovers fully after about 1 s (tau_rec_s ~ 0.3 s); there should be
% articles about measurements in slices showing this recovery

U   = 0.5; % Portion of available fraction of resources that is utilized in response to an action potential
U_s = 0.7; % Same as U, only for the thalamo-cotical synapses that convey the sensory input

% Connection strengths:
Jee = 6*factor/NE;    % exc2exc
Jei = -4*factor/NI;   % inh2exc
Jie = 0.5*factor/NE;  % exc2inh
Jii = -0.5*factor/NI; % inh2inh

Jee_1 = 0.045*factor_1/NE;  % exc2exc, between neighboring columns
Jie_1 = 0.0035*factor_1/NE; % exc2inh, between neighboring columns
Jee_2 = 0.015*factor_2/NE;  % exc2exc, from one column to its 2nd neighbor
Jie_2 = 0.0015*factor_2/NE; % exc2inh, from one column to its 2nd neighbor

% Simulation conditions:
t_eq      = 5; % The time given to reach equilibrium (in seconds). It is important to allow enough time, or else the response to the first stimulus is distorted.
dt        = 0.0001; % Time-step (in seconds)
post_stim = 0.100; % This is defined here so the simulations keeps runnng for 2*post_stim after the last stimulus offset

%% Stimulus Parameters
A = 5; % Peak magnitude of the input (the magnitude felt by neurons most sensitive to that input)
A_max = 50; % This is used for tuning curves at high amplitudes (comes into effect if A_max > alpha)

ISI  = 0.3; % Inter-stimulus interval (in seconds)
duration = 0.050; % Total duration of each stimulus (in seconds)
ramp_dur = 0.005; % durations of the ramps at the beginning and end of each stimulus (in seconds)

n_stim = 100; % Total no. of stimuli (Best take a product of 10)
%t_prot = n_stim*(duration + ISI); % Total time of the oddball protocol

%% Localization of Sensory Input Effect:
lambda_c = 5; %0.25; % Base value of the lambda_s localization parameter
alpha    = 101; %2; % Magnitude threshold for increase in localization parameter
delta    = 5; % Slope of the increase in lambda_s for magnitudes greater than alpha (determines the degree of localization of an input at high sound levels)

% Option for asymmetric tuning curves:
delta_right = delta; 
delta_left  = delta;

%% Mixing cells with different frequency preferences in each column:
mix_tcs = 1; % 1 or 0, mixing or segregating the input to individual columns

part_2f = 16; % 1 over the fraction of cells in each column that are sensitive to the frequency of the column 2 columns forward is about 1/12
part_2b = 16; % Same, but for the column 2 columns back
part_1f = 8; % 1 over the fraction of cells in each column that are sensitive to the frequency of the next column is about 1/6
part_1b = 8; % Same, but for 1 column back

%%
act_thresh     = 0; % Defines a neuron as active or non-active. In Loebel et al. 2007, this was chosen as 0.
%% 
cancel_syns = 0;
pr_syn_canc = 0.5; % Probability of synapse cancelation
overlap     = 0; % Whether or not there may be neurons responsive to both F2 and F1
zipper      = 1; % This means that canceled synapses in the specified columns corresponding to stimulus frequencies follow a zipper pattern
first_canc  = 1; % In case of zipper pattern, the first synapse canceled in F1
        
%% Sensory Input
lambda_c = 2.5/(log(2)^0.5);

lambda_s_right = lambda_c + (A_max > alpha)*(A_max - alpha)/delta_right; % Determines localiztion of input effect over the cortical sheet
lambda_s_left  = lambda_c + (A_max > alpha)*(A_max - alpha)/delta_left;

if ring_net
    Distances = toeplitz([0:floor(M/2) wrev(1:(floor(M-1)/2))]); % Distance of each frequency from all the others.
else
    Distances = toeplitz(0:(M-1));
end

Right_Dec_Args = Distances./lambda_s_right;
Left_Dec_Args  = Distances./lambda_s_left;

curve_type = 'linear'; % Type of the tuning curve; see cases for options and details.

switch(curve_type)
    case 'exponential' % this one is the original type, used in Loebel et al. 2007. 
        h_outline = (triu(exp(-Right_Dec_Args)) + tril(exp(-Left_Dec_Args),-1)); % Magnitudes of input received by each column for each presented frequency
    
    case 'biased_exponential' % this is similar, but would be sharper and have inhibition from all farther off frequencies
        curve_bias = 1/3;
        h_outline  = (1 + curve_bias)*(triu(exp(-Right_Dec_Args)) + tril(exp(-Left_Dec_Args),-1)) - curve_bias; % Ensures a magnitude of A at the tuning curve's peak
    
    case 'linear' % This type is triangular, following Eli's request that some columns do not receive input at all.
        h_outline = (triu(1-Right_Dec_Args) + tril(1-Left_Dec_Args,-1));
        h_outline(h_outline < 0) = 0;
        
    case 'mexican hat' % My own addition; has support in STRFs such as in deCharms et al. 1998; Loebel et al. 2007 state that lateral supression arises in the model due to synaptic depression 
        mex_dec = 1;
        mex_fre = 1;
        h_outline = (triu(cos(mex_fre*Right_Dec_Args).*exp(-mex_dec*Right_Dec_Args)) + tril(cos(mex_fre*Left_Dec_Args).*exp(-mex_dec*Left_Dec_Args),-1));

    case 'Gaussian' 
        h_outline = (triu(exp(-Right_Dec_Args.^2)) + tril(exp(-Left_Dec_Args.^2),-1)); % Magnitudes of input received by each column for each presented frequency
    
end

h_outline = h_outline(1:freq_stretch:end,:); % This gives the columns best frequencies that are spread evenly across the range.
h_outline(P+1:end,:) = []; % Cropping the rows of h_outline to the no. of columns, just in case.

h_outline = reshape(h_outline,[P 1 M]); % *) Aligning the tuning curves to the neurons.
h = repmat(h_outline,[1 NE 1]);

inp2inh = 0; % This states whether there is sensory input to inhibitory neurons

%% Plotting the Tuning Curve 
if plot_TC %&& all_plots
    tc_column = Sel_Columns(3); %floor(P/2);
    
    figure(1)
        hold on
        plot(reshape(h_outline(tc_column,1,:),1,M),'-b.')

        title(['Typical Tuning Curve (Column ' num2str(tc_column) ')'])
        xlabel('Frequency')
        ylabel('Input as Fraction of the Input at the Best Frequency')
        xlim([0 22])
end

%% Mixing the Responses of Neurons in Each Column
if mix_tcs
   h_unmixed = h; % Original h is saved for FRA finding in unmixed case.
   shif = freq_stretch; % The basic shift in BF.
   placer1 = floor(NE/part_1f);
   placer2 = placer1 + floor(NE/part_1b);
   placer3 = placer2 + floor(NE/part_2f);
   placer4 = placer3 + floor(NE/part_2b);
   
   for Q = 1:P % Doing an independent scramble for each column
        Mixed   = randperm(NE);
        
        h(Q,Mixed(1:placer1),:)              = h(Q,Mixed(1:placer1),[((M-shif+1):M) 1:(M-shif)]); % forward means to the right, i.e. increase in 3rd dimension
        h(Q,Mixed(1:placer1),((M-shif+1):M)) = ring_net*h(Q,Mixed(1:placer1),((M-shif+1):M));
        
        h(Q,Mixed((placer1+1):placer2),:)      = h(Q,Mixed((placer1+1):placer2),[(shif+1):M (1:shif)]);
        h(Q,Mixed((placer1+1):placer2),1:shif) = ring_net*h(Q,Mixed((placer1+1):placer2),1:shif);
        
        h(Q,Mixed((placer2+1):placer3),:)              = h(Q,Mixed((placer2+1):placer3),[((M-2*shif+1):M) 1:(M-2*shif)]);
        h(Q,Mixed((placer2+1):placer3),(M-2*shif+1):M) = ring_net*h(Q,Mixed((placer2+1):placer3),(M-2*shif+1):M);
        
        h(Q,Mixed((placer3+1):placer4),:)        = h(Q,Mixed((placer3+1):placer4),[(2*shif+1):M (1:2*shif)]);
        h(Q,Mixed((placer3+1):placer4),1:2*shif) = ring_net*h(Q,Mixed((placer3+1):placer4),1:2*shif);
    end
end

%% Plotting the Best Frequencies to Show Mixing
[h_Peaks Best_Freqs] = max(h,[],3);

if plot_mixing && all_plots
    figure(2)
        imagesc(Best_Freqs)
        title('Best Frequencies of All Neurons')
        xlabel('No. of Neuron in Column')
        ylabel('Column')
        colorbar
end

%% Initializing Variables:
%t_prot    = n_stim*(duration + ISI) + 2*post_stim;

e_step = (bg_high_E - bg_low_E)/(NE - 1); % *) How about input from a non-uniform (e.g. Gaussian) distribution?
i_step = (bg_high_I - bg_low_I)/(NI - 1);
Inp_E  = ones(P,1)*(bg_low_E:e_step:bg_high_E); % These are the inputs to all the excitatory neurons, hence the P columns and NE rows
Inp_I  = ones(P,1)*(bg_low_I:i_step:bg_high_I); % These are the inputs to all the inhibitory neurons

tmax         = t_eq;% + t_prot*stim; % Maximum time that the simulation will reach (in seconds) 
num_steps_eq = floor(tmax/dt); % Total number of steps in the simulation

E   = zeros(P,NE); % Acitivity of all excitatory neurons (in Hz)
I   = zeros(P,NI); % Activity of all inhibitory neurons (in Hz)
x   = zeros(P,NE); % The fractions of resources available for synaptic transmission in all excitatory neurons
y   = zeros(P,NI); % Same as x, only for inhibitory neurons
z   = 0;% Only one synapse is needed here; zeros(P,NE,M_aug); % Same as x, only for the thalamo-cortical synapses

% *) Why not initialize x, y at 1? - I checked; it just creates an
% initial peak which goes to ~0 immediately for the rest of t_eq.
% z should be started at 1, or else there is no repsonse to the first stimuli.

Gain_E = zeros(P,NE); % The change in E in a certain time-step
Gain_I = zeros(P,NI); % The change in I in a certain time-step

E_act  = zeros(P,NE,floor(t_eq/dt)); % Activity of all excitatory neurons during the time allowed for reaching equilibrium, for identification of the active neurons.
E_mean = zeros(P,num_steps_eq); % Mean excitatory neuron activities in all columns, in all time-steps (in Hz)
%I_mean = E_mean; % Mean inhibitory neuron activities in all columns, in all time-steps (in Hz)

EUx     = zeros(P,1);
EUx_1   = EUx;
EUx_2   = EUx;
IUy     = EUx;
E_sum   = EUx;
E_sum_1 = EUx;
E_sum_2 = EUx;

%% The Dynamic Loop: Initial run to obtain steady-state and find initial values
% First, the network is allowed to reach a steady-state with no sensory input. 
% Tha activity of all neurons is assessed to determine which will be presnted with input.

for i = 1:floor(t_eq/dt)
    % Pre-calculation for the inter-column exc2exc gain calculations:
    EUx    = diag(E*(U.*x)'); % Intra-column input from excitatory synapses, taking synaptic depression into account
    EUx_1 = [EUx(2:P) ; ring_net*EUx(1)] + [ring_net*EUx(P); EUx(1:P-1)]; % Excitatory input from neighboring column
    EUx_2 = [EUx(3:P) ; ring_net*EUx(1:2)] + [ring_net*EUx(P-1:P); EUx(1:P-2)]; % Excitatory input from the neighboring column once removed

    % Pre-calculation for the intra-column inhibitory gain: 
    IUy = diag(I*(U.*y)'); % Intra-column input from inhibitory synapses, taking synaptic depression into account
    
    % Pre-calculation for the inter-column exc2inh gain:
    E_sum   = sum(E,2); % Total intra-column input in each column
    E_sum_1 = [E_sum(2:P); ring_net*E_sum(1)] + [ring_net*E_sum(P); E_sum(1:P-1)]; % Total input to each column from its nearest neighbors
    E_sum_2 = [E_sum(3:P); ring_net*E_sum(1:2)] + [ring_net*E_sum(P-1:P); E_sum(1:P-2)]; % Same, for 2nd-nearest neighbors.
     
    % The gain each cell receives:
    Gain_E = bsxfun(@plus,Inp_E,(Jee.*EUx + Jei.*IUy + Jee_1.*EUx_1 + Jee_2.*EUx_2));
    Gain_I = bsxfun(@plus,Inp_I,(Jie.*E_sum + Jii.*sum(I,2) + Jie_1.*E_sum_1 + Jie_2.*E_sum_2));
    %Gain_E = (Jee.*EUx + Jei.*IUy + Jee_1.*EUx_1 + Jee_2.*EUx_2)*ones(1,NE) + Inp_E;
    %Gain_I = (Jie.*E_sum + Jii.*sum(I,2) + Jie_1.*E_sum_1 + Jie_2.*E_sum_2)*ones(1,NI) + Inp_I;
    
    % Implementing the non-linearity:
    Gain_E(Gain_E < 0)   = 0;
    Gain_E(Gain_E > 300) = 300;
    Gain_I(Gain_I < 0)   = 0;
    Gain_I(Gain_I > 300) = 300;
    
    % The variables' dynamics:
    E = E + (dt/tau_E)*(-E + Gain_E.*(1 - tau_ref_E*E));
    x = x + dt*((1 - x)./tau_rec - U.*E.*x);         
    I = I + (dt/tau_I)*(-I + Gain_I.*(1 - tau_ref_I*I));
    y = y + dt*((1 - y)./tau_rec - U.*I.*y); 
    z = z + dt*((1 - z)./tau_rec_s);

    E_act(:,:,i)    = E; % Tracking the activity of all neurons
    E_mean(:,i)   = mean(E,2); % Tracking mean excitatory activity
    %I_mean(:,i)   = mean(I,2); % Tracking mean inhibitory activity
end

% Recording equilibrium conditions:
E_eq = E;
I_eq = I;
x_eq = x; 
y_eq = y;
z_eq = z;

Gain_E_eq = Gain_E;
Gain_I_eq = Gain_I;

%% Finding Active Neurons for Targeting the Thalamo-Cortical Connections
% In Loebel et al. 2007, thalamo-cortical connections were targeted only at
% neurons with spontaneous activity above 0, as this reproduces the
% experimentally pbserved phenomena of hypersensitive locking suppression.
Neuron_mean    = mean(E_act,3);

Active_Neurons = Neuron_mean > act_thresh;

if all_plots
    figure(4)
        subplot(2,1,1)                
            imagesc(Neuron_mean)

            title('Mean Activity of All Neurons')                
            xlabel('No. of Neuron in Column')                
            ylabel('Column')                
            colorbar

        subplot(2,1,2)                
            imagesc(Active_Neurons)

            title(['Active Neurons (Threshold = ' num2str(act_thresh) ')'])                
            xlabel('No. of Neuron in Column')                
            ylabel('Column')
            colorbar
end

%% Canceling Some of the Synapses 
% This is supposed to make the cells in each column different in their inputs.

if cancel_syns
    mid_column  = floor(mean([mean(F2) mean(F1)]));
    Spec_Cols   = [mid_column];
    if overlap
        h(rand(size(h)) > (pr_syn_canc)) = 0;
        
    elseif zipper
        Canc_Synapses = zeros(1,NE);
        Canc_Synapses(first_canc:2:end) = 1;
        h(Spec_Cols,~~Canc_Synapses,F1)  = 0;
        h(Spec_Cols,~Canc_Synapses,F2) = 0;
        
    else % Basically, this is valid only if the probability of cancelation is 0.5.
        num_cancel    = floor(NE*pr_syn_canc);
        Canc_Synapses = [ones(1,num_cancel) zeros(1,NE - num_cancel)];
        Canc_Synapses = Canc_Synapses(randperm(NE));
        h(Spec_Cols,~Canc_Synapses,F1)  = 0;
        h(Spec_Cols,~~Canc_Synapses,F2) = 0;
    end    
    
    % Note: The cancelation of some of the synapses was introduced while
    % trying to get SSA in a single column not connected to other columns.
    % It was found that the identity of canceled synapses is important,
    % because of the different background input to different cells.
    % However, the problem lied elsewhere and this feature was no longer
    % used. If needed, I suggest using it while allowing overlap.
    
    % Plotting the active and canceled synapese:
    if all_plots
        figure(3)
        for spec_col = 1:length(Spec_Cols)
            h_rearr = zeros(M_aug,NE);
            for syn_row = 1:M
                h_rearr(syn_row,:) = h(Spec_Cols(spec_col),:,syn_row);
            end

            subplot(length(Spec_Cols),1,spec_col)    
                imagesc(h_rearr)

                title(['Synapses of Column ' num2str(Spec_Cols(spec_col))])
                xlabel('Neuron')
                ylabel('Frequency')
        end
    end
end

%% Targeting of Active Neurons
h(repmat(~Active_Neurons,[1 1 M])) = 0; % Non-active neurons do not get sensory input. This is achieved by canceling their tuning curves.

%% Saving
E_act = [];% just to save memory

filename = [ '/TYLT_Prot_' term '_Network_' num2str(net_num) '_Par.mat'];

save(filename) 
save([ '/TYLT_filename_' term '.mat'],'filename')

disp('Initialization done')
