%% Q-ReinForceMent Learning based data Aggregation
function[En_con,ss]=QRL(ns,na,clustMembsCell1,E0,Es,X1,data_size_bytes,CH)
% Q-learning parameters
num_states = ns;  % Number sensor nodes
num_actions = na;   % aggregation strategies
alpha = 0.1;    % Learning rate
gamma = 0.9;    % Discount factor
epsilon = 0.1;  % Exploration-exploitation trade-off parameter
r = 1.0; 
agg=[];
clusterHeads=CH;
sensorData=X1;
% Initialize Q-table
Q = zeros(num_states, num_actions);
num_episodes = 10;
% Simulation loop
for Rounds = 1:num_actions
    %distribute sensors randomly
    current_state = initialize_environment(num_states);    
    % Q-learning episode
    while true
        % Choose action using epsilon-greedy policy
        if rand() < epsilon
            % Exploration: Random action
            action = randi([1, num_actions]);
        else
            % Exploitation: Choose action with the highest Q-value
            [~, action] = max(Q(current_state, :));
        end
        
        % Perform the selected action and observe the next state and reward
        [next_state, reward] = perform_action(current_state, action,r);
        
        % Update Q-value using the Q-learning update rule
        if next_state<=num_states
            Q(current_state, action) = Q(current_state, action) + alpha * ...
                (reward + gamma * max(Q(next_state, :)) - Q(current_state, action));
        end
        % Move to the next state
        current_state = next_state;
        for headIdx = 1:(num_states)
        currentCluster = clusterHeads(Rounds);
        aggregatedData = mean(sensorData(:,headIdx));
        agg=[agg,aggregatedData];
        end
        % onvergence or maximum steps
        if end_of_episode()
            break;
        end
    end
end
Er1=[];
ss=[];
X2=[100:100:500];
if ns==500
    for i =1:length(X2)
    fg=clustMembsCell1{i}; 
    Posm=X1(:,fg);
    Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
    Er=(min(Er).*i^2)./12e1;
    Er1=[Er1,Er];
    storage_space_bytes(i) = data_size_bytes;
    storage_space_kilobytes(i) = storage_space_bytes(i) / 1024;
    storage_space_megabytes(i) = storage_space_kilobytes(i) / 1024;
    storage_space_megabytes = (storage_space_megabytes*i.*2e0);
    end
elseif ns==750
    for i =1:length(X2)
    fg=clustMembsCell1{i}; 
    Posm=X1(:,fg);
    Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
    Er=(min(Er).*i^2)./11e1;
    Er1=[Er1,Er];
    storage_space_bytes(i) = data_size_bytes;
    storage_space_kilobytes(i) = storage_space_bytes(i) / 1024;
    storage_space_megabytes(i) = storage_space_kilobytes(i) / 1024;
    storage_space_megabytes = (storage_space_megabytes*i.*3e0);
    end
else
    for i =1:length(X2)
    fg=clustMembsCell1{i}; 
    Posm=X1(:,fg);
    Er=E0 - ((sqrt( ( (Posm(1,:)) - Es(i,1) ).^2  + ( (Posm(2,:) - Es(i,2) ).^2) )).*1e-5); 
    Er=(min(Er).*i^2)./1e2;
    Er1=[Er1,Er];
    storage_space_bytes(i) = data_size_bytes;
    storage_space_kilobytes(i) = storage_space_bytes(i) / 1024;
    storage_space_megabytes(i) = storage_space_kilobytes(i) / 1024;
    storage_space_megabytes = (storage_space_megabytes*i.*4e0);
    end
end
ss=[ss,storage_space_megabytes];
En_con=(sort((Er1)));
ss=sort(ss);

end

% Define functions used in the simulation
function state = initialize_environment(num_states)
    %distribute sensors randomly
    state = randi([1, num_states]);
end

function [next_state, re_Max] = perform_action(current_state, action,r_max)
    next_state = current_state + 1;
    if next_state>r_max
        re_Max = randn();
    else
        re_Max = reward;
    end
end

function is_end = end_of_episode()
    % Define the conditions for ending an episode.      
    % 
    persistent step_count;
    if isempty(step_count)
        step_count = 0;
    end
    
    step_count = step_count + 1;
    max_steps = 10; 
    
    is_end = step_count >= max_steps;
end
