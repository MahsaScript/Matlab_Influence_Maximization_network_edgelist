% Read the data from the text file
file_path = 'higgs-activity_time.txt';  % Replace with the actual path to your file
columns = {'UserA', 'UserB', 'Timestamp', 'Interaction'};

df = readtable(file_path, 'Delimiter', ' ', 'ReadVariableNames', false);
df.Properties.VariableNames = columns;

% Display the first few rows of the table
disp(head(df));

% Create a directed graph
G = digraph();

% Add nodes from 'UserA' and 'UserB' columns
G = addnode(G, unique([df.UserA; df.UserB]));

% Add edges based on interactions
for index = 1:height(df)
    user_a = df.UserA(index);
    user_b = df.UserB(index);
    interaction = df.Interaction(index);
    
    % Depending on the interaction type, add edges
    if strcmp(interaction, 'RT') || strcmp(interaction, 'MT') || strcmp(interaction, 'RE')
        G = addedge(G, user_a, user_b, interaction);
    end
end

num_nodes = numnodes(G);
num_edges = numedges(G);

% Propagation graph creation

% Sample 500 nodes
sampled_nodes = unique([df.UserA; df.UserB]);
sampled_nodes = sampled_nodes(randperm(length(sampled_nodes), min(500, length(sampled_nodes))));

% Create a subgraph with sampled nodes
subgraph = subgraph(G, sampled_nodes);

% Draw the subgraph with different edge colors based on interaction types
edge_colors = containers.Map({'RT', 'MT', 'RE'}, {'r', 'g', 'b'});

edge_colors_list = cell(1, numedges(subgraph));
for i = 1:numedges(subgraph)
    interaction = subgraph.Edges.Weight{i};
    edge_colors_list{i} = edge_colors.isKey(interaction) * edge_colors(interaction) + ~edge_colors.isKey(interaction) * 'k';
end

% Draw the subgraph
figure;
plot(subgraph, 'EdgeColor', edge_colors_list, 'NodeColor', 'k', 'MarkerSize', 10);
legend(edge_colors.keys, 'Location', 'best');
title('Subgraph Visualization');

% Linear Threshold model

% Linear Threshold model function : takes a graph and seed nodes and returns activated nodes
function activated_nodes = linear_threshold_model(graph, seed_nodes, max_iterations)
    activated_nodes = unique(seed_nodes);
    new_nodes_activated = true;
    iteration = 0;

    while new_nodes_activated && iteration < max_iterations
        new_nodes_activated = false;
        for node = 1:numnodes(graph)
            if ~ismember(graph.Nodes.Name{node}, activated_nodes)
                neighbors = neighbors(graph, graph.Nodes.Name{node}, 'all');  % Incoming neighbors
                activated_neighbors = intersect(neighbors, activated_nodes);
                weights_sum = sum(graph.Edges.Weight(ismember(graph.Edges.EndNodes(:, 1), activated_neighbors) & ismember(graph.Edges.EndNodes(:, 2), graph.Nodes.Name{node})));

                if weights_sum >= graph.Nodes.Threshold(node)
                    activated_nodes = [activated_nodes; graph.Nodes.Name{node}];
                    new_nodes_activated = true;
                end
            end
        end
        iteration = iteration + 1;
    end
end

% Linear threshold model with random weights and threshold.

% Add attributes to the graph (replace 'weight' and 'threshold' with your actual attribute names)
for i = 1:numedges(subgraph)
    subgraph.Edges.Weight(i) = rand();  % Assign random weights
end

for i = 1:numnodes(subgraph)
    subgraph.Nodes.Threshold(i) = rand();  % Assign random thresholds
end

% Seed nodes (replace with your actual seed nodes)
seeds = sampled_nodes(randperm(length(sampled_nodes), 5));

% Run Linear Threshold Model
activated_nodes = linear_threshold_model(subgraph, seeds, 1000);

% Visualize the result
figure;
node_colors = cell(1, numnodes(subgraph));
for i = 1:numnodes(subgraph)
    if ismember(subgraph.Nodes.Name{i}, activated_nodes)
        node_colors{i} = 'r';
    else
        node_colors{i} = 'b';
    end
end
plot(subgraph, 'NodeColor', node_colors, 'MarkerSize', 10);
title('Activated Nodes Visualization');

% Create a table for the nodes in subgraph
subgraph_nodes = subgraph.Nodes.Name;
subgraph_nodes_table = table(subgraph_nodes, 'VariableNames', {'User'});

% Merge subgraph_nodes_table with df to get the corresponding rows from the original table
result_table = innerjoin(subgraph_nodes_table, df, 'Keys', 'UserA');

% Display the resulting table
disp(result_table);

% Dictionary to store propagation graphs for each action type
propagation_graphs = containers.Map;

% Convert 'Timestamp' column to datetime format
result_df.TimestampDate = datetime(result_df.Timestamp, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Iterate over unique action types
unique_action_types = unique(result_df.Interaction);
for action_type = unique_action_types'
    % Create a directed graph for propagation for each action type
    propagation_graph = digraph;

    % Add nodes and edges based on propagation conditions
    action_rows = result_df(strcmp(result_df.Interaction, action_type), :);
    for index = 1:height(action_rows)
        vi = action_rows.UserA(index);
        vj = action_rows.UserB(index);
        ti = action_rows.TimestampDate(index);

        if ismember(vi, subgraph2.Edges.EndNodes(:, 1)) && ismember(vj, subgraph2.Edges.EndNodes(:, 2))
            % Check if there's a previous action (vi, a, ti) and (vj, a, tj)
            prev_actions_vj = result_df(strcmp(result_df.UserB, vj) & strcmp(result_df.Interaction, action_type) & result_df.TimestampDate > ti, :);
            if ~isempty(prev_actions_vj)
                % Check the condition (iii) T(vi, vj) <= ti
                delta_t = ti - min(prev_actions_vj.TimestampDate);
                propagation_graph = addedge(propagation_graph, vi, vj, delta_t);
            end
        end
    end

    % Store the propagation graph for the current action type
    propagation_graphs(char(action_type)) = propagation_graph;
    fprintf('%d %d\n', numnodes(propagation_graph), numedges(propagation_graph));
end

n_nodes_RT = numnodes(propagation_graphs('RT'));
n_edges_RT = numedges(propagation_graphs('RT'));

% Display propagation graphs
disp(propagation_graphs);

dx = result_df(1:5000, :);

% Create a directed graph
G = digraph;

% Add nodes from 'UserA' and 'UserB' columns
G = addnode(G, unique([dx.UserA; dx.UserB]));

% Add edges based on interactions
for index = 1:height(dx)
    user_a = dx.UserA(index);
    user_b = dx.UserB(index);
    interaction = dx.Interaction(index);

    % Depending on the interaction type, add edges
    if strcmp(interaction, 'RT') || strcmp(interaction, 'MT') || strcmp(interaction, 'RE')
        G = addedge(G, user_a, user_b, interaction);
    end
end

n_nodes_G = numnodes(G);
n_edges_G = numedges(G);

% Probability Calculation
nodes_list = subgraph2.Nodes.Name;

% Initialize a matrix to store probabilities
probability_matrix = zeros(length(nodes_list), length(nodes_list));

% Iterate over propagation graphs
for action = keys(propagation_graphs)
    propagation_graph = propagation_graphs(action{1});
    for i = 1:length(nodes_list)
        for j = 1:length(nodes_list)
            if i ~= j  % Exclude self-loops
                % Calculate the number of actions influenced from node v to node u
                if ismember(nodes_list{i}, propagation_graph.Nodes.Name)
                    if ismember(nodes_list{j}, successors(propagation_graph, nodes_list{i}))
                        Av2u = 1;
                        % Calculate the total number of actions performed by node v
                        Av = outdegree(propagation_graph, nodes_list{i});
                        % Add up the probabilities from each graph
                        probability_matrix(i, j) = probability_matrix(i, j) + (Av2u / Av);
                    end
                end
            end
        end
    end
end

% Display the combined probability matrix
disp('Combined Probability Matrix:');
disp(probability_matrix);

unique_probabilities = unique(probability_matrix);

% "Linear Threshold Model with propagation probability"

% Required MATLAB dependencies:
% Uses built-in functions for random numbers, graph operations, and plotting.
% No additional toolbox imports are necessary beyond the default MATLAB graph and plotting functions.

function main
    %%% Define or load probability_matrix and subgraph2 as needed %%%
    % For demonstration purposes, create a sample probability_matrix and subgraph2.
    % NOTE: Replace these with your actual data.
    numNodes = 50; % example number of nodes
    probability_matrix = rand(numNodes); % sample probability matrix with values between 0 and 1
    % Create a directed graph for subgraph2 with numNodes nodes and random edges.
    % Here, we use an undirected graph since MATLAB's graph creates undirected graphs 
    % by default. To mimic directed behavior, one could use digraph.
    s = randi(numNodes, numNodes, 1); 
    t = randi(numNodes, numNodes, 1);
    subgraph2 = graph(s, t);
    
    % Initialize global thresholds vector to store thresholds for each node.
    global thresholds
    thresholds = zeros(numnodes(subgraph2), 1);
    
    % Linear threshold model with probability matrix
    % Assign random thresholds to each node.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for node in subgraph2.nodes():
    %     subgraph2.nodes[node]['threshold'] = random.random()  % Assign random thresholds
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In MATLAB, we assign thresholds to the Nodes table and also update the global thresholds.
    numSubNodes = numnodes(subgraph2);
    randomThresholds = rand(numSubNodes, 1);  % Generate random thresholds for each node
    subgraph2.Nodes.threshold = randomThresholds;
    thresholds = randomThresholds;  % Update global thresholds so that the model function can access them

    % Seed nodes (replace with your actual seed nodes)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % seeds = [random.choice(list(subgraph2.nodes)) for _ in range(100)]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In MATLAB, randomly choose 100 seed nodes from the available nodes.
    seeds = randi(numSubNodes, 100, 1);

    % Run Linear Threshold Model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % activated_nodes = linear_threshold_model(probability_matrix, seeds)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    activated_nodes = linear_threshold_model(probability_matrix, seeds);

    % Visualize the result
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pos = nx.spring_layout(subgraph2)
    % nx.draw(subgraph2, pos, with_labels=False, font_weight='bold', node_size=10, ...
    %         node_color=['red' if node in activated_nodes else 'blue' for node in subgraph2.nodes()])
    % plt.show()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In MATLAB, we compute a force-directed layout using the 'force' option,
    % then plot the graph and color the nodes: red if activated, blue otherwise.
    figure;
    p = plot(subgraph2, 'Layout', 'force');
    hold on;
    % Determine node colors based on activation status.
    nodeColors = repmat([0 0 1], numSubNodes, 1);   % blue for non-activated nodes
    for i = 1:length(activated_nodes)
        % Since MATLAB nodes are 1-indexed and activated_nodes contains indices, 
        % we set the activated node color to red.
        nodeColors(activated_nodes(i), :) = [1 0 0];  % red for activated nodes
    end
    % Overplot the nodes with the specified colors.
    scatter(p.XData, p.YData, 10, nodeColors, 'filled');
    hold off;
    title('LTM with probability matrix');
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Threshold Model Function with propagation probability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function activated_nodes = linear_threshold_model_with_probabilities(probability_matrix, seed_nodes, max_iterations)
    % Import global thresholds to mimic node['threshold'] from the original Python code.
    global thresholds

    % num_nodes = len(probability_matrix)
    num_nodes = size(probability_matrix, 1);
    
    % activated_nodes = set(seed_nodes)
    activated_nodes = unique(seed_nodes);
    
    % new_nodes_activated = True
    new_nodes_activated = true;
    
    % iteration = 0
    iteration = 0;
    
    % while new_nodes_activated and iteration < max_iterations:
    while new_nodes_activated && (iteration < max_iterations)
        % new_nodes_activated = False
        new_nodes_activated = false;
        
        % for node in range(num_nodes):
        %     if node not in activated_nodes:
        for node = 1:num_nodes
            if ~ismember(node, activated_nodes)
                % Get the incoming neighbors of the node
                % neighbors = [i for i in range(num_nodes) if probability_matrix[i, node] > 0]
                neighbors = [];
                for i = 1:num_nodes
                    if probability_matrix(i, node) > 0
                        neighbors(end+1) = i; %#ok<AGROW>
                    end
                end

                % Calculate the sum of probabilities from activated neighbors
                % weights_sum = sum(probability_matrix[neighbor, node] for neighbor in neighbors if neighbor in activated_nodes)
                weights_sum = 0;
                for j = 1:length(neighbors)
                    neighbor = neighbors(j);
                    if ismember(neighbor, activated_nodes)
                        weights_sum = weights_sum + probability_matrix(neighbor, node);
                    end
                end

                % Compare the sum with the threshold
                % if weights_sum >= node['threshold']:  
                %     activated_nodes.add(node)
                %     new_nodes_activated = True
                % In MATLAB, we mimic node['threshold'] with global thresholds.
                if weights_sum >= thresholds(node)
                    activated_nodes = unique([activated_nodes; node]);  %#ok<AGROW>
                    new_nodes_activated = true;
                end

            end
        end
        
        % iteration += 1
        iteration = iteration + 1;
    end
    
    % return activated_nodes
    % In MATLAB, the output variable activated_nodes is returned.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wrapper function to preserve original function call name in Python code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function activated_nodes = linear_threshold_model(probability_matrix, seed_nodes)
    % Calls the linear_threshold_model_with_probabilities with default max_iterations = 1000
    activated_nodes = linear_threshold_model_with_probabilities(probability_matrix, seed_nodes, 1000);
end

% Execute main function
main
