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
