import networkx.*;

function S = dhicm_dynamic_p(G, k)
    S = [];
    V = G.nodes();  % Assuming G is a NetworkX graph or a similar graph representation
    ddv = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    
    for v = V
        ddv(v) = G.in_degree(v);
    end

    for i = 1:k
        u = max(setdiff(V, S), 'ComparisonMethod', 'real', 'Key', @(v) ddv(v));
        S = [S, u];

        for v = G.predecessors(u)
            if ~ismember(v, S)
                if ismember(u, G.predecessors(v))
                    common_neighbors = intersect(G.predecessors(v), G.predecessors(u));
                    p = 0.01 + ((G.in_degree(u) + G.in_degree(v)) / G.number_of_nodes()) + (length(common_neighbors) / G.number_of_nodes());
                    p = p * G.get_edge_data(u, v).weight; % accounting for edge weight
                    ddv(v) = G.in_degree(v) - 1 - ((G.in_degree(u) - 1) * p);
                end
            end
        end
    end
end

function S = degree_centrality(G, k)
    S = [];
    V = G.nodes();  % Assuming G is a NetworkX graph or a similar graph representation
    ddv = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    
    for v = V
        ddv(v) = G.in_degree(v);
    end

    for i = 1:k
        u = max(setdiff(V, S), 'ComparisonMethod', 'real', 'Key', @(v) ddv(v));
        S = [S, u];
    end
end

function [mean_spread, A] = modified_ICM(graph_object, S, mc)
    spread = zeros(mc, 1);
    for i = 1:mc
        new_active = S;
        A = S;

        while ~isempty(new_active)
            targets = modified_propagate(graph_object, new_active);
            new_ones = [];

            rng(i); % Set random seed
            for j = 1:length(targets)
                curr_target = targets(j, 1);
                p = targets(j, 2);
                if rand < p
                    new_ones = [new_ones, curr_target];
                end
            end

            new_active = setdiff(new_ones, A);
            A = unique([A, new_active]);
        end

        spread(i) = length(A);
    end

    mean_spread = mean(spread);
end
x = [0  20 40 60 80 100];
y = [15 25 35 45 55 65];
plot(x,y,'-red','LineWidth',2);
hold on;
y = [28 38 48 58 68 78];
plot(x,y,'-blue','LineWidth',2);
title('Influence Speed vs. number of seeds for the PPO agent with and without the GCN');
ylabel('Influence Speed');
xlabel('Number of seed nodes');
legend('PPO with GCN','PPO without GCN', 'Location','southeast');
ylim([0 80]);
xlim([-5 110]);
function targets = modified_propagate(g, new_active)
    targets = [];
    for node = new_active
        for neighbor = g.predecessors(node)
            common_neighbors = intersect(g.predecessors(node), g.predecessors(neighbor));
            p = 0.01 + ((g.in_degree(node) + g.in_degree(neighbor)) / g.number_of_nodes()) + (length(common_neighbors) / g.number_of_nodes());
            p = p * g.get_edge_data(neighbor, node).weight; % accounting for edge weight
            targets = [targets; neighbor, p];
        end
    end
end

% creating graph
fileID = fopen('higgs-retweet_network.edgelist', 'r');
Lines = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);
edge_list = [];

for i = 1:length(Lines{1})
    u_v = str2double(strsplit(strtrim(Lines{1}{i})));
    edge_list = [edge_list; u_v(1), u_v(2), u_v(3)];
end

G = digraph(edge_list(:,1), edge_list(:,2), edge_list(:,3));  % using a list of edge tuples
total_edge_weight = sum(G.Edges.Weight);

total_edge_weight = 0;
for curr_node = G.Nodes
    for pred = G.predecessors(curr_node)
        total_edge_weight = total_edge_weight + G.Edges.Weight(G.findedge(pred, curr_node));
    end
end

%getting seed nodes using dhicm
start_time = tic;
s_dhicm = dhicm_dynamic_p(G, 1000);
s_degree = degree_centrality(G, 1000);
fprintf("Time getting seed set: %.2f\n", toc(start_time));

start_time = tic;
[mean_spread, A] = modified_ICM(G, s_dhicm, 10);
fprintf("Mean spread for dynamic p (DHICM): %.2f\n", mean_spread);
fprintf("Duration: %.2f\n", toc(start_time));

start_time = tic;
[mean_spread, A] = modified_ICM(G, s_degree, 10);
fprintf("Mean spread for dynamic p (Degree Centrality): %.2f\n", mean_spread);
fprintf("Duration: %.2f\n", toc(start_time));