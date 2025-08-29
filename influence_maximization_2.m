function S = select_seed_nodes(G, k, M, p)
    S = [];
    V = 1:numnodes(G);  % Assuming G is a graph object
    ddv = degree(G); 

    for i = 1:k
        [~, idx] = max(ddv(setdiff(V, S)));
        u = setdiff(V, S)(idx);
        S = [S, u];

        neighbors_u = neighbors(G, u);
        for v = neighbors_u'
            if ~ismember(v, S)
                ddv(v) = degree(G, v) - 1 - (degree(G, u) - 1) * p;
            end
        end
    end
end

G = random_graph(10000, 0.5);  % Create a random graph

S = select_seed_nodes(G, 15, 0, 0.5);

Gk = karate_club_graph();

s = select_seed_nodes(Gk, 10, 0, 0.5);

function [A, mean_spread] = IC(g, S, p, mc)
    spread = zeros(mc, 1);
    activated_nodes = [];

    for i = 1:mc
        new_active = S;
        A = S;

        while ~isempty(new_active)
            new_ones = [];
            for node = new_active
                rng(i);  % Set random seed
                neighbors_node = neighbors(g, node);
                success = rand(length(neighbors_node), 1) < p;
                new_ones = [new_ones; neighbors_node(success)];
            end
            new_active = setdiff(new_ones, A);
            A = unique([A, new_active]);
        end

        spread(i) = length(A);
        activated_nodes = A;
    end

    mean_spread = mean(spread);
end

G = random_graph(100, 0.5);

[spread, mean_spread] = IC(Gk, s, 0.5, 1000);

node_colors = arrayfun(@(node) ifelse(ismember(node, spread), 'red', 'blue'), Gk.Nodes, 'UniformOutput', false);

pos = layout(Gk);  % Layout the graph

plot(Gk, 'NodeColor', node_colors, 'NodeSize', 500);

num_nodes = numnodes(Gk);

G = random_graph(10000, 0.2);

num_edges = numedges(G);
num_nodes = numnodes(G);

plot(G, 'NodeSize', 10);

s = select_seed_nodes(G, 1000, 0, 0.5);

[spread, mean_spread] = IC(G, s, 0.5, 100);

node_colors = arrayfun(@(node) ifelse(ismember(node, spread), 'red', 'blue'), G.Nodes, 'UniformOutput', false);

pos = layout(G);  % Layout the graph

plot(G, 'NodeColor', node_colors, 'NodeSize', 10);

mean_spread;