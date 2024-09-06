function result = classify_nodes(A)

    row_norms = sum(abs(A), 2);

    col_norms = sum(abs(A), 1);

    N = size(A, 1);

    [~, row_indices] = sort(row_norms);
    row_ranks = zeros(N, 1);
    for i = 1:N
        row_ranks(row_indices(i)) = i / N;
    end

    [~, col_indices] = sort(col_norms);
    col_ranks = zeros(1, N);
    for i = 1:N
        col_ranks(col_indices(i)) = i / N;
    end

    result = cell(N, 1);

    for node = 1:N
        row_rank = row_ranks(node);
        col_rank = col_ranks(node);

        sink_value = sqrt(2) - norm([(row_rank - col_rank) - 1, 1/N]);

        source_value = sqrt(2) - norm([(row_rank - col_rank) - 1/N, 1]);

        if source_value > sink_value
            result{node} = 'source';
        else
            result{node} = 'sink';
        end
    end
end

