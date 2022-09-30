% Example : Input is 4-dim matrix
% output = [in(:, :, 1, 1)        0        0          0        0         0          ]
%          [       0       in(:, :, 2, 1)  0          0        0         0          ]
%          [       0              0       ...         0        0         0          ]
%          [       0              0        0   in(:, :, 1, 2)  0         0          ]
%          [       0              0        0          0       ...        0          ]
%          [       0              0        0          0        0  in(:, :, end, end)]
function output = blockdiag(mat)
mat_size = size(mat);
mat_buf = repmat({zeros(mat_size(1), mat_size(2))}, prod(mat_size(3:end)));
mat_buf((1:prod(mat_size(3:end)):numel(mat_buf)) + (0:prod(mat_size(3:end))-1))...
  = num2cell(mat, [1 2]);

output = cell2mat(mat_buf);
end