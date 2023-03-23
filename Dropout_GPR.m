function [m_post, m_pred, c_pred] = dropout_gpr(Y, x, yt, sn, sp, p)

% training data size
n_subj = numel(x);

[~, ys] = size(Y);
if ys ~= n_subj
    Y = Y';
    yt = yt';
end

% prior covariance
K = Y' * (sp^2 * Y);

% dropout regularization parameter
alpha = 1 - p;

% approximate covariance
Kd = alpha * K + (1 - alpha) * (sp^2 * eye(n_subj));

% inverse of the approximate covariance
Kd_inv = inv(Kd);

% mean
m_post = (sp^2 / sn^2) * Kd_inv * Y * x;

% evaluation of the predictive distribution
% mean
m_pred = yt' * m_post;

% covariance
c_pred = diag(yt' * (sp^2 * yt) - yt' * Kd_inv * yt);
end