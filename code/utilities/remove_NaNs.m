function [vec_RDM1, vec_RDM2] = remove_NaNs(RDMoi1, RDMoi2)

h = height(RDMoi1.RDM);
vec_RDM1 = RDMoi1.RDM(tril(true(h),-1)); % Extract lower triangular part, without the diagonal
vec_RDM2 = RDMoi2.RDM(tril(true(h),-1));

% check for missing values and remove from vector
nans = isnan(vec_RDM1);
if sum(nans) > 0
    vec_RDM1 = vec_RDM1(~nans);
    vec_RDM2 = vec_RDM2(~nans);

    if numel(vec_RDM1) ~= numel(vec_RDM2)
        error('Size of matrices not equal due to missing values')
    end
end

end