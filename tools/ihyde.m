function [result]=ihyde(parameter)

%
% -------------------------------------------------------------------------
% INPUT
% -------------------------------------------------------------------------
% Phi: matrix containing the regressors. Each row is a datapoint.
%
% y: column vector containing the output datapoints.
%
% s: maximal number of subsystems
%
% lambda: tuning parameter
%
% MAXITER:
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% idmodes: structure containing all information on the identified models.
%
%   idmodes.par{i}: is the PVs of the i-th mode.
%   idmodes.regions(i): is the region of the i-th mode. Regions are
%   subsets of idpar.Regressor_set.
%   idmodes.cov{i}; is the covariance of the PV of the i-th mode.
%   idmodes.Regressor_set: automatically computed polytope defining the
%   regressor set. This field is present only if idpar.Regressor_set is
%   empty or nonexistent.
%   idmodes.s: meaningful number of modes estimated during regression (it
%   might differ from idpar.s).
%   idmodes.regions_sim(i): is the region for simulating the i-th mode.
%   These regions are subsets of idpar.Regressor_set_sim and are defined
%   only if idpar.Regressor_set_sim is defined.
%   idmodes.stat_reattr: fraction of points reattribted to modes, if
%   reattribution has been performed.
%   idmodes.clust_valid: structure containing information about the
%   clustering results.
%   idmodes.pattern_rec_valid: structure containing information about the
%   pattern recognition results.
%   idmodes.Weight_primal: weights used in creating LDs (usually it is a
%   vector of 1's)
%   idmodes.adjacences: each row is a pair (i,j), i<j indicating that the
%   regions i and j are adjacent The constraint i<j avoids storing both
%   pairs (i,j) and (j,i).
% F: structure containing information about the mode datasets (the
% classified datapoints that are also inliers, i.e. not discarded during
% regression).
%
%   F.X{i}: matrix of regressors assigned to the i-th mode dataset.
%   F.y: cell array: F.y{i} vector of outputs assigned to the i-th mode
%   dataset.
%   F.pos{i}: indexes of the points assigned to the i-th mode dataset.
%   F.pos{i}(1)=5 means that Xid(5,:) and yid(5) are the first points
%   composing the i-th mode dataset.
%
% xi: structure containing information about the xi-points (they are either
% FVs or LPVs).
%
%	xi.points{i}: xi-point based on the i-th LDs.
%   xi.IR{i}: INVERSE of the covariance of the i-th xi-point.
%   xi.weights{i}: scalar confidence measure of the i-th feature vector.
%
% LDs: structure containing information about local datasets (LDs)
% and local models
%
%	LDs.X{i}:  matrix of regressors belonging to the i-th local dataset
%	(each row is a point).
%   LDs.y{i}: vector of outputs belonging to the i-th local dataset.
%   LDs.pos{i}: vector of indexes of datapoints in the i-th local
%   dataset, e.g. Xid(LDs.pos{i}(1),:) is the first regressor in the i-th
%   LD.
%   LDs.weights{i} weight associated to the i-th datapoint used in
%   weighted LS for computing mode PVs.
%   LDs.meanX{i} average of regressors in the i-th local dataset.
%   LDs.models{i} parameters of the i-th local model.
%   LDs.models_var{i} INVERSE variance of the i-th local model.
%   LDs.X_ivar{i} INVERSE of the variance of the the regressors in the
%   i-th LD.
%
% inliers: structure containing information on the inliers
%
%	inliers.pos(j): index of the j-th inlier in Xid and yid
%	(i.e.Xid(inliers.pos(j),:) and yid(inliers.pos(j)) are the j-th
%	inliers).
%   inliers.class(j): classification of the j-th inlier.

% Copyright is with the following author:
%
% (C) 2016 Ye Yuan,
Phi = parameter.Phi;
lambda = parameter.lambda;
epsilon = parameter.epsilon;
MAXITER = parameter.MAXITER;
y = parameter.y;
result = parameter;
deleted_idx = [];
for k = 1 : parameter.max_s

    T = null(Phi')';

    y1 = T*y;

 
    
%         ty = max(abs(y1));
%         y1 = y1/ty;
  

    Omega = (T*T')^(-1);

    % Normalization
%     for i = 1 : size(T,2)
%         t(i) = norm(T(:,i),2);
%         if t(i)==0
%             t(i)=1;
%         end
%         T(:,i) = T(:,i)/t(i);
%     end
%     
%     w =  tac_reconstruction(y1, T, lambda(1), MAXITER, Omega,1);
%     z = w(:,end)./t'*ty;
    solver_parameter.y = y1;
    solver_parameter.A = T;
    solver_parameter.lambda = lambda(1);
    solver_parameter.MAXITER = MAXITER;
    solver_parameter.Omega = Omega;
    solver_parameter.y_normalize = 1;
    solver_parameter.if_print = 1;
    solver_parameter.whole_index = 1;
    fprintf(2,'Finding %dth subsystem\n', k);
    z = solver( solver_parameter);
%     error_l2 = norm(y1 - T*theta,2)
%     error_l2_l1 = (lambda(1)*norm( theta, 1)+ 0.5*(TT*theta-y1)'*Omega*(TT*theta-y1))
    

    
    n =  size(y,1);
    

    delete_idx = [];


    result.theta(:,k) = abs(z);    %save theta in check for tunning lambda and epsilon
    for i = 1 : n
        if abs(z(i)) < epsilon(1)  %%compared with error,it is more convenient to find a better epsilon by choosing  theta
            delete_idx = [delete_idx i];
        end
    end
    
    
    %
    if(k > 1)
        idx_sys{k} = setdiff(delete_idx, deleted_idx'); %   only save the kth delete_idx
        
    else
        idx_sys{k} = delete_idx;
    end
    
    
    
    
    if size(idx_sys{k},2) == 0
        result.s = k;
        break;
    end
    [delete_idx,    alone]= find_idx( idx_sys{k},0);%Find the longest continuous sequence
    result.longgest_idx{k} = delete_idx;
    
    
    
    
    
    
    T = Phi(delete_idx,:);
    solver_parameter.y = y(delete_idx);
    solver_parameter.A = T;
    solver_parameter.lambda = lambda(2);
    solver_parameter.MAXITER = MAXITER;
    solver_parameter.y_normalize = 1;
    solver_parameter.Omega = 1;
    solver_parameter.if_print = 1;
    solver_parameter.whole_index = 1;
    fprintf(2, 'Identifing the %dth subsystem\n',k);
    sys(:,k) = solver( solver_parameter);
    
%     parameter.delete_state = parameter.state(delete_idx);
    
    result.sys(:,k) = sys(:,k);
    error = abs(y - Phi * sys(:,k));
    result.sys_allerror(:,k) = abs(parameter.y - parameter.Phi * sys(:,k));
    result.error(:,k) = error;
    delete_idx = [];
    for i = 1 : n
        if error(i) < epsilon(2)  %%compared with error,it is more convenient to find a better epsilon by choosing  theta
            delete_idx = [delete_idx i];
        end
    end
    
    %
    if(k > 1)
        idx_sys{k} = setdiff(delete_idx, deleted_idx'); %   only save the kth delete_idx
        
    else
        idx_sys{k} = delete_idx;
    end
    
    result.idx(k,1:size(idx_sys{k},2)) = idx_sys{k}; %save kth idx_sys
    result.idx_sys{k} = idx_sys{k};
    
    y(delete_idx) = 0;
    Phi(delete_idx, :) = 0;
    
    
    no_delete_idx = find(y~=0); %save all the deleted_idx
    result.no_delet_idx = no_delete_idx;
    
    if size(no_delete_idx',2)~=0
        [~,    alone]= find_idx( no_delete_idx',0);%if alone = 1,it mearns the length of the longest consecutive sequence is less than the threshold
    end
    
    result.deleted_idx = find(y==0);
    
    deleted_idx = find(y==0); %save all the deleted_idx
    size(deleted_idx )
    result.deleted_size(k) = size(deleted_idx ,1);
    fprintf( 2,'After the identification of %dth subsystem,there are %d elements have been deleted.\n\n\n\n',k,size(deleted_idx,1 ));
    if size(deleted_idx,1 )==size(y,1) || (k>1&&result.deleted_size(k-1) == size(deleted_idx ,1))||alone==1%if all idx have been deleted ,break the loop and save the number of subsystems
        
        result.s = k;
        break;
    end

    result.s = k;
end



