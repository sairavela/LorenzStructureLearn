function [params, var_out, mderr, par_all, var_all, runT] = term_elim(pars,var_in, new_t, x0,alph0)
%% Do Parameter Estimation to further eliminate
merr = inf;
runT = 0; cxx = eye(2);mderr = inf;


par_all = zeros(length(alph0),100);
var_all = zeros(length(alph0),100);

tic
sels = zeros(1,1000); 
noise1 = randn(length(alph0),100); 
noise2 = randn(length(alph0),100);

old_t = (pars ~= 0); 
max_var = 100;
if max(var_in)~= 0
    max_var = max(var_in)
end

alph1 = repmat(pars(:),[1 100])+10*sqrt(max_var)*noise1.*new_t + 10*old_t.*noise2.*sqrt(var_in);
sel = (pars(:)~=0) + new_t;
sel0 = sel;
selN = sum(sel);
sel_i = sel.*(1:27)';
% alph1 = alph1.*repmat(sel,[1 100]);
alphst = mean(alph1,2);
runT = 0;

prev_var = var(alph1,0,2);
while (runT<500 && mderr> 1e-6 && mean(prev_var)>1e-7)
    runT = runT + 1;
    %True system
    xthis = x0+randn(length(x0),1);
    temp =Lorenz_xnp1(xthis,alph0);
    xtrue = temp(:,end);  %+randn(2,1)*1e-12;
    
    %ensemble system
    xens = repmat(xthis,[1 100]);%+randn(2,100)*1e-16;
    for ens = 1:size(alph1,2)
        temp =Lorenz_xnp1(xens(:,ens),alph1(:,ens));
        xensout(:,ens) = temp(:,end);
    end
    xensout = xensout + 1e-5*randn(3, 100);
    % Error truth
%     merr = mean(abs(mean(alph1,2)-alph0));
    % Error on performance not on parameter
    derr = (xtrue - xensout);
    mderr = (mean(abs(derr(:))));
   
    % Graphs
    subplot(121); bar(abs((mean(alph1,2)-alph0)./alph0));
    xlabel('Parameters'); ylabel('Percent Error');
    set(gca,'yscale','log');
    ylim([10e-8 10e0]);
    title('Convergence');
    drawnow;
    
    %Variable Selection and Updates
    
    errens = xtrue - xensout;
    cax = (alph1 - mean(alph1,2))*(xensout - mean(xensout,2))'/99;
    cxx = cov(xensout');
    %Need to make this  efficient
    % note -- not the enkf update
    incrmnt = cax*pinv(cxx) *errens;
%     sel = ones(size(incrmnt,1),1);
    %sel = FisherInfSel(alph1,exper);
    %sel = FisherInfSel(alph1);
    
%     sel_i_nz = nonzeros(sel_i); 
%     dense_a = alph1(sel_i_nz,:); %reshape(nonzeros(alph1), [length(nonzeros(alph1))/100,100]);
%     sel1 = MutInfSel(dense_a,errens);
%     sel2 = zeros(length(alph0),1);
%     for i = 1:length(sel1)
%         if sel1(i) ==1
%             sel2(sel_i_nz(i)) = 1;
%         end
%     end
%     sel = sel2;
%     sels(runT) = sum(sel);
    alph1 = alph1+ sel.*incrmnt;
    nused(runT) = sum(sel);
    
    %Graphs
    subplot(122); plot(alph0,'*k','LineWidth',3); hold on;
    plot(alphst,'+','LineWidth',2);
    plot(alph1,'o','MarkerSize',10,'LineWidth',2);
    %plot(alph0,'k','LineWidth',3);
    hold off;
    xlabel('Parameters'); ylabel('Solution Ensemble');
    ylim([-0.1,0.2]);
    legend('truth','First Guess','solution'); title('Ensemble Parameter Estimation');
    xlim([0 27])
    par_all(:,runT) = mean(alph1,2);
    var_all(:,runT) = var(alph1,0,2);
    prev_var = var(alph1,0,2);
end

params = mean(alph1,2);
var_out = var_all(:,runT);

end