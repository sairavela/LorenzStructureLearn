clear all;
Lorenz_params;
hold off;

%for exper = 1:6

% True parameters
alph0 = [W(:)];

%Guess a parameter ensemble
%At first all terms will be 0 
I = eye(2);
alph1 = zeros(length(alph0),100)*10;
alphst = mean(alph1,2);
runTT = 0 
vi_order = zeros(length(alph0),1);
%% Rank terms by MI between error and the terms
selT = 0;

elim = 0*alph0;
mderr = inf;
par = alph1(:,1);
c = 0;
var_in = 10*ones(length(alph0),1);

errors = zeros(10,1);
pars_all = zeros(27, 1000);
vars_all = zeros(27, 1000);
elim_all = zeros(27, 10);

while (mderr > (1e-5) && c< 5)
    [vi, mean_err] = term_select(par, elim, alph0, x0)
    if c ==0
        errors(1) = mean(abs(mean_err));
    end
    c = c+1;
    vi_in = zeros(length(alph0),1);
    for i = vi
        vi_in(i) = 1;
    end

    [par_out, var_out, mderr, par_all, var_all, runT] = term_elim(par,var_in,vi_in,x0,alph0);
    for i = 1:length(alph0)
        if (par_out(i)~=0) && (abs(par_out(i))<1e-3) && (var_out(i)<1e-6)
            elim(i) = 1;
            par_out(i) = 0;
        end
        
    end
    elim_all(:,c) = elim;
    par = par_out;
    var_in = var_out;
    errors(c+1) = mderr; 
    pars_all(:,(runTT+1):(runTT + runT)) = par_all(:,1:runT);
    vars_all(:,(runTT+1):(runTT + runT)) = var_all(:,1:runT);
    runTT = runTT + runT; 
end
    
% %True system
% xthis = x0 + randn(length(x0),1000);
% temp =Lorenz(xthis,alph0);
% xterms = Lorenz_terms(xthis);
% xtrue = temp;
% 
% %Zero system gives error equal to true system
% %Find MI between terms and error
% 
% errens = xtrue;%-mean(xtrue,2);
% dterms = xterms - mean(xterms,2);
% cae = dterms*errens'/99;
% evar = mean(errens.^2,2); %Ask question about why this differs from MulInfSel
% tvar = var(dterms,0,2);
% rh = cae./sqrt(tvar*evar');
% 
% mutinf = -0.5*log(1-rh.^2);
% [vs,vi] = sort(mutinf(:),'descend');


% %% Previous Parameter Estimation
% merr = inf;
% runT = 0; cxx = eye(2);mderr = inf;
% 
% 
% par_all = zeros(length(alph0),1000);
% var_all = zeros(length(alph0),1000);
% 
% tic
% sels = zeros(1,1000); 
% 
% while (runT<3 && mderr>1e-6)
%     runT = runT + 1;
%     
%     %True system
%     xthis = x0+randn(length(x0),1);
%     temp =Lorenz_xnp1(xthis,alph0);
%     xtrue = temp(:,end);  %+randn(2,1)*1e-12;
%     
%     %ensemble system
%     xens = repmat(xthis,[1 100]);%+randn(2,100)*1e-16;
%     for ens = 1:size(alph1,2)
%         temp =Lorenz_xnp1(xens(:,ens),alph1(:,ens));
%         xensout(:,ens) = temp(:,end);
%     end
%     % Error truth
%     merr = mean(abs(mean(alph1,2)-alph0));
%     % Error on performance not on parameter
%     derr = (xtrue - xensout);
%     mderr = (mean(abs(derr(:))));
%     
%     % Graphs
%     subplot(121); bar(abs(mean(alph1,2)-alph0));
%     xlabel('Parameters'); ylabel('Absolute Error');
%     set(gca,'yscale','log');
%     ylim([10e-16 10e0]);
%     title('Convergence');
%     drawnow;
%     
%     %Variable Selection and Updates
%     
%     errens = xtrue - xensout;
%     cax = (alph1 - mean(alph1,2))*(xensout - mean(xensout,2))'/99;
%     cxx = cov(xensout');
%     %Need to make this  efficient
%     % note -- not the enkf update
%     incrmnt = cax*pinv(cxx) *errens;
%     sel = ones(size(incrmnt,1),1);
%     %sel = FisherInfSel(alph1,exper);
%     %sel = FisherInfSel(alph1);
% %     sel = MutInfSel(alph1,errens);
%     sels(runT) = sum(sel);
%     alph1 = alph1+ sel.*incrmnt;
%     nused(runT) = sum(sel);
%     
%     %Graphs
%     subplot(122); plot(alph0,'*k','LineWidth',3); hold on;
%     plot(alphst,'+','LineWidth',2);
%     plot(alph1,'o','MarkerSize',10,'LineWidth',2);
%     %plot(alph0,'k','LineWidth',3);
%     hold off;
%     xlabel('Parameters'); ylabel('Solution Ensemble');
%     ylim([-0.1,0.2]);
%     legend('truth','First Guess','solution'); title('Ensemble Parameter Estimation');
%     xlim([0 27])
%     par_all(:,runT) = mean(alph1,2);
%     var_all(:,runT) = var(alph1,0,2);
%     
% end
% %runiter(exper) = runT;
% %end
% ttimean = toc
% text(1.5,1,sprintf('iterations=%d',runT));
% sol = mean(alph1,2);
% disp([ alphst sol alph0 std(alph1,[],2)])
% save anerondat ttimean alph1 runT;