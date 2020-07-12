function [t_out, mean_err] = term_select(alph1, elim, alph0, x0)
    %params are the current parameters (can all be zero)
    %elim is a matrix where =1 if a term has been killed by enkf
    
    %True system
    xthis = x0 + 2*randn(length(x0),10000);
    temp =Lorenz_xnp1(xthis,alph0);
    
    temp2 = Lorenz_xnp1(xthis,alph1);
    xterms = Lorenz_terms(xthis);
    xtrue = temp;
    
    elim2 = elim + (alph1~=0);
    
    %how many terms left in the "bag"
    t_curr = length(alph1) -sum((alph1(:)~=0))-sum(elim(:));

    %Zero system gives error equal to true system
    %Find MI between terms and error

    errens = temp2 -temp;%-mean(xtrue,2);
    dterms = xterms - mean(xterms,2);
    cae = dterms*(errens-mean(errens,2))'/(length(dterms-1));
    evar = var(errens,0,2); 
    tvar = var(xterms,0,2);
    rh = cae./sqrt(tvar*evar');

    mutinf = -0.5*log(1-rh.^2);
    mutinf = mutinf';
    [vs,vi] = sort(mutinf(:),'descend')
    vis = zeros(length(vs),1);
    vss = zeros(length(vs),1);
    c = 1;
    for i = vi'
        if elim2(i) == 0
            vis(c) = vi(c);
            vss(c) = vs(c);
        end
        c = c+1;
    end
    vss = nonzeros(vss);
    vis = nonzeros(vis);
    nn = length(vss);
    svs = sum(vss);
    vsum = ((svs-cumsum(vss))/svs)+((1:nn)'/nn);
    n =find(vsum==min(vsum));
%    sel = zeros(size(alph1,1),1)
%    sel(vi(1:n)) = 1;
%    sel = sel.*(1:27)';
%    t_out = nonzeros(sel);

%     [vs,vi]
    out = zeros(ceil(t_curr/2),1);
    c = 0;
    i = 0;
    alph1s = alph1(:);
    elims = elim(:);
    while (i < length(alph1(:))) && (c<n)
        i = i+1;
        j = vi(i);
        if (alph1s(j)==0) && (elims(j) == 0)
            c = c+1;
            out(i) = j;
        end

    end
    t_out = nonzeros(out);
    mean_err = mean(abs(errens),2);

end