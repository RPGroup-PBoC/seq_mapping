function val = beta_conditions(beta, E, e_bj)
    expect_e = sum(sum(e_bj.*exp(-beta.*e_bj))./sum(exp(-beta.*e_bj)));
    %expect_esq = prod(sum(e_bj.^2.*exp(-beta*e_bj)))/Z;
    %expect_desq = expect_esq - expect_e^2;
    val = expect_e - E;
    %deriv = - expect_desq;
end