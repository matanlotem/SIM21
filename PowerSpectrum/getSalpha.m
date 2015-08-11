%
function Sa = getSalpha(Tk,z)

Sa=exp(-.803*(Tk).^(-2/3)*(1e-6*getTauGPn(z)).^(1/3));

end