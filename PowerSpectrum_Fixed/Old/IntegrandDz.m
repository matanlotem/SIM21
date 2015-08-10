function integr = IntegrandDz(a)
    global omm
    integr = a.^(3/2)./((1-omm)*a.^3+omm).^(3/2);
end