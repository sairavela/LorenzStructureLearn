function bag = Lorenz_terms(xnp)

    x2 = [xnp(1,:).*xnp(2,:); xnp(1,:).*xnp(3,:); xnp(2,:).*xnp(3,:); xnp(1,:).^2; xnp(2,:).^2; xnp(3,:).^2];
    x1 = [xnp; x2];
    bag = x1; 
end