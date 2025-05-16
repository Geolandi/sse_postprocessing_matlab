function [m,Cm] = invert_comp(G,d,Cd,m0,Cm0)
    G_t = G';
    B = inv(G*Cm0*G_t + Cd);
    m = m0 + (Cm0*G_t / (G*Cm0*G_t + Cd) * (d - G*m0));
    Cm = Cm0 - Cm0 * G_t * B * G * Cm0;

end