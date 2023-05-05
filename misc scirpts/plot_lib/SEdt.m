% Finit diff i rummet, fjärde ordningen på andra derivatan

function [dpsi] = SEdt(psi, U, hbar, m, dx)
    
    Nx = length(psi);
    C = 1i*hbar/(2*m);
    
    dpsi = zeros(Nx,1);
    
    dpsi(2) = C * (psi(1)-2*psi(2)+psi(3))/(dx^2) - 1i*U(2)*psi(2)/hbar;
    dpsi(Nx-1) = C * (psi(Nx-2)-2*psi(Nx-1)+psi(Nx))/(dx^2) - ... 
        1i*U(Nx-1)*psi(Nx-1)/hbar;
    
    for i=3:Nx-2
      dpsi(i)= C * (-psi(i+2)+16*psi(i+1)-30*psi(i)+16*psi(i-1)-psi(i-2))/(12*dx^2) - ... 
          1j*U(i)*psi(i)/hbar;
    end
end