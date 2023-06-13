function [eigenvalues] = model(Thetas)

%MODEL: Eigenvalues of 2 x 2 square matrix

eigenvalues(:,1) = 0.5*((Thetas(:,1) + 2*Thetas(:,2)) + (Thetas(:,1).^2 + 4*(Thetas(:,2).^2)).^0.5);

eigenvalues(:,2) = 0.5*((Thetas(:,1) + 2*Thetas(:,2)) - (Thetas(:,1).^2 + 4*(Thetas(:,2).^2)).^0.5);

end