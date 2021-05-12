classdef SetMembership < handle
    properties
        Omega
        W
        theta_hat
        theta_bounds
        AB_hat
        np
        ny
        AB0
        ABi
    end
    
    methods
        function obj = SetMembership(Omega,W, ABi, AB0)
            obj.Omega = Omega;
            obj.W = W;
            obj.ABi = ABi;
            obj.AB0 = AB0;
            obj.np = size(obj.ABi,3);
            obj.ny = size(obj.ABi,1);
            aux = obj.Omega.outerApprox;
            obj.theta_hat = sum(aux.A.*aux.b)'/2;
        end
        
        function [Omega, D] = update(obj,xp,x,u)
            phixu = obj.get_phixu(x,u);
            
            % Compute non falsified set
            D = Polyhedron(-obj.W.A*phixu, obj.W.b - obj.W.A*(xp - obj.AB0*[x;u]));
            
            % Compute intersection
            if D.isBounded() %&& obj.Omega.intersect(D)
                Omega = obj.Omega.intersect(D);
                
                obj.Omega = Omega.minHRep;
            else
                Omega = obj.Omega;
            end
            % Point Estimate (center of bounding box)
            aux = Omega.outerApprox;
            bounds_mat = aux.A.*aux.b;
            obj.theta_hat = sum(bounds_mat)'/2;
            obj.theta_bounds = nonzeros(bounds_mat);
        end
        
        function phixu = get_phixu(obj,x,u)
            phixu = [];
            % define phi(x,u) = [[A_1, B_1][x;u],... ]
            for i=1:obj.np
                phixu = [phixu, squeeze(obj.ABi(:,:,i))*[x;u]];
            end
        end
        
        function AB = get_AB(obj)

            %aux = obj.Omega.chebyCenter;
            %obj.theta_hat = aux.x;
            AB = obj.AB0;
            for i=1:obj.np
                AB = AB + obj.ABi(:,:,i)*obj.theta_hat(i);
            end
            obj.AB_hat = AB;
        end
        
        function W_theta = estimateW_theta(obj,X,U)
            aux = X*U; %product polyhedron
            Vxu = aux.V; %vertices
            W_theta = Polyhedron(zeros(0,obj.ny),[]);
            for i=1:size(Vxu,1)
                phixuv(:,:,i) = obj.get_phixu(Vxu(i,1:obj.ny)',Vxu(i,obj.ny+1:end)');
                aux = phixuv(:,:,i)*(obj.Omega - obj.theta_hat);
                W_theta = Polyhedron('V',[W_theta.V; aux.V]).minVRep;
            end
        end
        
        function reset(obj,Omega)
            obj.Omega = Omega;
            
            % Point Estimate (center of bounding box)
            aux = obj.Omega.outerApprox;
            obj.theta_hat = sum(aux.A.*aux.b)'/2;
        end
    end
end

