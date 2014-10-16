classdef fem2d
    

    properties(Constant, SetAccess = private)
        % Reference triangle (K_t):
        %
        %                 xi_2
        %                 ^  
        %     (-1,1)      |
        %       +         +--> xi_1 
        %       |\
        %       | \
        %       |  \
        %       |   \
        %       +----+
        % (-1,-1)     (1,-1)

        % vertices of the reference triangle
        vt_ref = [-1, -1;
                   1, -1;
                  -1,  1];
    
    end
    properties(SetAccess = private)
        degree
        xi
        quadW
        nr_qnodes
        mesh_nodes
        mesh_elem
        nr_nodes
        nr_elem
        jmat
        ijmat
        jacob
        Q_prj
        % xi   : Gauss points on K_t
        % quadW: Quadrature weights for LGL points
    end
    
    methods
        function obj = fem2d(degree, mesh_nodes, mesh_elem)

            obj.degree = degree;
            [obj.xi, obj.quadW] = fem2d.Gauss_data(degree);
            obj.nr_qnodes  = length(obj.quadW);
            obj.mesh_nodes = mesh_nodes;
            obj.mesh_elem  = mesh_elem;
            obj.nr_elem    = size(mesh_elem, 1);
            obj.nr_nodes   = size(obj.mesh_nodes, 1);
        end 

        function obj = init_solver(obj)
            obj = obj.quad_prj();
            obj = obj.calc_jacobian();
        end


        function [u_sample, x_sample] = interp_elem(obj, u, xi)
            lxi = size(xi, 1);
            S   = zeros(lxi, 3);
            S(:, 1) = fem2d.varphi1(xi);
            S(:, 2) = fem2d.varphi2(xi);
            S(:, 3) = fem2d.varphi3(xi);
            if (length(u) ~= obj.nr_nodes),
                error('length of input vector is incorrect.')
            end
            u_sample = zeros( obj.nr_elem*lxi, 1);
            x_smaple = zeros( obj.nr_elem*lxi, 2);
            index = (1:lxi);
            for i=1:obj.nr_elem,
                i1 = obj.mesh_elem(i,1);
                i2 = obj.mesh_elem(i,2);
                i3 = obj.mesh_elem(i,3);

                x_sample(index, :) =  S(:, 1)*obj.mesh_nodes(i1, :)...
                                    + S(:, 2)*obj.mesh_nodes(i2, :)...
                                    + S(:, 3)*obj.mesh_nodes(i3, :); 
                                    
                u_sample(index, 1) =  S(:, 1)*u(i1)...
                                    + S(:, 2)*u(i2)...
                                    + S(:, 3)*u(i3); 

                index = lxi+index;
            end
        end

        function r = norm_L2(obj, u)
            [u_s, ~] = interp_elem(obj, u, obj.xi);
            u_sabs = u_s.*conj(u_s);
            index = (1:obj.nr_qnodes);
            r = 0;
            for i=1:obj.nr_elem,
                r = r + sum(obj.quadW.*u_sabs(index,1))*obj.jacob(i);
                index = obj.nr_qnodes + index;
            end
            r = sqrt(r);
        end

    end


    methods(Access = private)

        function obj = quad_prj(obj)
            obj.Q_prj = zeros(3, obj.nr_qnodes);
            obj.Q_prj(1, :) = transpose(obj.quadW.*fem2d.varphi1(obj.xi));
            obj.Q_prj(2, :) = transpose(obj.quadW.*fem2d.varphi2(obj.xi));
            obj.Q_prj(3, :) = transpose(obj.quadW.*fem2d.varphi3(obj.xi));
        end

        function obj = calc_jacobian(obj)
            obj.jmat  = zeros(obj.nr_elem, 4);
            obj.ijmat = zeros(obj.nr_elem, 4);
            obj.jacob = zeros(obj.nr_elem, 1);

            for i=1:obj.nr_elem,
                x1 = obj.mesh_nodes(obj.mesh_elem(i,1),:);
                x2 = obj.mesh_nodes(obj.mesh_elem(i,2),:);
                x3 = obj.mesh_nodes(obj.mesh_elem(i,3),:);
                J1 = 0.5*(x2-x1);
                J2 = 0.5*(x3-x1); 
                M = [J1', J2'];
                iM = inv([J1', J2']);
                obj.jmat(i,:)  = M(:); 
                obj.ijmat(i,:) = iM(:);
                obj.jacob(i)   = abs(det(M));
            end
        end
    end    

    methods(Static=true, Access = private)
    % Static methods do not require a specific instance of the class.
    % These ordinary functions associated with the class.

        function [xi, quadW] = Gauss_data(degree)
        % Gauss quadrature data on the reference triangular element (K_t).
        % xi    - Gauss quadrature nodes (in K_t)
        % quadW - associated weights

            nr_nodes = rule_full_size(degree);
            [ xi, quadW ] = triasymq( degree,...
                                      fem2d.vt_ref(1, :),...
                                      fem2d.vt_ref(2, :),...
                                      fem2d.vt_ref(3, :),...
                                      nr_nodes);
            xi    = transpose(xi);
            quadW = transpose(quadW);
        end

        function r = varphi1(xi)
            r = -(xi(:, 1)+xi(:, 2))/2;
        end
        function r = varphi2(xi)
            r = (1+xi(:, 1))/2;
        end
        function r = varphi3(xi)
            r = (1+xi(:, 2))/2;
        end

    end
    methods(Static=true)
        function [p, t] = rect_mesh(centroid, dim, h)
            p0 = centroid -0.5*dim;
            dx = [h,0];
            dy = [0,h];
            m  = 1+dim(1)/h;
            n  = 1+dim(2)/h;
            index = (0:m-1);
            p = zeros(m*n, 2);
            shift = 1;
            for j=0:n-1
                p(index+shift,1) = p0(1)+index*dx(1);
                p(index+shift,2) = p0(2)+index*dx(2);
                shift = shift + m;
                p0 = p0 + dy;
            end
            t = zeros(2*(m-1)*(n-1),3);
            for j=0:n-2
                for i=0:m-2
                    t(2*(i+(m-1)*j)+1,1) = i+m*j+1;
                    t(2*(i+(m-1)*j)+1,2) = i+1+m*j+1;
                    t(2*(i+(m-1)*j)+1,3) = i+m*(j+1)+1;

                    t(2*(i+(m-1)*j)+2,1) = i+1+m*j+1;
                    t(2*(i+(m-1)*j)+2,2) = i+1+m*(j+1)+1;
                    t(2*(i+(m-1)*j)+2,3) = i+m*(j+1)+1;
                end
                
            end

        end



    end
end
