classdef test_fem2d < matlab.unittest.TestCase
% test = test_fem1d
% test.run
% run(test, 'test_GP_ABC1a_CQ(testCase')

    methods(Test)
        function test_interpolation1(testCase)
            tol = 1e-6;
            h = 1/100;
            dim = [4,2];
            centroid = [0,0];
            [p,t] = fem2d.rect_mesh(centroid, dim, h);
            degree = 4;
            obj = fem2d(degree, p, t);
            
            obj = obj.init_solver();
            mypoly = @(x,y) x.^2+y.^2;
            u = mypoly(p(:,1), p(:,2));
            
            [u_sample, x_sample] = obj.interp_elem(u, obj.xi);
            u_actual = mypoly(x_sample(:,1), x_sample(:,2));
            e = norm(u_sample-u_actual)/norm(u_actual);
            disp(e)
            testCase.verifyLessThan( e, tol );
        end
        function test_interpolation2(testCase)
            tol = 1e-6;
            load('mesh_data_rect')

            degree = 4;
            obj = fem2d(degree, p, t);
            
            obj = obj.init_solver();
            gaussian = @(x,y)exp(-2*x.^2-4*y.^2);
            u = gaussian(p(:,1), p(:,2));
            
            [u_sample, x_sample] = obj.interp_elem(u, obj.xi);
            u_actual = gaussian(x_sample(:,1), x_sample(:,2));
            e = norm(u_sample-u_actual)/norm(u_actual);
            testCase.verifyLessThan( e, tol );
        end
        function test_interpolation3(testCase)
            tol = 1e-6;
            load('mesh_data')
            degree = 4;
            obj = fem2d(degree, p, t);
            
            obj = obj.init_solver();
            gaussian = @(x,y)exp(-2*x.^2-4*y.^2);
            u = gaussian(p(:,1), p(:,2));
            
            [u_sample, x_sample] = obj.interp_elem(u, obj.xi);
            u_actual = gaussian(x_sample(:,1), x_sample(:,2));
            e = norm(u_sample-u_actual)/norm(u_actual);
            testCase.verifyLessThan( e, tol );
        end
        function test_norm1(testCase)
            tol = 1e-6;
            h = 1/100;
            dim = [4,2];
            centroid = [0,0];
            [p,t] = fem2d.rect_mesh(centroid, dim, h);

            degree = 4;
            obj = fem2d(degree, p, t);
            
            obj = obj.init_solver();
            gaussian = @(x,y)exp(-2*x.^2-4*y.^2);
            u = gaussian(p(:,1), p(:,2));
            
            u_norm = obj.norm_L2(u);
            u_norm_exact = sqrt(pi/sqrt(32));
            e = norm(u_norm-u_norm_exact)/norm(u_norm_exact);
            testCase.verifyLessThan( e, tol );
        end
        function test_norm2(testCase)
            tol = 1e-6;
            load('mesh_data_best')
            degree = 4;
            obj = fem2d(degree, p, t);
            
            obj = obj.init_solver();
            a = 8; b = 8;
            gaussian = @(x,y)exp(-a*x.^2-b*y.^2);
            u = gaussian(p(:,1), p(:,2));
            %trimesh(t, p(:,1), p(:,2), u);

            u_norm = obj.norm_L2(u);
            u_norm_exact = sqrt(pi/sqrt(4*a*b));
            e = norm(u_norm-u_norm_exact)/norm(u_norm_exact);
            testCase.verifyLessThan( e, tol );
        end
    end
end

