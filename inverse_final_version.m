clear
clc

%initial joint configuration
q0=[pi -pi/2 pi/2]';

% weighted matrix
%Please input flag when use W;
W=[4 0 0;0 1 0;0 0 1];
[qq]=ik(q0,2,W);  

%%% inverse loop
%%% q0 is initial joint configuration
%%% flag == 1   std solution (default method)
%%%      == 2   maxmize manipulability
%%%      == 3   limit joints' angle
%%% W weighted matrix default as eye(3)
%%% Please set a flag when use W
%%% plot_3r working for simulation is located at the end of the loop
%%% Please uncomment plot_error when plot norm(error)
%%% the yellow circle is the desired trajactory
%%% the black one is the real trajactory

% e=plot_manipubility(qq,qqq);
      
function [q_buffer]=ik(q0,flag,W)

    if nargin<2
        flag=1;
    end
    if nargin <3
        W=eye(3);
    end

    %initial joint configuration  and tip position with respect to base.
    q=q0;
    
    %initialize count,timer and delta t
    count=0;
    t=0;
%     dt=0.001; %1ms
    dt=0.1;  %100ms
    
    %break when erro < zero
    zero=1e-6;
    % default gain matrix
%     K=diag([500,500]);
    
    %gain for dt=0.1s
    K=diag([5,5]);
    
    % store error, tip position and joint configuration per loop
    error=[];
    tip_buffer=[];
    q_buffer=[];
    
    %inverse loop
    while true
       
        if count> 5000 %no more than 5s 
            break
        end
       %current tip position
       tip=fk(q);
       
       %store current tip position and joint configuration
       tip_buffer=[tip_buffer;tip'];
       q_buffer=[q_buffer;q'];
       
        %desired tip position with repect to base
        xd=[0.25*(1-cos(pi*t));0.25*(1-sin(pi*t))]+[1;1];
       
        
        %compute error
        e=xd-tip;
        if norm(e)<zero
            break;
        end
        error=[error;norm(e)];

        %desired tip velocity
        vd=[0.25*pi*sin(pi*t);-0.25*pi*cos(pi*t)];
        
        %compute total Cartisian velocity;
        v_total = vd+K*e;
        %compute jacobian
        J=Jacobian(q);
        %compute pesudo inverse
        Jpinv=rpinv(J,W);
        %compute null space
        JN=null(J);
        
        %compute joint velocity 
        %%% Please change the gain matrix manually when you change dt
        if flag == 1
        %%% std solution
        vq=Jpinv*v_total;  
        elseif flag == 2
        %%% manipulability solution
        vq=Jpinv*v_total+2.5*JN*[0;sin(q(2))*cos(q(2));sin(q(3))*cos(q(3))];
        elseif flag == 3
        %%% joints limited solution
        %joints limits  min     max
        q1lim=pi/2;   % -pi/2 ~ 0
        q2lim=pi;     % 0       pi
        q3lim=pi;     % 0       pi

        vq=Jpinv*v_total+2.5*JN*[-(q(1)+pi/4)/q1lim/q1lim;-(q(2)-pi/2)/q2lim/q2lim;-(q(3)-pi/2)/q3lim/q3lim];    
        end
        %update joint position
        q=q+vq*dt;
        
        %update count and timer
        count=count+1;
        t=t+dt;
    end
    loop=count;
    
    %simulation 
    %tracking trajactory
     plot_3r(tip_buffer,q_buffer);
    function f=plot_3r(tip_buffer,q_buffer)
        %initialize
        f = figure;
        axis equal
        axis([-3 3 -3 3])
        hold on
        
        %prepare for ellipse
        ttt = linspace(0, 2*pi, 50);
        % unit circle before transformation
        unit_circle = [cos(ttt);sin(ttt)];
        
        % tip trajactory
        traj=animatedline('color','k','LineWidth',3);
        %desired trajactory
        circle=animatedline('color','y','LineWidth',1); 
        
        for j=1:loop   
        %trajactory tracking
        figure(f)
        q1=q_buffer(j,1);
        q2=q_buffer(j,2);
        q3=q_buffer(j,3);
        tipx=tip_buffer(j,1);
        tipy=tip_buffer(j,2);
        link1 =line([0 cos(q1)],[0 sin(q1)],'Marker','o','LineWidth',1.5,'color','r');
        link2 =line([cos(q1) cos(q1)+cos(q1+q2)],[sin(q1) sin(q1)+sin(q1+q2)],'Marker','o','LineWidth',1.5,'color','g');
        link3 = line([cos(q1)+cos(q1+q2) cos(q1)+cos(q1+q2)+cos(q1+q2+q3)],[sin(q1)+sin(q1+q2) sin(q1)+sin(q1+q2)+sin(q1+q2+q3)],'Marker','o','LineWidth',1.5,'color','b');
        
        cx=0.25*(1-cos(pi*dt*j))+1;
        cy=0.25*(1-sin(pi*dt*j))+1;
        
        addpoints(circle,cx,cy);
        addpoints(traj,tipx,tipy);
        
        ellp=plot_ellipse(q_buffer(j,:),tip_buffer(j,:),f);
        drawnow
        
        
        pause(dt);
        delete(link1);
        delete(link2);
        delete(link3);
        delete(ellp);
        end
        
        function h=plot_ellipse(q,tip,f)
            figure(f);
            J_temp = Jacobian(q);
            JJ= inv(J_temp*J_temp');
            
            center=[tip(1),tip(2)];
            ellipse=sqrtm(JJ)*unit_circle;
            
            x = ellipse(1,:)+center(1); 
            z = ellipse(2,:)+center(2);
            h=plot(x,z);
        end
    end


    %plot error to time
%     plot_error(error);
    function f=plot_error(error)        
        %initialize
        %please change title and label manully
        f=figure;
        axis([0 5 0 1e-4]);
        title('1000% Gain')
        xlabel('t');
        ylabel('|error|');
        err=animatedline('color','b');

        %plot norm(error)
        for i=1:loop
            figure(f);
            addpoints(err,i/1000,error(i));
        end
    end
end
      
        function e=plot_manipubility(qq,qqq)
            e=figure;
            hold on
            title('manipulability measure')
            std=animatedline('color','b');
            max=animatedline('color','g');
            legend('std','max')
            n=size(qq,1);
            for i=1:n
                J_std=Jacobian(qq(i,:));
                m=sqrt(det(J_std*J_std'));
                addpoints(std,i/1000,m);
                
                J_max=Jacobian(qqq(i,:));
                mm=sqrt(det(J_max*J_max'));
                addpoints(max,i/1000,mm);
            end
        end
  %plot joint angles 
  % qq for standard q_buffer
  % qqq for joints limited q_buffer
    function f=plot_joints(qq,qqq)
        f1=figure;
        axis([0 5 -3.5 3.5]);
        title('joint 3')
        xlabel('t');
        ylabel('rad/s');
        joint_std=animatedline('color','b');
        joint_limited=animatedline('color','y');
        
%         joint1=animatedline('color','r');       
%         joint2=animatedline('color','g');
%         joint3=animatedline('color','b');

        maxlim=animatedline('LineStyle',':');
        minlim=animatedline('LineStyle',':');
        legend('std','limited');
        for i=1:5001
            figure(f1);
            addpoints(joint_std,i/1000,qq(i,1));  
            addpoints(joint_limited,i/1000,qqq(i,1));
            
%             addpoints(joint3,i/1000,qq(i,3));
            addpoints(maxlim,i/1000,pi);
            addpoints(minlim,i/1000,0);
        end
    end

 
%%%forward kinematic
function tip = fk(q)

xtip = cos(q(1))+cos(q(1)+q(2))+cos(q(1)+q(2)+q(3));

ytip = sin(q(1))+sin(q(1)+q(2))+sin(q(1)+q(2)+q(3));

tip=[xtip;ytip];
end

%%%Jacobian
function J = Jacobian(q)

J11=-sin(q(1))-sin(q(1)+q(2))-sin(q(1)+q(2)+q(3));
J12=-sin(q(1)+q(2))-sin(q(1)+q(2)+q(3));
J13=-sin(q(1)+q(2)+q(3));

J21=cos(q(1))+cos(q(1)+q(2))+cos(q(1)+q(2)+q(3));
J22=cos(q(1)+q(2))+cos(q(1)+q(2)+q(3));
J23=cos(q(1)+q(2)+q(3));

J=[J11 J12 J13;
   J21 J22 J23];
end

%%% pseudo inverse based on SVD
function Jpinv=rpinv(J,W)
    if nargin<2
        [U,S,V] = svd(J);
        T=S;
        T(find(S~=0)) = 1./S(find(S~=0));
        Jpinv = V*T'*U';
    else
    %when w is a diagnalized,this comment should be more efficient
    %     V=W;
    %     V(find(W~=0)) = 1./W(find(W~=0));

        Jpinv =inv(W)*J'*inv(J*inv(W)*J');
    end
end

%%%caculate nullspace of J
%%% Please check homework 2 for derivation.
function N=null(J)
    [U,S,V]=svd(J);
    N=eye(3)-V*[1 0 0;0 1 0;0 0 0]*V';
end
