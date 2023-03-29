
% ����˳�����ࡢ���ٶȵȳ�ʼ������
l1 = 100;  % �������� mm
l2 = 350; l3 = 350; l4 = 500;
m = 300; % ֧��A��Cˮƽ����
S0 = 400; % ֧��A��C��ֱ����





% ������˽Ƕȡ����ٶȣ��Ǽ��ٶ�
%��AB
alpha1= 0:10:360; % �Ա��� ʱ�� s
omega1=5; % �������ٶ� rad/s
acceleration1=0; % �����Ƚ��ٶ���ת
n = length(alpha1); % �Ա�������

%���廬��λ�ơ��ٶȡ����ٶ�
ySlide = zeros(1,n); % �����λ��
velocitySlide = zeros(1,n); % ������ٶ� mm/s
accelerationSlide = zeros(1,n); % �����ٶ� mm/s

for iterTime = 1:n
% ��BD
syms alpha2 ;
syms omega2;
syms acceleration2;

% ��CD
syms alpha3 ;
syms omega3;
syms acceleration3;

% ��DE
syms alpha4 ;
syms omega4;
syms acceleration4;



%�˵�λ�Ʒ���
% ʸ������(a)��l4+l3=S;
% ʸ������(b): m+l1+l2+l3=S0
    
syms S;

% ʸ�����̵�ͶӰ
% ����(a)
% ��x��
q1 = l4*cosd(alpha4)+l3*cosd(alpha3);
% % ��y��
q2 = l4*sind(alpha4)+l3*sind(alpha3)-S;
%����(b)
% ��x��
q3 = m+l1*cosd(alpha1(iterTime))+l2*cosd(alpha2)+l3*cosd(alpha3);
% ��y��
q4 = l1*sind(alpha1(iterTime))+l2*sind(alpha2)+l3*sind(alpha3)-S0;

%���λ�Ʒ��̣��ó� alpha2��alpha3��alpha4��S
%���õ�������ķ�������ǰһʱ�̵�λ����ȷ����һʱ�̵�λ��
if iterTime == 1
   x0 = [415 99 171 46];
end

   T = vpasolve (q1,q2,q3,q4,x0);
   alpha2 = T.alpha2;
   alpha3 = T.alpha3;

   alpha4 = T.alpha4;
    ySlide(iterTime) = T.S;
   if iterTime > 1
   x0 = [ySlide((iterTime-1)) alpha2 alpha3 alpha4];
   end
   
   
   
%�˵�λ�Ʒ��̶�ʱ����һ�׵�
syms V;
% ��x��
q5 = -l4*sind(alpha4)*omega4-l3*sind(alpha3)*omega3;
% ��y��
q6 = l4*cosd(alpha4)*omega4+l3*cosd(alpha3)*omega3-V;
%����(b)
% ��x��
q7 = (-1)*l1*sind(alpha1(iterTime))*omega1-l2*sind(alpha2)*omega2-l3*sind(alpha3)*omega3;
% ��y��
q8 = l1*cosd(alpha1(iterTime))*omega1+l2*cosd(alpha2)*omega2+l3*cosd(alpha3)*omega3;
%����ٶȷ��̣��ó� omega2��omega3��omega4��velocitySlide
P = vpasolve (q5,q6,q7,q8);
velocitySlide(iterTime) = double(P.V);
omega2 = P.omega2;
omega3 = P.omega3;
omega4 = P.omega4;




%�˵�λ�Ʒ��̶�ʱ������׵�
syms B;
%����(a)
% ��x��
q9 = -l4*cosd(alpha4)*(omega4)^2-l4*sind(alpha4)*acceleration4-l3*cosd(alpha3)*(omega3)^2-l3*sind(alpha3)*acceleration3;
% ��y��
q10 = -l4*sind(alpha4)*(omega4)^2+l4*cosd(alpha4)*acceleration4-l3*sind(alpha3)*(omega3)^2+l3*cosd(alpha3)*acceleration3-B;
%����(b)
% ��x��
q11 = (-1)*l1*cosd(alpha1(iterTime))*(omega1)^2-l1*sind(alpha1(iterTime))*acceleration1-l2*cosd(alpha2)*(omega2)^2-l2*sind(alpha2)*acceleration2-l3*cosd(alpha3)*(omega3)^2-l3*sind(alpha3)*acceleration3;
% ��y��
q12 = (-1)*l1*sind(alpha1(iterTime))*(omega1)^2+l1*cosd(alpha1(iterTime))*acceleration1-l2*sind(alpha2)*(omega2)^2+l2*cosd(alpha2)*acceleration2-l3*sind(alpha3)*(omega3)^2+l3*cosd(alpha3)*acceleration3;
%�����ٶȷ��̣��ó� acceleration2,acceleration3,acceleration4,accelerationSlide
Q = vpasolve (q9,q10,q11,q12);
accelerationSlide(iterTime) = Q.B;

clearvars  alpha2 alpha3 alpha4 S omega2 omega3 omega4 V acceleration2 acceleration3 acceleration4 B T P Q;
   end
figure(1)
plot(alpha1, ySlide,'-*');
title('����������� ����λ��');
xlabel('alpha1(��)');
ylabel('λ��(mm)');


figure(2)
plot(alpha1, velocitySlide,'-o');
title('����������� �����ٶ�');
xlabel('alpha1(��)');
ylabel('�ٶ�(mm/s)');

figure(3)
plot(alpha1, accelerationSlide,'-^');
title('����������� ������ٶ�');
xlabel('alpha1(��)');
ylabel('�ٶ�(mm/s^2)');