#include "su_rov.h"
#include <cmath>
#include <iostream>



SU_ROV::SU_ROV(QObject *parent) : QObject(parent)
{
    X_protocol = new x_protocol("kx_pult.conf", "x",X);
    K_protocol = new Qkx_coeffs("kx_pult.conf","k");

    X[1][0]=32;
    connect(&timer, &QTimer::timeout,[this](){
        X[2][0]=K[32];
        this->tick(0,0,0,0,0,0,0.01);
    });
    resetModel();
    timer.start(100);

    m = 20;
    cv1[1] = 10.9; cv1[2] = 95.0; cv1[3] = 63.3;
    cv2[1] = 10.9; cv2[2] = 114; cv2[3] = 76;
    cw1[1] = 228.6; cw1[2] = 366; cw1[3] = 366; // kak v rabote Egorova
    cw2[1] = 2.29; cw2[2] = 36.6; cw2[3] = 36.6;
    //Vt[1] = 1; Vt[2] = 1; Vt[3] = 1; Vt[4] = 0; Vt[5] = 0; Vt[6] = 0; // скорость течения
    //Wv[1] = 0; Wv[2] = 0; Wv[3] = 0; Wv[4] = 0; Wv[5] = 0; Wv[6] = 0; //внешние возмущения, лин. скорости([1]-[3], угловые скорости - [4]-[6])
    //h[1]= ; h[2]= ; h[3]= ; // радиус-вектор координат центра водоизмещения
    lambda[1][1] = 50; lambda[2][2] = 101; lambda[3][3] = 101;
    lambda[4][4] = 50; lambda[5][5] = 50; lambda[6][6] = 50;
    Ta[1][1] = 0; Ta[1][2] = 0; Ta[1][3] = 1; Ta[1][4] = 1; Ta[1][5] = 0; Ta[1][6] = 0;
    Ta[2][1] = 0.766; Ta[2][2] = -0.766; Ta[2][3] = 0; Ta[2][4] = 0; Ta[2][5] = -0.766; Ta[2][6] = 0.766;
    Ta[3][1] = -0.64; Ta[3][2] = -0.64; Ta[3][3] = 0; Ta[3][4] = 0; Ta[3][5] = -0.64; Ta[3][6] = -0.64;
    Ta[4][1] = 117.7; Ta[4][2] = -117.7; Ta[4][3] = 0; Ta[4][4] = 0; Ta[4][5] = 6.28; Ta[4][6] = -6.28;
    Ta[5][1] = 112.5; Ta[5][2] = 112.5; Ta[5][3] = 0; Ta[5][4] = 0; Ta[5][5] = -112.5; Ta[5][6] = -112.5;
    Ta[6][1] = 134; Ta[6][2] = -134; Ta[6][3] = 118.5; Ta[6][4] = -118.5; Ta[6][5] = 134; Ta[6][6] = -134;
    //матрица сил и моментов инерции (проверить вторую матрицу, пока я вбила из мат модели, но кажется там я ошиблась она же не симметрична относительно оси, что странно)
    C[1][1] = 0; C[1][2] = (m+lambda[2][2])*a[20]; C[1][3] = -(m + lambda[3][3])*a[19]; C[1][4] = 0; C[1][5] = 0; C[1][6] = 0;
    C[2][1] = -(m + lambda[1][1])*a[20]; C[2][2] = 0; C[2][3] = (m + lambda[3][3])*a[18]; C[2][4] = 0; C[2][5] = 0; C[2][6] = 0;
    C[3][1] = (m + lambda[1][1])*a[19]; C[3][2] = -(m+lambda[2][2])*a[18]; C[3][3] = 0; C[3][4] = 0; C[3][5] = 0; C[3][6] = 0;
    C[4][1] = 0; C[4][2] = 0; C[4][3] = 0; C[4][4] = 0; C[4][5] = -(J[3]+lambda[6][6])*a[20]; C[4][6] = (J[2]+lambda[5][5])*a[19];
    C[5][1] = 0; C[5][2] = 0; C[5][3] = 0; C[5][4] = (J[3]+lambda[6][6])*a[20]; C[5][5] = 0; C[5][6] = -(J[1]+lambda[4][4])*a[18];
    C[6][1] = 0; C[6][2] = 0; C[6][3] = 0; C[6][4] = -(J[2]+lambda[5][5])*a[19]; C[6][5] = (J[1]+lambda[4][4])*a[18]; C[6][6] = 0;
    J[1] = 4; J[2] = 19.8; J[3] = 19.8; //moment inercii apparata vdol sootvetstvuushih osei
    kd = 3; //koefficient usilenija dvizhitelei
    Td = 0.15; //postojannaya vremeni dvizhitelei
    depth_limit=50;
    max_depth=50;
}

void SU_ROV::model(const float Upl,const float Upp,const float Usl,const float Usp, const float Uzl, const float Uzp) {
    int limit1, limit2;
    double G;

    //модули упоров движителей
    Ppl = a[7];  // передний левый(1)
    Ppp = a[8];  // передний правый(2)
    Psl = a[9];  // стредний левый(3)
    Psp = a[10];  // средний правый(4)
    Pzl = a[11];  // задний левый(5)
    Pzp = a[12];  // задний правый(6)

    //проекции упоров движителей на продольную ось апарата X
    Ppl_x = Ppl*Ta[1][1];
    Ppp_x = Ppp*Ta[1][2];
    Psl_x = Psl*Ta[1][3];
    Psp_x = Psp*Ta[1][4];
    Pzl_x = Pzl*Ta[1][5];
    Pzp_x = Pzp*Ta[1][6];

    //проекции упоров движителей на продольную ось апарата Y
    Ppl_y = Ppl*Ta[2][1];
    Ppp_y = Ppp*Ta[2][2];
    Psl_y = Psl*Ta[2][3];
    Psp_y = Psp*Ta[2][4];
    Pzl_y = Pzl*Ta[2][5];
    Pzp_y = Pzp*Ta[2][6];

    //проекции упоров движителей на продольную ось апарата Z
    Ppl_z = Ppl*Ta[3][1];
    Ppp_z = Ppp*Ta[3][2];
    Psl_z = Psl*Ta[3][3];
    Psp_z = Psp*Ta[3][4];
    Pzl_z = Pzl*Ta[3][5];
    Pzp_z = Pzp*Ta[3][6];

    //момент создаваемый движетельным комплексом вокруг оси X
    Mpl_x = Ppl*Ta[4][1];
    Mpp_x = Ppp*Ta[4][2];
    Msl_x = Psl*Ta[4][3];
    Msp_x = Psp*Ta[4][4];
    Mzl_x = Pzl*Ta[4][5];
    Mzp_x = Pzp*Ta[4][6];

    //момент создаваемый движетельным комплексом вокруг оси Y
    Mpl_y = Ppl*Ta[5][1];
    Mpp_y = Ppp*Ta[5][2];
    Msl_y = Psl*Ta[5][3];
    Msp_y = Psp*Ta[5][4];
    Mzl_y = Pzl*Ta[5][5];
    Mzp_y = Pzp*Ta[5][6];

    //момент создаваемый движетельным комплексом вокруг оси Z
    Mpl_z = Ppl*Ta[6][1];
    Mpp_z = Ppp*Ta[6][2];
    Msl_z = Psl*Ta[6][3];
    Msp_z = Psp*Ta[6][4];
    Mzl_z = Pzl*Ta[6][5];
    Mzp_z = Pzp*Ta[6][6];

    double g = 9.81;
    G = m*g; //вес аппарата
    Fa = 200;
    Farx[0] = 0; Farx[1] = 0; Farx[2] = -Fa;

    //obnulenie verticalnoi polozhitelnoi skorosti apparata pri dostizhenii poverhnosti
    limit1 = limit2 = 0;
    if (a[17] >= max_depth) {
      a[17] = max_depth;
        if (a[3] <= 0) {
          a[3] = 0;
          limit1 = 1;
      }
    };

    //obnulenie verticalnoi polozhitelnoi skorosti apparata pri dostizhenii dna
    if (a[17] <= 0)
    {
      a[17] = 0;
        if (a[3] >= 0)
      {
          a[3] = 0;
          limit2 = 1;
      }
    };

    Fdx = Ppl_x + Ppp_x + Psl_x + Psp_x + Pzl_x + Pzp_x; // вектор сил и моментов, создаваемых движительным комплексом
    Fgx = -cv1[1] * a[1] * fabs(a[1]) - cv2[1] * a[1]; //произведение D1*Vx
    FloatageX = -sin(a[5]) * (G + Farx[2]);
    Fcx = C[1][1]*a[1] + C[1][2]*a[2]+C[1][3]*a[3]+C[1][4]*a[18]+C[1][5]*a[19] + C[1][6]*a[20];
    //FloatageX = 0; //обнуление плавучести
    da[1] = (1/(m + lambda[1][1])) * (Fdx + Fgx + Fcx + FloatageX + Wv[1]); //vx'

    Fdy = Ppl_y + Ppp_y + Psl_y + Psp_y + Pzl_y + Pzp_y; // вектор сил и моментов, создаваемых движительным комплексом
    Fgy = -cv1[2] * a[2] * fabs(a[2]) - cv2[2] * a[2]; //произведение D1*Vy
    FloatageY = cos(a[5]) * sin(a[4]) * (G + Farx[2]);
    Fcy = C[2][1]*a[1] + C[2][2]*a[2]+C[2][3]*a[3]+C[2][4]*a[18]+C[2][5]*a[19] + C[2][6]*a[20];
    //FloatageY = 0; //обнуление плавучести
    da[2] = (1/(m + lambda[2][2])) * (Fdy + Fgy + Fcy + FloatageY + Wv[2]); //vy'

    Fdz = Ppl_z + Ppp_z + Psl_z + Psp_z + Pzl_z + Pzp_z; // вектор сил и моментов, создаваемых движительным комплексом
    Fgz = -cv1[3] * a[3] * fabs(a[3]) - cv2[3] * a[3]; //произведение D1*Vz
    FloatageZ = cos(a[4]) * cos(a[5]) * (G + Farx[2]);
    Fcz = C[3][1]*a[1] + C[3][2]*a[2]+C[3][3]*a[3]+C[3][4]*a[18]+C[3][5]*a[19] + C[3][6]*a[20];
    //FloatageZ = 0; //обнуление плавучести
    da[3] = (1/(m + lambda[3][3])) * (Fdz + Fgz + Fcz + FloatageZ + Wv[3]); //vz'

// da[4-6] -> производная угла крена, дифферента, курса
//следующие 3 уравнения это Кинематические уравнения для углов Эйлера-Крылова
//описывающее преобразование вектора угловых скоростей относительно осей НПА Ox,Oy,Oz в вектор
//угловых скоростей  по курсу, дифференту и крену соответственно.

    da[4] = a[18] + (1/cos(a[5]) * ((a[19]) * sin(a[4]) * sin(a[5])  + sin(a[5]) * cos(a[4]) * a[20])) + Vt[4];  //proizvodnaya krena

    da[5] = a[19] * cos(a[4]) - sin(a[4]) * a[20] + Vt[5];  //proizvodnaya differenta

    da[6] = (1/cos(a[5])) * (a[19] * sin(a[4]) + cos(a[4]) * (a[20])) + Vt[6]; //proizvodnaya kursa
 // Из матмодели имеем
 //K_двi - усредненный коэффициент усиления i-го движителя; T_двi=J_i/K_v1i  – наибольшее значение постоянной времени i-го ВМА
    da[7] = (1/Td) * (kd * (double)Upl - Ppl);  // передний нижний правый(1)
    da[8] = (1/Td) * (kd * (double)Upp - Ppp);  // передний нижний левый(2)
    da[9] = (1/Td) * (kd * (double)Usl - Psl);  // задний нижний левый(3)
    da[10] = (1/Td) * (kd * (double)Usp - Psp); //задний нижний правый(4)
    da[11] = (1/Td) * (kd * (double)Uzl - Pzl); // передний верхний правый(5)
    da[12] = (1/Td) * (kd * (double)Uzp - Pzp); // передний верхний левый(6)

    double alfa[4][4]; //матрица перевода из связанной СК в глобальную СК
    alfa[1][1] = cos(a[5])*cos(a[6]);
    alfa[2][1] = sin(a[6])*cos(a[5]);
    alfa[3][1] = -sin(a[5]);
    alfa[1][2] = cos(a[6])*sin(a[5])*sin(a[4])-cos(a[4])*sin(a[6]);
    alfa[2][2] = cos(a[6])*cos(a[4])+sin(a[4])*sin(a[5])*sin(a[6]);
    alfa[3][2] = sin(a[4])*cos(a[5]);
    alfa[1][3] = sin(a[6])*sin(a[4])+cos(a[6])*cos(a[4])*sin(a[5]);
    alfa[2][3] = sin(a[5])*sin(a[6])*cos(a[4])-cos(a[6])*sin(a[4]);
    alfa[3][3] = cos(a[5])*cos(a[4]);

    da[15] = alfa[1][1] * a[1] + alfa[1][2] * a[2] + alfa[1][3] * a[3] + Vt[1];
    //dx_global

    da[16] = alfa[2][1] * a[1] + alfa[2][2] * a[2] + alfa[2][3] * a[3] + Vt[2];
    //dy_global

    da[17] = alfa[3][1] * a[1] + alfa[3][2] * a[2] + alfa[3][3] * a[3] + Vt[3];
    //dz_global

    double Fax = -sin(a[5])*Fa;
    double Fay = sin(a[4])*cos(a[5])*Fa;
    double Faz = cos(a[5])*cos(a[4])*Fa;

    Mdx = Mpl_x + Mpp_x + Msl_x + Msp_x + Mzl_x + Mzp_x;
    Mgx = -cw1[1] * a[18] * fabs(a[18]) - cw2[1] * a[18];
    Max = -h[2]*Faz + h[3]*Fay;
    //Max = 0; //obnulenie momenta ot sily Arhimeda
    Mcx = C[4][1]*a[1] + C[4][2]*a[2]+C[4][3]*a[3]+C[4][4]*a[18]+C[4][5]*a[19] + C[4][6]*a[20];
    da[18] = (1/(J[1] + lambda[4][4])) * (Mdx + Mcx + Mgx + Max + Wv[4]);

    Mdy = Mpl_y + Mpp_y + Msl_y + Msp_y + Mzl_y + Mzp_y;
    Mgy = -cw1[2] * a[19] * fabs(a[19]) - cw2[2] * a[19];
    May = -Faz*h[1] + Fax*h[3];
    //May = 0; //obnulenie momenta ot sily Arhimeda
    Mcy = C[5][1]*a[1] + C[5][2]*a[2]+C[5][3]*a[3]+C[5][4]*a[18]+C[5][5]*a[19] + C[5][6]*a[20];
    da[19] = (1/(J[2] + lambda[5][5])) * (Mdy + Mcy + Mgy + May + Wv[5]);

    Mdz = Mpl_z + Mpp_z + Msl_z + Msp_z + Mzl_z + Mzp_z;
    Mgz = -cw1[3] * a[20] * fabs(a[20]) - cw2[3] * a[20];
    Maz = -h[1]*Fay + h[2]*Fax;
    //Maz = 0; //obnulenie momenta ot sily Arhimeda
    Mcz = C[6][1]*a[1] + C[6][2]*a[2]+C[6][3]*a[3]+C[6][4]*a[18]+C[6][5]*a[19] + C[6][6]*a[20];
    da[20] = (1/(J[3] + lambda[6][6])) * (Mdz + Mcz + Mgz + Maz + Wv[6]);

    da[21] = a[1];
    da[22] = a[2];
    da[23] = a[3];

}

void SU_ROV::resetModel(){
    for (int i=0;i<ANPA_MOD_CNT;i++) {a[i] = 0.0f; da[i]=0.0f;}   //f на конце означает число с плавающей точкой
    for (int i=0; i<7;i++){
        Wv[i]=0;
        Vt[i]=0;
        h[i]=0;   //потом исправить на реальное значение сверху
    }
}

void SU_ROV::tick(const float Upl,const float Upp,const float Usl,const float Usp, const float Uzl, const float Uzp,const float Ttimer){

    runge(Upl, Upp, Usl, Usp, Uzl, Uzp,Ttimer,Ttimer);
}

SU_ROV::~SU_ROV(){

}

void SU_ROV::runge(const float Upl,const float Upp,const float Usl,const float Usp, const float Uzl, const float Uzp, const float Ttimer, const float dt) {
    const double Kc = 180/M_PI;
    double a1[24], y[24];
    int i;
    const double H1 = dt;
    const int n = ANPA_MOD_CNT;
    model(Upl, Upp, Usl, Usp, Uzl, Uzp);
    for (i = 1; i < n; i++) {
      a1[i] = a[i];
      y[i] = da[i];
      a[i] = a1[i] + 0.5 * H1 * da[i];
    }

    model(Upl, Upp, Usl, Usp, Uzl, Uzp);
    for (i = 1; i < n; i++)
    {
      y[i] = y[i]+ 2 * da[i];
      a[i] = a1[i] + 0.5 * H1 * da[i];
    }

    model(Upl, Upp, Usl, Usp, Uzl, Uzp);
    for (i = 1; i < n; i++) {
      y[i] = y[i] + 2 * da[i];
      a[i] = a1[i] + H1 * da[i];
    }

    model(Upl, Upp, Usl, Usp, Uzl, Uzp);
    for (i = 1; i < n; i++) {
      a[i] = a1[i] + (H1 / 6) * (y[i] + da[i]);
    }

    //данные в СУ ( с преобразованием координат)

    x_global = a[15]; //koordinata apparata v globalnoi SK
    y_global = a[16];  //koordinaty apparata v globalnoi SK (преобразование координат)
    z_global = a[17]; //otstojanie ot dna otnositelno repernoi tochki, kotoraja na dne
    cur_depth = max_depth + z_global;  //tekush"aya glubina SPA
    Wx = a[18] * Kc; //uglovye skorosti SPA v svyazannyh osyah v gradus/sekunda
    Wy = a[19] * Kc;
    Wz = a[20] * Kc;

    vx_local = a[1]; vy_local = a[2]; vz_local = a[3];  //lineinye skorosti SPA v svyazannyh osyah
    vx_global = da[15]; vy_global = da[16]; vz_global = da[17];  // lineinye skorosti SPA v globalnyh osyah

    Gamma_g = a[4] * Kc; // ugol krena
    Tetta_g = a[5] * Kc; // ugol differenta
    Psi_g = a[6] * Kc; // ugol kursa (преобразование координат)

    W_Gamma_g = da[4] * Kc; // proizvodnaya ugla krena
    W_Tetta_g = da[5] * Kc; // proizvodnaya ugla differenta
    W_Psi_g = da[6] * Kc; // proizvodnaya ugla kursa

    N = fabs(Psi_g / 360);
    if (Psi_g >= 360) Psi_gi = Psi_g - N * 360; // ugol kursa na indikaciu
    if (Psi_g <= -360) Psi_gi = Psi_g + N * 360;

    deltaSx = vx_local * Ttimer; //prirash"enie koordinaty X dlya SVS (v svyazannoi s SPA SK)
    sumX += deltaSx;

    deltaSz = vz_local * Ttimer; //prirash"enie koordinaty Z dlya SVS (v svyazannoi s SPA SK)
    sumZ += deltaSz;


    X[10][0]=Wx;
    X[11][0]=Wy;
    X[12][0]=Wz;

    X[13][0]=vx_local;
    X[14][0]=vy_local;
    X[15][0]=vz_local;

    X[16][0]=W_Gamma_g;
    X[17][0]=W_Tetta_g;
    X[18][0]=W_Psi_g;

    X[19][0]=x_global;
    X[20][0]=y_global;
    X[21][0]=z_global;

    X[22][0]=Ppl;
    X[23][0]=Ppp;
    X[24][0]=Psl;
    X[25][0]=Psp;
    X[26][0]=Pzl;
    X[27][0]=Pzp;

}


