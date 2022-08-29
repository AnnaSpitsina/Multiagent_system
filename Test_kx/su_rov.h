#ifndef SU_ROV_H
#define SU_ROV_H

#include <QObject>
#include <QVector>
#include "qkx_coeffs.h"
#include "kx_protocol.h"
#include <QTimer>

#define ANPA_MOD_CNT 24

extern double X[2000][2];
extern QVector<double> K;

struct InitData {
    double m;
    double Fa;
    double Farx[3];
    double cv1[4];
    double cv2[4];
    double cw1[4];
    double cw2[4];
    double lambda[6][6];
    double Ta[6][8];
    double J[3];
    double kd;
    double h[4]; //радиус-вектор координат центра водоизмещения НПА в связанной СК, [1]-x [2]-y [3]-z
    double Td;
    double depth_limit;
    double max_depth;
}; //struct InitData

class SU_ROV : public QObject {
    Q_OBJECT
public:
    explicit SU_ROV(QObject *parent = 0);
    virtual ~SU_ROV();
signals:

public slots:
private:
    void start();

public:
    void model(const float Upl,const float Upp,const float Usl,const float Usp, const float Uzl, const float Uzp);
    void runge(const float Upl,const float Upp,const float Usl,const float Usp, const float Uzl, const float Uzp,const float Ttimer,const float dt=0.01);

    double a[ANPA_MOD_CNT];
    double da[ANPA_MOD_CNT];
    //константы
    double m;
    double Fa;
    double Farx[3];
    double g;
    double G;
    double cv1[4];
    double cv2[4];
    double cw1[4];
    double cw2[4];
    double lambda[7][7];
    double Ta[7][9];
    double C[7][7];
    double Vt[7]; //влиянние течения
    double Wv[7]; //вектор силы моментов, вызванных внешними возмущениями
    double J[4];
    double kd;
    double h[4]; //радиус-вектор координат центра водоизмещения НПА в связанной СК [1]-x [2]-y [3]-z
    double Td;
    double depth_limit;
    double max_depth;
    //переменные
    double sumX, sumZ;
    double cur_depth, Wx, Wy, Wz;
    double Psi_g, Gamma_g, Tetta_g;

    double Psi_gi, W_Psi_g, W_Gamma_g, W_Tetta_g;
    int N;
    double deltaSx, deltaSz;

    double Ppl, Ppp, Psl, Psp, Pzl, Pzp;
    double Ppl_x, Ppp_x, Psl_x, Psp_x, Pzl_x, Pzp_x;
    double Ppl_y, Ppp_y, Psl_y, Psp_y, Pzl_y, Pzp_y;
    double Ppl_z, Ppp_z, Psl_z, Psp_z, Pzl_z, Pzp_z;
    double Upl, Upp, Usl, Usp, Uzl, Uzp; //напряжения движителей

    double FloatageX, FloatageY, FloatageZ, Fdx, Fdy, Fdz, Fgx, Fgy, Fgz, Fcx, Fcy, Fcz;
    double Mdx, Mdy, Mdz, Mgx, Mgy, Mgz, Mcx, Mcy, Mcz;
    double Mpl_x, Mpp_x, Msl_x, Msp_x, Mzl_x, Mzp_x;
    double Mpl_y, Mpp_y, Msl_y, Msp_y, Mzl_y, Mzp_y;
    double Mpl_z, Mpp_z, Msl_z, Msp_z, Mzl_z, Mzp_z;
    double Max,May,Maz; // моменты от силы Архимеда

    double x_global, y_global, z_global;
    double vx_local,  vy_local, vz_local;  //lineinye skorosti SPA v svyazannyh osyah
    double vx_global, vy_global, vz_global;

public:
    void resetModel();
    void tick(const float Upl,const float Upp,const float Usl,const float Usp, const float Uzl, const float Uzp,const float Ttimer);
    float Fx,Fy,Fz; //total forces for XYZ-axis
    float Mx,My,Mz; //total moments for XYZ-axis
protected:
    Qkx_coeffs * K_protocol=nullptr;
    x_protocol * X_protocol=nullptr;
    QTimer timer;
};

#endif // SU_ROV_H
