#include <QCoreApplication>
#include "su_rov.h"

double X[2000][2];

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    SU_ROV rov_Model;
    return a.exec();
}
