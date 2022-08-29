#ifndef K_DESCRIPTION_H
#define K_DESCRIPTION_H



enum KDescription {

    k_kren = 4,
    k_dif,
    k_kurs,

    Ppnp = 7,
    Ppnl,
    Pznp,
    Pznl,
    Ppvp,
    Ppvl,
    Pzvl,
    Pzvp,

    Umarsh_limit = 15, //f Ограничение максимального сигнала СУ по маршу
    Upsi_limit, //f Ограничение максимального сигнала СУ по курсу
    Uteta_limit, //f Ограничение максимального сигнала СУ по дифференту
    Ugamma_limit, //f Ограничение максимального сигнала СУ по крену

    Wx = 19,
    Wy,
    Wz,

    Upnp_limit = 23, //f Ограничение максимального напряжения на ВМА
    Upnl_limit, //f Ограничение максимального напряжения на ВМА
    Uznp_limit, //f Ограничение максимального напряжения на ВМА
    Uznl_limit, //f Ограничение максимального напряжения на ВМА
    Upvp_limit,
    Upvl_limit,
    Uzvl_limit,
    Uzvp_limit,

    Upnp_scale, //f Коэффициент масштабирующий сигнал из СУ в напряжение на ВМА МВЛ
    Upnl_scale, //f Коэффициент масштабирующий сигнал из СУ в напряжение на ВМА МНЛ
    Uznp_scale, //f Коэффициент масштабирующий сигнал из СУ в напряжение на ВМА МВП
    Uznl_scale, //f Коэффициент масштабирующий сигнал из СУ в напряжение на ВМА МНП
    Upvp_scale,
    Upvl_scale,
    Uzvl_scale,
    Uzvp_scale,

    Kurs_ruchnoi_scale=40, //f Коэффициент усиления по курсу в ручном режиме
    Kurs_otladka=43, //n Коэффициент для настройки контура скорости
    Kurs_K1, //f Коэффициент К1 контура курса
    Kurs_K2, //f Коэффициент К2 контура курса
    Limit_Upsi=49, //f Ограничение максимального управляющего сигнала по курсу





};

#endif // K_DESCRIPTION_H
