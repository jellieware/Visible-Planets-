#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PI 3.141592653589793238462643383279502884
#define RAD (PI / 180.0)
#define DEG (180.0 / PI)

#define AU_KM       149597870.7
#define C_AU_DAY    173.1446326847

typedef struct {
    double ra;
    double dec;
    double dist;
} Ephem;

typedef struct {
    double a, e, i, L, lp, ln;
    double ca, ce, ci;
    double cL, clp, cln;
    double b, c, s, f;
    int perturb;
} Orbit;

static const char *names[] = {
    "Mercury", "Venus", "Sun", "Mars", "Jupiter",
    "Saturn", "Uranus", "Neptune", "Moon"
};

static Orbit planets[8] = {
    {0.38709843,0.20563661,7.00559432,252.25166724,77.45771895,48.33961819,
     0.00000000,0.00002123,-0.00590158,149472.67486623,0.15940013,-0.12214182,
     0,0,0,0,0},

    {0.72332102,0.00676399,3.39777545,181.97970850,131.76755713,76.67261496,
     -0.00000026,-0.00005107,0.00043494,58517.81560260,0.05679648,-0.27274174,
     0,0,0,0,0},

    {1.00000018,0.01673163,-0.00054346,100.46691572,102.93005885,-5.11260389,
     -0.00000003,-0.00003661,-0.01337178,35999.37306329,0.31795260,-0.24123856,
     0,0,0,0,0},

    {1.52371243,0.09336511,1.85181869,-4.56813164,-23.91744784,49.71320984,
     0.00000097,0.00009149,-0.00724757,19140.29934243,0.45223625,-0.26852431,
     0,0,0,0,0},

    {5.20248019,0.04853590,1.29861416,34.33479152,14.27495244,100.29282654,
     -0.00002864,0.00018026,-0.00322699,3034.90371757,0.18199196,0.13024619,
     -0.00012452,0.06064060,-0.35635438,38.35125000,1},

    {9.54149883,0.05550825,2.49424102,50.07571329,92.86136063,113.63998702,
     -0.00003065,-0.00032044,0.00451969,1222.11494724,0.54179478,-0.25015002,
     0.00025899,-0.13434469,0.87320147,38.35125000,1},

    {19.18797948,0.04685740,0.77298127,314.20276625,172.43404441,73.96250215,
     -0.00020455,-0.00001550,-0.00180155,428.49512595,0.09266985,0.05739699,
     0.00058331,-0.97731848,0.17689245,7.67025000,1},

    {30.06952752,0.00895439,1.77005520,304.22289287,46.68158724,131.78635853,
     0.00006447,0.00000818,0.00022400,218.46515314,0.01009938,-0.00606302,
     -0.00041348,0.68346318,-0.10162547,7.67025000,1}
};

static double rev360(double x)
{
    x = fmod(x, 360.0);
    if (x < 0.0) x += 360.0;
    return x;
}

static double rev180(double x)
{
    x = rev360(x);
    if (x > 180.0) x -= 360.0;
    return x;
}

static double get_jd(int y, int m, int d, double h_ut)
{
    if (m <= 2) {
        y--;
        m += 12;
    }

    int A = y / 100;
    int B = 2 - A + (A / 4);

    return floor(365.25 * (y + 4716))
         + floor(30.6001 * (m + 1))
         + d + B - 1524.5
         + h_ut / 24.0;
}

static double centuries(double jd)
{
    return (jd - 2451545.0) / 36525.0;
}

static double get_obliquity(double T)
{
    double eps =
        23.439291111
        - 0.0130041667 * T
        - 0.000000164 * T * T
        + 0.0000005036 * T * T * T;

    return eps * RAD;
}

static double solve_kepler(double M_deg, double e)
{
    double M = rev180(M_deg) * RAD;

    double E = M + e * sin(M) * (1.0 + e * cos(M));

    for (int k = 0; k < 20; k++) {
        double f = E - e * sin(E) - M;
        double fp = 1.0 - e * cos(E);
        double delta = f / fp;

        E -= delta;

        if (fabs(delta) < 1.0e-14)
            break;
    }

    return E;
}

static void get_heliocentric(
    Orbit p, double jd,
    double *X, double *Y, double *Z)
{
    double T = centuries(jd);

    double a = p.a + p.ca * T;
    double e = p.e + p.ce * T;
    double inc = (p.i + p.ci * T) * RAD;
    double L = p.L + p.cL * T;
    double lp = p.lp + p.clp * T;
    double node = (p.ln + p.cln * T) * RAD;

    double omega = (lp - (p.ln + p.cln * T)) * RAD;

    double M = L - lp;

    if (p.perturb) {
        M += p.b * T * T
           + p.c * cos(p.f * T * RAD)
           + p.s * sin(p.f * T * RAD);
    }

    double E = solve_kepler(M, e);

    double xp = a * (cos(E) - e);
    double yp = a * sqrt(fmax(0.0, 1.0 - e * e)) * sin(E);

    double v = atan2(yp, xp);
    double r = hypot(xp, yp);
    double u = v + omega;

    *X = r * (cos(node) * cos(u) -
              sin(node) * sin(u) * cos(inc));

    *Y = r * (sin(node) * cos(u) +
              cos(node) * sin(u) * cos(inc));

    *Z = r * sin(u) * sin(inc);
}

static void ecliptic_to_equatorial(
    double x, double y, double z, double eps,
    double *X, double *Y, double *Z)
{
    *X = x;
    *Y = y * cos(eps) - z * sin(eps);
    *Z = y * sin(eps) + z * cos(eps);
}

static Ephem get_planet_ephem(int id, double jd)
{
    Ephem out = {0.0, 0.0, 0.0};

    if (id == 2) {
        double ex, ey, ez;

        get_heliocentric(planets[2], jd, &ex, &ey, &ez);

        double gx = -ex;
        double gy = -ey;
        double gz = -ez;

        double x, y, z;

        ecliptic_to_equatorial(
            gx, gy, gz,
            get_obliquity(centuries(jd)),
            &x, &y, &z
        );

        out.ra = rev360(atan2(y, x) * DEG);
        out.dec = atan2(z, hypot(x, y)) * DEG;
        out.dist = hypot(hypot(gx, gy), gz);

        return out;
    }

    int p_idx = id;
    if (id > 2) p_idx = id - 1;

    double ex, ey, ez;

    get_heliocentric(
        planets[2], jd,
        &ex, &ey, &ez
    );

    double px, py, pz;

    get_heliocentric(
        planets[p_idx], jd,
        &px, &py, &pz
    );

    double gx = px - ex;
    double gy = py - ey;
    double gz = pz - ez;

    double distance = hypot(hypot(gx, gy), gz);

    for (int k = 0; k < 4; k++) {
        double light_time = distance / C_AU_DAY;
        double planet_jd = jd - light_time;

        get_heliocentric(
            planets[p_idx], planet_jd,
            &px, &py, &pz
        );

        gx = px - ex;
        gy = py - ey;
        gz = pz - ez;

        distance = hypot(hypot(gx, gy), gz);
    }

    double x, y, z;

    ecliptic_to_equatorial(
        gx, gy, gz,
        get_obliquity(centuries(jd)),
        &x, &y, &z
    );

    out.ra = rev360(atan2(y, x) * DEG);
    out.dec = atan2(z, hypot(x, y)) * DEG;
    out.dist = distance;

    return out;
}

static Ephem get_moon_ephem(double jd)
{
    double T = centuries(jd);

    double Lp = rev360(
        218.3164477 +
        481267.88123421 * T -
        0.0015786 * T * T
    );

    double D = rev360(
        297.8501921 +
        445267.1114034 * T -
        0.0018819 * T * T
    );

    double M = rev360(
        357.5291092 +
        35999.0502909 * T -
        0.0001536 * T * T
    );

    double Mp = rev360(
        134.9633964 +
        477198.8675055 * T +
        0.0087414 * T * T
    );

    double F = rev360(
        93.2720950 +
        483202.0175233 * T -
        0.0036539 * T * T
    );

    double lon = Lp;
    double lat = 0.0;
    double distance_km = 385000.56;

    struct {
        int d, m, mp, f;
        double lon;
        double dist;
    } terms[] = {
        { 0,  0,  1,  0,  6.289, -20905.355 },
        { 2,  0, -1,  0,  1.274, -3699.111 },
        { 2,  0,  0,  0,  0.658, -2955.968 },
        { 0,  0,  2,  0,  0.214, -569.925 },
        { 0,  1,  0,  0, -0.186, 48.888 },
        { 0,  0,  0,  2, -0.114, -3.149 },
        { 2,  0, -2,  0,  0.059, 246.158 },
        { 2, -1, -1, 0,  0.057, -152.138 },
        { 2,  0,  1,  0,  0.053, -170.733 },
        { 0, -1,  1,  0,  0.046, -204.586 },
        { 4,  0, -1,  0,  0.041, -129.620 },
        { 1,  0,  0,  0,  0.041, 104.755 },
        { 0,  0,  1,  2,  0.035, -10.321 },
        { 0,  0,  1,  0,  0.031, 0.0 },
        { 4,  0, -1,  0,  0.015, -34.782 }
    };

    int count = (int)(sizeof(terms) / sizeof(terms[0]));

    for (int i = 0; i < count; i++) {
        double arg =
            (terms[i].d * D +
             terms[i].m * M +
             terms[i].mp * Mp +
             terms[i].f * F) * RAD;

        lon += terms[i].lon * sin(arg);
        distance_km += terms[i].dist * cos(arg);
    }

    lat +=
        5.128 * sin(F * RAD)
        + 0.280 * sin((Mp + F) * RAD)
        + 0.277 * sin((Mp - F) * RAD)
        + 0.173 * sin((2.0 * D - F) * RAD)
        + 0.055 * sin((2.0 * D - Mp - F) * RAD)
        + 0.046 * sin((2.0 * D - Mp + F) * RAD)
        + 0.033 * sin((2.0 * D + F) * RAD)
        + 0.017 * sin((2.0 * Mp + F) * RAD);

    double distance_au = distance_km / AU_KM;

    double x =
        distance_au *
        cos(lon * RAD) *
        cos(lat * RAD);

    double y =
        distance_au *
        sin(lon * RAD) *
        cos(lat * RAD);

    double z =
        distance_au *
        sin(lat * RAD);

    double X, Y, Z;

    ecliptic_to_equatorial(
        x, y, z,
        get_obliquity(T),
        &X, &Y, &Z
    );

    Ephem out;

    out.ra = rev360(atan2(Y, X) * DEG);
    out.dec = atan2(Z, hypot(X, Y)) * DEG;
    out.dist = distance_au;

    return out;
}

static Ephem get_geocentric(int id, double jd)
{
    if (id == 8)
        return get_moon_ephem(jd);

    return get_planet_ephem(id, jd);
}

static double get_gmst(double jd)
{
    double T = (jd - 2451545.0) / 36525.0;

    return rev360(
        280.46061837 +
        360.98564736629 * (jd - 2451545.0) +
        0.000387933 * T * T -
        T * T * T / 38710000.0
    );
}

static double get_horizon(int id)
{
    if (id == 2)
        return -0.8333;

    if (id == 8)
        return 0.125;

    return -0.5667;
}

static double get_altitude(
    int id,
    double jd,
    double lat,
    double lon)
{
    Ephem e = get_geocentric(id, jd);

    double hour_angle =
        rev180(
            get_gmst(jd) +
            lon -
            e.ra
        );

    double s =
        sin(lat * RAD) *
        sin(e.dec * RAD)
        +
        cos(lat * RAD) *
        cos(e.dec * RAD) *
        cos(hour_angle * RAD);

    if (s > 1.0) s = 1.0;
    if (s < -1.0) s = -1.0;

    return asin(s) * DEG;
}

static int find_event(
    int id,
    double start_jd,
    double lat,
    double lon,
    double target_altitude,
    int rising,
    double *event_jd)
{
    const double step = 10.0 / 1440.0;

    double t0 = start_jd;

    double a0 =
        get_altitude(
            id, t0, lat, lon
        ) - target_altitude;

    for (int k = 1; k <= 216; k++) {

        double t1 =
            start_jd + k * step;

        double a1 =
            get_altitude(
                id, t1, lat, lon
            ) - target_altitude;

        int crossing = 0;

        if (rising) {
            if (a0 <= 0.0 && a1 > 0.0)
                crossing = 1;
        } else {
            if (a0 >= 0.0 && a1 < 0.0)
                crossing = 1;
        }

        if (crossing) {

            double lo = t0;
            double hi = t1;

            for (int j = 0; j < 35; j++) {

                double mid = (lo + hi) * 0.5;

                double am =
                    get_altitude(
                        id, mid, lat, lon
                    ) - target_altitude;

                if (rising) {
                    if (am > 0.0)
                        hi = mid;
                    else
                        lo = mid;
                } else {
                    if (am < 0.0)
                        hi = mid;
                    else
                        lo = mid;
                }
            }

            *event_jd = (lo + hi) * 0.5;

            return 1;
        }

        t0 = t1;
        a0 = a1;
    }

    return 0;
}

static void jd_to_local(
    double jd,
    double timezone,
    int *hour,
    int *minute)
{
    double h =
        (
            jd -
            floor(jd - 0.5) -
            0.5
        ) * 24.0 + timezone;

    while (h < 0.0)
        h += 24.0;

    while (h >= 24.0)
        h -= 24.0;

    *hour = (int)h;

    *minute =
        (int)(
            (h - *hour) * 60.0 +
            0.5
        );

    if (*minute >= 60) {
        *minute = 0;
        (*hour)++;

        if (*hour >= 24)
            *hour = 0;
    }
}

int main(void)
{
    double lat = 42.3;
    double lon = -71.0;

    time_t raw = time(NULL);

    struct tm g = *gmtime(&raw);
    struct tm l = *localtime(&raw);

    double hour_ut =
        g.tm_hour +
        g.tm_min / 60.0 +
        g.tm_sec / 3600.0;

    double hour_local =
        l.tm_hour +
        l.tm_min / 60.0 +
        l.tm_sec / 3600.0;

    double timezone =
        hour_local - hour_ut;

    if (timezone > 12.0)
        timezone -= 24.0;

    if (timezone < -12.0)
        timezone += 24.0;

    double jd =
        get_jd(
            g.tm_year + 1900,
            g.tm_mon + 1,
            g.tm_mday,
            hour_ut
        );

    /*
     * Julian date corresponding to today's local midnight.
     *
     * This is used to determine whether a rise/set event
     * belongs to the next local calendar day.
     */
    double local_today_key =
        floor(jd + timezone / 24.0 - 0.5);

    printf("\n");

    printf(
        "Report Run Time: %04d-%02d-%02d @ %02d:%02d:%02d Local\n",
        l.tm_year + 1900,
        l.tm_mon + 1,
        l.tm_mday,
        l.tm_hour,
        l.tm_min,
        l.tm_sec
    );

    printf(
        "=========================================================\n"
    );

    printf(
        "%-10s | %-8s | %-9s | %-10s | %-10s\n",
        "BODY",
        "STATUS",
        "ALTITUDE",
        "NEXT RISE",
        "NEXT SET"
    );

    printf(
        "=========================================================\n"
    );

    for (int i = 0; i < 9; i++) {

        double altitude =
            get_altitude(
                i,
                jd,
                lat,
                lon
            );

        double rise_jd;
        double set_jd;

        int have_rise =
            find_event(
                i,
                jd,
                lat,
                lon,
                get_horizon(i),
                1,
                &rise_jd
            );

        int have_set =
            find_event(
                i,
                jd,
                lat,
                lon,
                get_horizon(i),
                0,
                &set_jd
            );

        char rise_string[32];
        char set_string[32];

        /*
         * NEXT RISE
         */
        if (have_rise) {

            int h, m;

            jd_to_local(
                rise_jd,
                timezone,
                &h,
                &m
            );

            /*
             * The event is after today's local
             * calendar day boundary.
             */
            int tomorrow =
                (floor(rise_jd + timezone / 24.0 - 0.5) >
                 local_today_key);

            sprintf(
                rise_string,
                "%02d:%02d%s",
                h,
                m,
                tomorrow ? "+" : ""
            );

        } else {

            strcpy(
                rise_string,
                "N/A"
            );
        }

        /*
         * NEXT SET
         */
        if (have_set) {

            int h, m;

            jd_to_local(
                set_jd,
                timezone,
                &h,
                &m
            );

            /*
             * The event is after today's local
             * calendar day boundary.
             */
            int tomorrow =
                (floor(set_jd + timezone / 24.0 - 0.5) >
                 local_today_key);

            sprintf(
                set_string,
                "%02d:%02d%s",
                h,
                m,
                tomorrow ? "+" : ""
            );

        } else {

            strcpy(
                set_string,
                "N/A"
            );
        }

        printf(
            "%-10s | %-8s | %7.2f° | %-10s | %-10s\n",
            names[i],
            altitude >= 0.0 ? "ABOVE" : "BELOW",
            altitude,
            rise_string,
            set_string
        );
    }

    printf(
        "=========================================================\n"
    );

    printf(
        "Note: '+' indicates an event occurring on the next local day.\n"
    );

    printf(
        "Planet model: JPL long-range elements with\n"
        "Jupiter-Saturn-Uranus-Neptune perturbation terms,\n"
        "iterative Kepler solution and planetary light-time correction.\n"
    );

    printf(
        "Rise/set: direct altitude search followed by bisection.\n\n\n"
    );

    return 0;
}
