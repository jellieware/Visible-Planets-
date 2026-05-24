#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846
#define RAD (PI / 180.0)
#define DEG (180.0 / PI)

typedef struct {
    double ra;
    double dec;
} Ephem;

typedef struct {
    double a, e, i, L, lp, ln;
    double cy_a, cy_e, cy_i, cy_L, cy_lp, cy_ln;
} Orbit;

// Explicitly segregated planetary elements. Earth/Sun Proxy separated entirely.
Orbit sun_elements = {1.00000011, 0.01671022, 0.00005, 100.46435, 102.94719, -11.26064, 0, -0.00003804, -0.01300, 35999.37244, 0.32327, -0.44541};

Orbit planets[7] = {
    {0.38709893, 0.20563069, 7.00487, 252.25084, 77.45645, 48.33167, 0, 0.00002040, -0.00594, 149472.67411, 0.15901, -0.12510}, // 0: Mercury
    {0.72333199, 0.00677323, 3.39471, 181.97973, 131.53298, 76.68069, 0, -0.00004776, -0.00079, 58517.81538, 0.00213, -0.27769}, // 1: Venus
    {1.52366231, 0.09341233, 1.85061, 355.45332, 336.04084, 49.57854, 0, 0.00011902, -0.00724, 19140.30268, 0.44388, -0.29277}, // 2: Mars
    {5.20336301, 0.04839266, 1.30530, 34.40438, 14.75385, 100.55615, 0.00060737, -0.00012880, -0.00415, 3034.74612, 0.19112, 0.20477}, // 3: Jupiter
    {9.53707032, 0.05415060, 2.48446, 49.94432, 92.43194, 113.71504, -0.00301530, -0.00036762, 0.00193, 1222.49436, -0.41897, -0.28867}, // 4: Saturn
    {19.19126393, 0.04716771, 0.76998, 313.23218, 170.96424, 74.22988, 0.00152025, -0.00019150, -0.00269, 428.48202, 0.40805, -0.20944}, // 5: Uranus
    {30.06896348, 0.00858587, 1.76917, 304.88003, 44.97135, 131.72169, -0.00125196, 0.00002514, -0.00334, 218.45945, -0.32241, -0.15049}  // 6: Neptune
};

const char* names[] = {"Mercury", "Venus", "Sun", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Moon"};

double rev360(double a) {
    double r = fmod(a, 360.0);
    return (r < 0) ? r + 360.0 : r;
}

double get_jd(int y, int m, int d, double h_ut) {
    if (m <= 2) { y--; m += 12; }
    int A = y / 100;
    int B = 2 - A + (A / 4);
    return floor(365.25 * (y + 4716)) + floor(30.6001 * (m + 1)) + d + B - 1524.5 + (h_ut / 24.0);
}

double get_gmst(double jd) {
    double T = (jd - 2451545.0) / 36525.0;
    return rev360(280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T * T - (T * T * T / 38710000.0));
}

void get_helio_p(Orbit p, double t, double *X, double *Y, double *Z) {
    double a = p.a + p.cy_a * t;
    double e = p.e + p.cy_e * t;
    double i = (p.i + p.cy_i * t) * RAD;
    double L = rev360(p.L + p.cy_L * t);
    double lp = rev360(p.lp + p.cy_lp * t) * RAD;
    double ln = rev360(p.ln + p.cy_ln * t) * RAD;

    double M = rev360(L - (lp * DEG)) * RAD;
    double E = M + e * sin(M) + (e * e / 2.0) * sin(2.0 * M);
    double xv = a * (cos(E) - e);
    double yv = a * (sqrt(1.0 - e * e) * sin(E));

    double v = atan2(yv, xv);
    double r = sqrt(xv * xv + yv * yv);
    double u = v + lp - ln;

    *X = r * (cos(ln) * cos(u) - sin(ln) * sin(u) * cos(i));
    *Y = r * (sin(ln) * cos(u) + cos(ln) * sin(u) * cos(i));
    *Z = r * (sin(u) * sin(i));
}

Ephem get_geocentric(int id, double jd) {
    double t = (jd - 2451545.0) / 36525.0;
    double obl = (23.439291 - 0.01300416 * t) * RAD;
    double ra = 0, dec = 0;

    if (id == 2) { // Sun Track
        double X, Y, Z;
        get_helio_p(sun_elements, t, &X, &Y, &Z);
        double ex = -X, ey = -Y * cos(obl) + Z * sin(obl), ez = -Y * sin(obl) - Z * cos(obl);
        ra = rev360(atan2(ey, ex) * DEG);
        dec = atan2(ez, sqrt(ex * ex + ey * ey)) * DEG;
    } else if (id == 8) { // Moon Track
        double d = jd - 2451543.5;
        double L = rev360(218.316 + 13.176396 * d);
        double M = rev360(134.963 + 13.064993 * d);
        double F = rev360(93.272 + 13.229350 * d);
        double le = L + 6.289 * sin(M * RAD);
        double lte = 5.128 * sin(F * RAD);

        double mx = cos(le * RAD) * cos(lte * RAD);
        double my = sin(le * RAD) * cos(lte * RAD);
        double mz = sin(lte * RAD);

        ra = rev360(atan2(my * cos(obl) - mz * sin(obl), mx) * DEG);
        dec = atan2(my * sin(obl) + mz * cos(obl), sqrt(mx * mx + (my * cos(obl) - mz * sin(obl)) * (my * cos(obl) - mz * sin(obl)))) * DEG;
    } else { // Clean, Isolated Planet Array Indices Lookups
        int p_idx = id;
        if (id > 2) p_idx = id - 1; // Map around Sun's display position entry

        double Xp, Yp, Zp, Xe, Ye, Ze;
        get_helio_p(planets[p_idx], t, &Xp, &Yp, &Zp);
        get_helio_p(sun_elements, t, &Xe, &Ye, &Ze);

        double gx = Xp - Xe, gy = Yp - Ye, gz = Zp - Ze;
        double ex = gx, ey = gy * cos(obl) - gz * sin(obl), ez = gy * sin(obl) + gz * cos(obl);
        ra = rev360(atan2(ey, ex) * DEG);
        dec = atan2(ez, sqrt(ex * ex + ey * ey)) * DEG;
    }
    Ephem out = {ra, dec};
    return out;
}

void unwind_coords(double *y1, double *y2, double *y3) {
    if ((*y2 - *y1) < -180.0) *y2 += 360.0;
    if ((*y2 - *y1) >  180.0) *y2 -= 360.0;
    if ((*y3 - *y2) < -180.0) *y3 += 360.0;
    if ((*y3 - *y2) >  180.0) *y3 -= 360.0;
}

double interpolate_meeus(double y1, double y2, double y3, double n) {
    double a = y2 - y1, b = y3 - y2, c = b - a;
    return y2 + (n / 2.0) * (a + b + n * c);
}

void calculate_next_events(int id, double jd_today, double lat, double lon, double tz, double current_local_hour,
                           int *rh, int *rm, int *sh, int *sm, int *r_plus, int *s_plus) {
    *rh = *rm = *sh = *sm = -1;
    *r_plus = *s_plus = 0;

    double h0 = -0.5667;
    if (id == 2) h0 = -0.8333;
    if (id == 8) h0 = 0.125;

    Ephem ep1 = get_geocentric(id, jd_today - 1.0);
    Ephem ep2 = get_geocentric(id, jd_today);
    Ephem ep3 = get_geocentric(id, jd_today + 1.0);

    unwind_coords(&ep1.ra, &ep2.ra, &ep3.ra);
    unwind_coords(&ep1.dec, &ep2.dec, &ep3.dec);

    double gmst0 = get_gmst(jd_today);
    double current_ut_hour = current_local_hour - tz;

    for (int step = 0; step < 24; step++) {
        double h_ut_start = current_ut_hour + step;
        double m0 = h_ut_start / 24.0;

        double ra  = interpolate_meeus(ep1.ra, ep2.ra, ep3.ra, m0);
        double dec = interpolate_meeus(ep1.dec, ep2.dec, ep3.dec, m0);
        double gmst = rev360(gmst0 + 360.98564736629 * m0);
        double ha = rev360(gmst + lon - ra);
        if (ha > 180.0) ha -= 360.0;
        double alt = asin(sin(lat * RAD) * sin(dec * RAD) + cos(lat * RAD) * cos(dec * RAD) * cos(ha * RAD)) * DEG;

        double h_ut_next = current_ut_hour + step + 1;
        double m1 = h_ut_next / 24.0;

        double ra_next  = interpolate_meeus(ep1.ra, ep2.ra, ep3.ra, m1);
        double dec_next = interpolate_meeus(ep1.dec, ep2.dec, ep3.dec, m1);
        double gmst_next = rev360(gmst0 + 360.98564736629 * m1);
        double ha_next = rev360(gmst_next + lon - ra_next);
        if (ha_next > 180.0) ha_next -= 360.0;
        double alt_next = asin(sin(lat * RAD) * sin(dec_next * RAD) + cos(lat * RAD) * cos(dec_next * RAD) * cos(ha_next * RAD)) * DEG;

        // Trace Rise Intersections
        if (alt < h0 && alt_next >= h0 && *rh == -1) {
            double fraction = (h0 - alt) / (alt_next - alt);
            double event_local_hour = current_local_hour + step + fraction;

            if (event_local_hour >= 24.0) {
                event_local_hour -= 24.0;
                *r_plus = 1;
            }
            if (event_local_hour < 0) event_local_hour += 24.0;

            *rh = (int)event_local_hour;
            *rm = (int)((event_local_hour - *rh) * 60.0);
        }

        // Trace Set Intersections
        if (alt >= h0 && alt_next < h0 && *sh == -1) {
            double fraction = (alt - h0) / (alt - alt_next);
            double event_local_hour = current_local_hour + step + fraction;

            if (event_local_hour >= 24.0) {
                event_local_hour -= 24.0;
                *s_plus = 1;
            }
            if (event_local_hour < 0) event_local_hour += 24.0;

            *sh = (int)event_local_hour;
            *sm = (int)((event_local_hour - *sh) * 60.0);
        }
    }
                           }

                           int main() {
                               double lat, lon;
                               int is_dst;

                               printf("==================================================\n");
                               printf("  CORRECTED PLANETARY GEOCENTRIC EPHEMERIS SYSTEM \n");
                               printf("==================================================\n");
                               printf("Enter Latitude (e.g. 40.7128): ");
                               if (scanf("%lf", &lat) != 1) return 1;
                               printf("Enter Longitude (e.g. -74.0060): ");
                               if (scanf("%lf", &lon) != 1) return 1;
                               printf("Is DST Active? (1=Yes, 0=No): ");
                               if (scanf("%d", &is_dst) != 1) return 1;

                               time_t raw = time(NULL);
                               struct tm *gmt = gmtime(&raw);
                               double hour_ut = gmt->tm_hour + (gmt->tm_min / 60.0) + (gmt->tm_sec / 3600.0);
                               struct tm *local = localtime(&raw);
                               double hour_local = local->tm_hour + (local->tm_min / 60.0) + (local->tm_sec / 3600.0);

                               double tz_offset = hour_local - hour_ut;
                               if (tz_offset > 12.0) tz_offset -= 24.0;
                               if (tz_offset < -12.0) tz_offset += 24.0;

                               double jd_now = get_jd(gmt->tm_year + 1900, gmt->tm_mon + 1, gmt->tm_mday, hour_ut);
                               double jd_today_ut = get_jd(gmt->tm_year + 1900, gmt->tm_mon + 1, gmt->tm_mday, 0.0);

                               printf("\n");
                               printf("Report Run Time: %04d-%02d-%02d @ %02d:%02d:%02d Local\n",
                                      local->tm_year + 1900, local->tm_mon + 1, local->tm_mday, local->tm_hour, local->tm_min, local->tm_sec);
                               printf("========================================================================\n");
                               printf("%-10s | %-14s | %-8s | %-9s | %-10s\n", "BODY", "HORIZON STATUS", "ALTITUDE", "NEXT RISE", "NEXT SET");
                               printf("========================================================================\n");

                               for (int i = 0; i < 9; i++) {
                                   Ephem now_ep = get_geocentric(i, jd_now);
                                   double gmst = get_gmst(jd_now);

double ha = rev360(gmst + lon - now_ep.ra);double alt = asin(sin(lat * RAD) * sin(now_ep.dec * RAD) + cos(lat * RAD) * cos(now_ep.dec * RAD) * cos(ha * RAD)) * DEG;const char* status = (alt >= 0.0) ? "ABOVE" : "BELOW";int rh, rm, sh, sm, r_plus, s_plus;calculate_next_events(i, jd_today_ut, lat, lon, tz_offset, hour_local, &rh, &rm, &sh, &sm, &r_plus, &s_plus);char r_str[16], s_str[16];if (rh == -1) sprintf(r_str, "N/A"); else sprintf(r_str, "%02d:%02d%s", rh, rm, r_plus ? "+" : "");if (sh == -1) sprintf(s_str, "N/A"); else sprintf(s_str, "%02d:%02d%s", sh, sm, s_plus ? "+" : "");printf("%-10s | %-14s | %6.2f°  | %-9s | %-10s\n", names[i], status, alt, r_str, s_str);}printf("========================================================================\n");printf("Note: '+' indicates an event occurring tomorrow (after the next local midnight).\n");return 0;}
