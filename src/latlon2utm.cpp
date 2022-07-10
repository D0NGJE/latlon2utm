#include <iostream>
#include <ros/ros.h>
#include <iomanip>
#include <sensor_msgs/NavSatFix.h>
#include <geometry_msgs/PoseStamped.h>
#include <string>
#include <cmath>
#include <vector>

class LatLon2UTM{
    public:
    LatLon2UTM();

    private:
    ros::NodeHandle nh_;
    ros::Subscriber gps_sub_;    

    void wgs2utm(double lat, double lon, int zone,double& east, double& north); // WGS84 to UTM
    void utm2wgs(double east, double north, int zone, double& lat, double& lon); //UTM to WGS84
    void gps_callback(const sensor_msgs::NavSatFix::ConstPtr& gps_msg);
};

LatLon2UTM::LatLon2UTM() : nh_("~"){
    gps_sub_ = nh_.subscribe("/gps_data", 1, &LatLon2UTM::gps_callback, this);
}

void LatLon2UTM::wgs2utm(double lat, double lon, int zone , double& east, double& north){
    double lat_rad = lat * M_PI/180;
    double lon_rad = lon * M_PI/180;

    double phi = lat_rad;
    double lambda = lon_rad;
    double lambda0 = (zone * 6 -183) * M_PI/180;
    double sm_a = 6378137;
    double sm_b = 6356752.31;

    double ep2 = (pow(sm_a, 2.0) - pow(sm_b, 2.0)) / pow(sm_b, 2.0);
    double nu2 = ep2*pow(cos(phi), 2.0);
    double N = pow(sm_a, 2.0) / (sm_b * sqrt(1 + nu2));
    double l = lambda - lambda0;
    double t = tan(phi);
    double t2 = t * t;

    double l3coef = 1 - t2 + nu2;
    double l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2);
    double l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2 - 58.0 * t2 * nu2;
    double l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2 - 330.0 * t2 * nu2;
    double l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2);
    double l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2);

    east = N * cos(phi) * l + 
        (N / 6.0 * pow(cos(phi), 3.0) * l3coef * pow(l, 3.0)) + 
        (N / 120.0 * pow(cos(phi), 5.0) * l5coef * pow(l, 5.0)) + 
        (N / 5040.0 * pow(cos(phi), 7.0) * l7coef * pow(l, 7.0));

    double n = (sm_a - sm_b) / (sm_a + sm_b);
    double alpha = ((sm_a + sm_b) / 2.0) * (1.0 + (pow(n, 2.0) / 4.0) + (pow(n,4.0) / 64.0));
    double beta = (-3.0 * n / 2.0) + (9.0 * pow(n, 3.0) / 16.0) + (-3.0 * pow(n, 5.0) / 32.0);
    double gamma = (15.0 * pow(n, 2.0) / 16.0) + (-15.0 * pow(n,4.0) / 32.0);
    double delta = (-35.0 * pow(n,3.0) / 48.0) + (105.0 * pow(n, 5.0) / 256.0);
    double epsilon = (315.0 * pow(n, 4.0) / 512.0);

    double ArcLengthMeridian = alpha * (phi + (beta * sin(2.0 * phi)) + (gamma * sin(4.0 * phi)) + (delta * sin(6.0  * phi)) + (epsilon * sin(8.0 * phi)));

    north = ArcLengthMeridian + 
            (t / 2.0 * N * pow(cos(phi), 2.0) * pow(l, 2.0)) + 
            (t / 24.0 * N * pow(cos(phi), 4.0) * l4coef * pow(l, 4.0)) + 
            (t / 720.0 * N * pow(cos(phi), 6.0) * l6coef * pow(l, 6.0)) + 
            (t / 40320.0 * N * pow(cos(phi), 8.0) * l8coef * pow(l, 8.0));
}

void LatLon2UTM::utm2wgs(double east, double north, int zone, double& lat, double& lon){
    double equatorial_radius = 6378137;
    double polar_radius = 6356752.3141;
    double flattening_3d = 0.00167922039780297994;
    double scale_factor = 0.9996;
    double polar_axis = 6356752.3141;
    double equ_rad = 6378137;
    char Dir = 'N';
    double false_easting = 500000;
    double central_meridian = 129;

    double n = flattening_3d;
    double a = equatorial_radius;
    double k0 = scale_factor;
    double e = sqrt(1 - (polar_radius / equatorial_radius) * (polar_radius / equatorial_radius));
    double xi_north = 0;

    double beta1 = (1/2)*n-(2/3)*pow(n,2)+(37/96)*pow(n,3) - (1/360)*pow(n,4) - (81/512)*pow(n,5) + (96199/604800)*pow(n,6) - (5406467/38707200)*pow(n,7) + (7944359/67737600)*pow(n,8) - (7378753979/97542144000)*pow(n,9) + (25123531261/804722688000)*pow(n,10);
    double beta2 = (1/48)*pow(n,2) + (1/15)*pow(n,3) - (437/1440)*pow(n,4) + (46/105)*pow(n,5) - (1118711/3870720)*pow(n,6) + (51841/1209600)*pow(n,7) + (24749483/348364800)*pow(n,8) - (115295683/1397088000)*pow(n,9) + (5487737251099/51502252032000)*pow(n,10);
    double beta3 = (17/480)*pow(n,3) - (37/840)*pow(n,4) - (209/4480)*pow(n,5) + (5569/90720)*pow(n,6) + (9261899/58060800)*pow(n,7) - (6457463/17740800)*pow(n,8) + (2473691167/9289728000)*pow(n,9) - (852549456029/20922789888000)*pow(n,10);
    double beta4 = (4397/161280)*pow(n,4) - (11/504)*pow(n,5) - (830251/7257600)*pow(n,6) + (466511/2494800)*pow(n,7) + (324154477/7664025600)*pow(n,8) - (937932223/3891888000)*pow(n,9) - (89112264211/5230697472000)*pow(n,10);
    double beta5 = (4583/161280)*pow(n,5) - (108847/3991680)*pow(n,6) - (8005831/63866880)*pow(n,7) + (22894433/124540416)*pow(n,8) + (112731569449/557941063680)*pow(n,9) - (5391039814733/10461394944000)*pow(n,10);
    double beta6 = (20648693/638668800)*pow(n,6) -  (16363163/518918400)*pow(n,7) - (2204645983/12915302400)*pow(n,8) + (4543317553/18162144000)*pow(n,9) + (54894890298749/167382319104000)*pow(n,10);
    double beta7 = (219941297/5535129600)*pow(n,7) - (497323811/12454041600)*pow(n,8) - (79431132943/332107776000)*pow(n,9) + (4346429528407/12703122432000)*pow(n,10);
    double beta8 = (191773887257/3719607091200)*pow(n,8) -  (17822319343/336825216000)*pow(n,9) - (497155444501631/1422749712384000)*pow(n,10);
    double beta9 = (11025641854267/158083301376000)*pow(n,9)  - (492293158444691/6758061133824000)*pow(n,10);
    double beta10 = (7028504530429620/72085985427456000)*pow(n,10);
    double AA = (a/(1+n))*(1+(1/4)*pow(n,2)+(1/64)*pow(n,4)+(1/256)*pow(n,6)+(25/16384)*pow(n,8)+(49/65536)*pow(n,10));

    if(Dir == 0){
        xi_north = north/(k0*AA);
    }
    else{
        xi_north = (10000000 - north)/(k0*AA);
    }

    double eta_east = (east - false_easting)/(k0*AA);
    double xi_prime = xi_north - (beta1*sin(2*xi_north)*cosh(2*eta_east)+beta2*sin(4*xi_north)*cosh(4*eta_east)+beta3*sin(6*xi_north)*cosh(6*eta_east)+beta4*sin(8*xi_north)*cosh(8*eta_east)+beta5*sin(10*xi_north)*cosh(10*eta_east)+beta6*sin(12*xi_north)*cosh(12*eta_east)+beta7*sin(14*xi_north)*cosh(14*eta_east));
    double eta_prime = eta_east - (beta1*cos(2*xi_north)*sinh(2*eta_east)+beta2*cos(4*xi_north)*sinh(4*eta_east)+beta3*cos(6*xi_north)*sinh(6*eta_east)+beta4*cos(8*xi_north)*sinh(8*eta_east)+beta5*cos(10*xi_north)*sinh(10*eta_east)+beta6*cos(12*xi_north)*sinh(12*eta_east)+beta7*cos(14*xi_north)*sinh(14*eta_east));
    double tau_prime = sin(xi_prime)/sqrt(pow(sinh(eta_prime),2)+pow(cos(xi_prime),2));
    double longr = atan(sinh(eta_prime)/cos(xi_prime));
    double Tau0 = tau_prime;
    double sigma0 = sinh(e*atanh(e*Tau0/sqrt(1+pow(Tau0,2))));
    double f_tau0 = Tau0*sqrt(1+pow(sigma0,2))-sigma0*sqrt(1+pow(Tau0,2))-tau_prime;
    double df_tau_dtau0 = (sqrt((1+pow(sigma0,2))*(1+pow(Tau0,2)))-sigma0*Tau0)*(1-pow(e,2))*sqrt(1+pow(Tau0,2))/(1+(1-pow(e,2))*pow(Tau0,2));
    double Tau1 = Tau0-f_tau0/df_tau_dtau0;
    double sigma1 = sinh(e*atanh(e*Tau1/sqrt(1+pow(Tau1,2))));
    double f_tau1 = Tau1*sqrt(1+pow(sigma1,2))-sigma1*sqrt(1+pow(Tau1,2))-tau_prime;
    double df_tau_dtau1 = (sqrt((1+pow(sigma1,2))*(1+pow(Tau1,2)))-sigma1*Tau1)*(1-pow(e,2))*sqrt(1+pow(Tau1,2))/(1+(1-pow(e,2))*pow(Tau1,2));
    double Tau2 = Tau1-f_tau1/df_tau_dtau1;
    double sigma2 = sinh(e*atanh(e*Tau2/sqrt(1+pow(Tau2,2))));
    double f_tau2 = Tau2*sqrt(1+pow(sigma2,2))-sigma2*sqrt(1+pow(Tau2,2))-tau_prime;
    double df_tau_dtau2 = (sqrt((1+pow(sigma2,2))*(1+pow(Tau2,2)))-sigma2*Tau2)*(1-pow(e,2))*sqrt(1+pow(Tau2,2))/(1+(1-pow(e,2))*pow(Tau2,2));
    double Tau3 = Tau2-f_tau2/df_tau_dtau2;
    double sigma3 = sinh(e*atanh(e*Tau3/sqrt(1+pow(Tau3,2))));
    double f_tau3 = Tau3*sqrt(1+pow(sigma3,2))-sigma3*sqrt(1+pow(Tau3,2))-tau_prime;
    double df_tau_dtau3 = (sqrt((1+pow(sigma3,2))*(1+pow(Tau3,2)))-sigma3*Tau3)*(1-pow(e,2))*sqrt(1+pow(Tau3,2))/(1+(1-pow(e,2))*pow(Tau3,2));
    double Tau4 = Tau3-f_tau3/df_tau_dtau3;
    double sigma4 = sinh(e*atanh(e*Tau4/sqrt(1+pow(Tau4,2))));
    double f_tau4 = Tau4*sqrt(1+pow(sigma4,2))-sigma4*sqrt(1+pow(Tau4,2))-tau_prime;
    double df_tau_dtau4 = (sqrt((1+pow(sigma4,2))*(1+pow(Tau4,2)))-sigma4*Tau4)*(1-pow(e,2))*sqrt(1+pow(Tau4,2))/(1+(1-pow(e,2))*pow(Tau4,2));
    double Tau5 = Tau4-f_tau4/df_tau_dtau4;
    double sigma5 = sinh(e*atanh(e*Tau5/sqrt(1+pow(Tau5,2))));
    double f_tau5 = Tau5*sqrt(1+pow(sigma5,2))-sigma5*sqrt(1+pow(Tau5,2))-tau_prime;
    double df_tau_dtau5 = (sqrt((1+pow(sigma5,2))*(1+pow(Tau5,2)))-sigma5*Tau5)*(1-pow(e,2))*sqrt(1+pow(Tau5,2))/(1+(1-pow(e,2))*pow(Tau5,2));
    double Tau6 = Tau5-f_tau5/df_tau_dtau5;
    double sigma6 = sinh(e*atanh(e*Tau6/sqrt(1+pow(Tau6,2))));
    double f_tau6 = Tau6*sqrt(1+pow(sigma6,2))-sigma6*sqrt(1+pow(Tau6,2))-tau_prime;

    lon = central_meridian + (longr * 180 / M_PI);
    lat = atan(Tau4) * 180 / M_PI;
}

void LatLon2UTM::gps_callback(const sensor_msgs::NavSatFix::ConstPtr& gps_msg){
    double lat = gps_msg->latitude;
    double lon = gps_msg->longitude;
    double east, north;
    double re_lat, re_lon;
    int zone;
    if(lon < 0){
        zone = (lon + 180) / 6 + 1;
    }
    else{
        zone = lon / 6 + 31;
    }
    wgs2utm(lat, lon, zone, east, north);
    utm2wgs(east, north, zone, re_lat, re_lon);

    double easting = east * 0.9996 + 500000;
    double northing = north * 0.9996;

    std::cout << std::setprecision(9);
    std::cout << "lat : " << lat << std::endl;
    std::cout << "lon : " << lon << std::endl;
    std::cout << "zone : " << zone << std::endl;
    std::cout << "E : " << easting << std::endl;
    std::cout << "N : " << northing << std::endl;
    std::cout << "re lat : " << re_lat << std::endl;
    std::cout << "re lon : " << re_lon << std::endl << std::endl;
    
}

int main(int argc, char** argv){
    ros::init(argc, argv, "latlon2utm");
    LatLon2UTM node;
    ros::spin();
    return 0;
}