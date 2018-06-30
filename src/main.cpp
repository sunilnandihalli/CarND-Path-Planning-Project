#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
using namespace std;
#include <deque>
#include <cmath>
#include "Eigen-3.3/Eigen/Dense"
using Eigen::MatrixXd;
using Eigen::VectorXd;
vector<double> JMT(vector<double> start, vector<double> end, double T) {
  MatrixXd A = MatrixXd(3, 3);
  A << T * T * T, T * T * T * T, T * T * T * T * T,
      3 * T * T, 4 * T * T * T, 5 * T * T * T * T,
      6 * T, 12 * T * T, 20 * T * T * T;
  MatrixXd B = MatrixXd(3, 1);
  B << end[0] - (start[0] + start[1] * T + .5 * start[2] * T * T),
      end[1] - (start[1] + start[2] * T),
      end[2] - start[2];
  MatrixXd Ai = A.inverse();
  MatrixXd C = Ai * B;
  vector<double> result = {start[0], start[1], .5 * start[2]};
  for (int i = 0; i < C.size(); i++) {
    result.push_back(C.data()[i]);
  }
  return result;
}

#include <functional>
#include <cmath>
#include <fstream>
struct WP {
  double x,xdot,xdotdot;
};

std::deque<WP> WPs,WPd;
const double max_s = 6945.554;
const double mps2mph = 2.23694;
const double dt = 0.02; // time between waypoints
const double time_horizon = 3.0; // time horizon in seconds
const int num_waypoints = time_horizon / dt;
const double lane_width = 4.0;
const double max_jerk = 5; // meters per sec^3
const double max_acc = 5; // meters per sec^2
const double max_vel_mph = 45;

const double max_acc_change_per_frame = max_jerk * dt;
const double sentinel = 999999999;
double jmax = max_jerk;
double amax = max_acc;
double max_vel = max_vel_mph / mps2mph; // meters per sec

double secant(const std::function<double(double)>& f,double min,double max,double tol,const char* call_name) {
  double x0(min),x1(max);
  double f0(f(x0)),f1(f(x1));
  while(true) {
    //std::cout<<" x0 : "<<x0<<" f0 : "<<f0<<" x1 : "<<x1<<" f1 : "<<f1<<std::endl;
    double xs = (f0*x1-f1*x0)/(f0-f1);
    double fs = f(xs);
    if(fabs(fs)<tol) {
      //std::cout<<call_name<<" xs "<<xs<<" fs "<<fs<<" x0 "<<x0<<" f0 "<<f0<<" x1 "<<x1<<" f1 "<<f1<<" return xs "<<std::endl;
      return xs;
    } else {
      if(f0*fs<0) {
        f1 = fs;
        x1 = xs;
      } else if (f1*fs<0) {
        f0 = fs;
        x0 = xs;
      } else {
        if(fabs(f0)<fabs(f1)) {
	  //std::cout<<call_name<<" xs "<<xs<<" fs "<<fs<<" x0 "<<x0<<" f0 "<<f0<<" x1 "<<x1<<" f1 "<<f1<<" return x0 "<<std::endl;
	  return x0;
        } else {
	  //std::cout<<call_name<<" xs "<<xs<<" fs "<<fs<<" x0 "<<x0<<" f0 "<<f0<<" x1 "<<x1<<" f1 "<<f1<<" return x1 "<<std::endl;
          return x1;
	}
      }
    }
  }
}

// assuming maximum jerk is constant
// minimum velocity change when changing acceleration
double deltaVelocity(double a0,double a1) {
  return (a0+a1)*0.5*fabs(a1-a0)/jmax;
}

// acceleration and velocity are initial values
double dist(double v,double a,double j,double t) {
  return v*t+a*t*t*0.5+j*t*t*t/6.0;
}

// returns <distChange,velChange,jerk_fn,time>
std::tuple<double,double,std::function<double(double)>,double> distDuringAccelChangeH(double v0,double a0,double a1) {
  double j = (a1>a0?1:-1)*jmax;
  double t1 = (a1-a0)/j;
  double d = dist(v0,a0,j,t1);
  double dv = deltaVelocity(a0,a1);
  auto jfn = [t1,j](double t) { if(t>=0 && t<=t1) return j;};
  return std::make_tuple(d,dv,jfn,t1);
}

std::tuple<double,double,std::function<double(double)>,double> distDuringAccelChange(double v0,double a0,double a1) {
      auto x = distDuringAccelChangeH(v0,a0,a1);
      static char buffer[1000];
      sprintf(buffer,"distDuringAccelChange : v0 %f a0 %f a1 %f returns dist %f velChange %f deltaT %f",v0,a0,a1,std::get<0>(x),
	      std::get<1>(x),std::get<3>(x));
      //std::cout<<buffer<<std::endl;
      return x;
}


// dist,jfn,totaltime
std::tuple<double,std::function<double(double)>,double> distDuringVelocityChangeH(double a0,double v0,double v1) {
  double dv = v1-v0;
  // change in velocities
  auto direct = distDuringAccelChange(v0,a0,0);
  auto toAmax = distDuringAccelChange(v0,a0,amax);
  auto fromAmax = distDuringAccelChange(v0+std::get<1>(toAmax),amax,0);
  auto toAmin = distDuringAccelChange(v0,a0,-amax);
  auto fromAmin = distDuringAccelChange(v0+std::get<1>(toAmin),-amax,0);
  auto h = [v0,v1](std::tuple<double,double,std::function<double(double)>,double> to,
						       std::tuple<double,double,std::function<double(double)>,double> from,
											      double apeak) {
    double t1(std::get<3>(to)),t3(std::get<3>(from)),d1(std::get<0>(to)),d3(std::get<0>(from));
    double u1(v0+std::get<1>(to)),u3(v1-std::get<1>(from));
    double t2 = (u3-u1)/apeak;
    double d2 = (u1+u3)*0.5*t2;
    auto jfn1 = std::get<2>(to);
    auto jfn3 = std::get<2>(from);
    auto jfn = [jfn1,jfn3,t1,t2,t3](double t) {
      if(t>=0 && t<t1) {
	return jfn1(t);
      } else if (t<=t1+t2) {
	return 0.0;
      } else if (t<=t1+t2+t3) {
	return jfn3(t-t1-t2);
      }
    };
    if(t2>=0) {
      return std::make_tuple(d1+d2+d3,jfn,t1+t2+t3);
    } else {
      throw("error");
    }
  };
  if(fabs(dv-std::get<1>(direct))<1e-6) {
    return std::make_tuple(std::get<0>(direct),std::get<2>(direct),std::get<3>(direct));
  } else if (std::get<1>(toAmax)+std::get<1>(fromAmax)<dv) {
    return h(toAmax,fromAmax,amax);
  } else if (std::get<1>(toAmin)+std::get<1>(fromAmin)>dv) {
    return h(toAmin,fromAmin,-amax);
  } else {
    auto f = [a0,dv](double apeak) {
      return deltaVelocity(a0,apeak)+deltaVelocity(apeak,0)-dv;
    };
    static char buffer[1000];
    sprintf(buffer,"apeak : a0 %f v0 %f ",a0,v0);
    double apeak = secant(f,-amax,amax,1e-6,buffer);
    auto toApeak = distDuringAccelChange(v0,a0,apeak);
    auto fromApeak = distDuringAccelChange(v0+std::get<1>(toApeak),apeak,0);
    return h(toApeak,fromApeak,apeak);
  }
}
// dist,jfn,totaltime
std::tuple<double,std::function<double(double)>,double> distDuringVelocityChange(double a0,double v0,double v1) {
      auto x = distDuringVelocityChangeH(a0,v0,v1);
      static char buffer[1000];
      sprintf(buffer,"distDuringVelocityChange : a0 %f v0 %f v1 %f returns dist %f totaltime %f",a0,v0,v1,std::get<0>(x),std::get<2>(x));
      //std::cout<<buffer<<std::endl;
      return x;
    }
std::tuple<std::function<double(double)>,double> achieveZeroAcelAndVel(double a0,double v0,double vmin,double vmax,double delta_d) {
  auto direct = distDuringVelocityChange(a0,v0,0);
  auto toVmin = distDuringVelocityChange(a0,v0,vmin);
  auto fromVmin = distDuringVelocityChange(0,vmin,0);
  auto toVmax = distDuringVelocityChange(a0,v0,vmax);
  auto fromVmax = distDuringVelocityChange(0,vmax,0);
  double directDist = std::get<0>(direct);
  double vminDist = std::get<0>(toVmin)+std::get<0>(fromVmin);
  double vmaxDist = std::get<0>(toVmax)+std::get<0>(fromVmax);
  auto h = [delta_d](std::tuple<double,std::function<double(double)>,double> to,std::tuple<double,std::function<double(double)>,double> from,double vpeak) {
    double t3 = std::get<2>(from);
    double d3 = std::get<0>(from);
    double t1 = std::get<2>(to);
    double d1 = std::get<0>(to);
    double d2 = delta_d-d1-d3;
    double t2 = d2/vpeak;
    auto jfn3 = std::get<1>(from);
    auto jfn1 = std::get<1>(to);
    auto jfn = [t1,t2,t3,jfn1,jfn3](double t) {
      if(t>=0 && t < t1) {
	return jfn1(t);
      } else if (t<=t1+t2) {
	return 0.0;
      } else if (t<=t1+t2+t3) {
	return jfn3(t-t1-t2);
      }
    };
    return std::make_tuple(jfn,t1+t2+t3);
  };
  if(fabs(delta_d-directDist) < 1e-5 || (delta_d<vminDist && !(vmin<0)) || (delta_d>vmaxDist && !(vmax>0))) {
    return std::make_tuple(std::get<1>(direct),std::get<2>(direct));
  } else if(delta_d < vminDist) {
    if(vmin<0) {
      return h(toVmin,fromVmin,vmin);
    } else {
      throw("should not come here");
    }
  } else if(delta_d > vmaxDist) {
    if(vmax>0) {
      return h(toVmax,fromVmax,vmax);
    } else {
      throw("should not come here");
    }
  } else {
    auto f = [a0,v0,delta_d](double vpeak) {
      return std::get<0>(distDuringVelocityChange(a0,v0,vpeak))+std::get<0>(distDuringVelocityChange(0,vpeak,0))-delta_d;
    };
    static char buffer[1000];
    sprintf(buffer,"peakVel a0 %f v0 %f delta_d %f",a0,v0,delta_d);
    double vpeak = secant(f,vmin,vmax,1e-5,buffer);
    auto toVpeak = distDuringVelocityChange(a0,v0,vpeak);
    auto fromVpeak = distDuringVelocityChange(0,vpeak,0);
    return h(toVpeak,fromVpeak,vpeak);
  }
}



// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() {
  return M_PI;
}
double deg2rad(double x) {
  return x * pi() / 180;
}
double rad2deg(double x) {
  return x * 180 / pi();
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}
int ClosestWaypoint(double x,
                    double y,
                    const vector<double>& maps_x,
                    const vector<double>& maps_y) {
  double closestLen = 100000; //large number
  int closestWaypoint = 0;
  for (int i = 0; i < maps_x.size(); i++) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x, y, map_x, map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }
  }
  return closestWaypoint;
}

int NextWaypoint(double x,
                 double y,
                 double theta,
                 const vector<double>& maps_x,
                 const vector<double>& maps_y) {
  int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);
  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];
  double heading = atan2((map_y - y), (map_x - x));
  double angle = fabs(theta - heading);
  angle = min(2 * pi() - angle, angle);
  if (angle > pi() / 4) {
    closestWaypoint++;
    if (closestWaypoint == maps_x.size()) {
      closestWaypoint = 0;
    }
  }
  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x,
                         double y,
                         double theta,
                         const vector<double>& maps_x,
                         const vector<double>& maps_y) {
  int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);
  int prev_wp;
  prev_wp = next_wp - 1;
  if (next_wp == 0) {
    prev_wp = maps_x.size() - 1;
  }
  double n_x = maps_x[next_wp] - maps_x[prev_wp];
  double n_y = maps_y[next_wp] - maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];
  // find the projection of x onto n
  double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
  double proj_x = proj_norm * n_x;
  double proj_y = proj_norm * n_y;
  double frenet_d = distance(x_x, x_y, proj_x, proj_y);
  //see if d value is positive or negative by comparing it to a center point
  double center_x = 1000 - maps_x[prev_wp];
  double center_y = 2000 - maps_y[prev_wp];
  double centerToPos = distance(center_x, center_y, x_x, x_y);
  double centerToRef = distance(center_x, center_y, proj_x, proj_y);
  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }
  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; i++) {
    frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
  }
  frenet_s += distance(0, 0, proj_x, proj_y);
  return {frenet_s, frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s,
                     double d,
                     const vector<double>& maps_s,
                     const vector<double>& maps_x,
                     const vector<double>& maps_y) {
  int prev_wp = -1;
  while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
    prev_wp++;
  }
  int wp2 = (prev_wp + 1) % maps_x.size();
  double heading =
      atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s - maps_s[prev_wp]);
  double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
  double seg_y = maps_y[prev_wp] + seg_s * sin(heading);
  double perp_heading = heading - pi() / 2;
  double x = seg_x + d * cos(perp_heading);
  double y = seg_y + d * sin(perp_heading);
  return {x, y};
}
bool
double_equals(double d1, double d2, double tol = 1e-9) {
  double delta = d1 - d2;
  return fabs(delta) < tol;
}

double dist(double s1, double s2) {
  if (double_equals(s1, sentinel) || double_equals(s2, sentinel))
    return sentinel;
  const double thr = 1000;
  double ret = s2 - s1;
  if (ret > thr) {
    ret -= max_s;
  } else if (ret < -thr) {
    ret += max_s;
  }
  if (abs(ret) > thr) {
    std::cout << "s1 : " << s1 << " s2 : " << s2 << " ret : " << ret
              << " thr : " << thr << std::endl;
//    assert(abs(ret) < thr);
  }
  return ret;
}

double
rescind_s(double s, double ds) {
  s -= ds;
  if (s < 0) {
    s += max_s;
  }
  return s;
}

double
advance_s(double s, double ds) {
  s += ds;
  if (s > max_s) {
    s -= max_s;
  }
  return s;
}

double
advance_s(double dt, double s, double v, double a = 0.0, double adot = 0.0) {
  if (double_equals(s, sentinel))
    return sentinel;
  double ret = s + v * dt + 0.5 * a * dt * dt + adot * dt * dt * dt / 6;
  if (ret > max_s)
    ret -= max_s;
  return ret;
}

bool
possibleToChange(int laneId,
                 json& sensor_fusion,
                 double car_s /*already at ref_dt*/,
                 double car_v,
                 double ref_dt,
                 double& target_s,
                 double& target_vel) {

   double nearest_car_front_s = sentinel;
   double nearest_car_front_vel = 1000;
   double nearest_car_back_s = sentinel;
   double nearest_car_back_vel = 1000;
   for (auto sf:sensor_fusion) {
     double d = sf[6];
     if (d < lane_width * (laneId + 1) && d > lane_width * laneId) {
       double vx = sf[3];
       double vy = sf[4];
       double v = sqrt(vx * vx + vy * vy);
       double s = sf[5];
       s = advance_s(ref_dt, s, v);
       if (dist(car_s, s) > 0) {
         if (double_equals(nearest_car_front_s, sentinel)
             || dist(nearest_car_front_s, s) < 0) {
           nearest_car_front_s = s;
           nearest_car_front_vel = v;
         }
       } else if (dist(car_s, s) < 0) {
         if (double_equals(nearest_car_back_s, sentinel)
             || dist(nearest_car_back_s, s) > 0) {
           nearest_car_back_s = s;
           nearest_car_back_vel = v;
         }
       }
     }
   }
   // assuming 0 acceleration
   double nearest_car_back_s_time_horizon =
       advance_s(time_horizon, nearest_car_back_s, nearest_car_back_vel);
   double nearest_car_front_s_time_horizon =
       advance_s(time_horizon, nearest_car_front_s, nearest_car_front_vel);
   bool no_vehicle_behind =
       double_equals(nearest_car_back_s_time_horizon, sentinel);
   bool no_vehicle_front =
       double_equals(nearest_car_front_s_time_horizon, sentinel);
   if (no_vehicle_behind && no_vehicle_front) {
     target_s = advance_s(time_horizon, car_s, max_vel);
     target_vel = max_vel;
     return true;
   } else if (no_vehicle_behind) {
     if (dist(car_s, nearest_car_front_s) > 1.5 * lane_width) {
       target_s = rescind_s(nearest_car_front_s_time_horizon, 1.5 * lane_width);
       target_vel = std::min(max_vel, nearest_car_front_vel);
       return true;
     }
   } else if (no_vehicle_front) {
     if (dist(nearest_car_back_s, car_s) > 2 * lane_width) {
       target_s = advance_s(nearest_car_back_s_time_horizon, 3 * lane_width);
       target_vel = max_vel;
       return true;
     }
   } else {
     if ((dist(nearest_car_back_s_time_horizon, nearest_car_front_s_time_horizon)
         > 4 * lane_width) &&
         (dist(nearest_car_back_s, car_s) > 3 * lane_width) &&
         (dist(car_s, nearest_car_front_s) > 1.5 * lane_width)) {
       target_s = rescind_s(nearest_car_front_s_time_horizon, 1.5 * lane_width);
       target_vel = std::min(max_vel, nearest_car_front_vel);
       return true;
     }
   }
   return false;
 }

 enum DIRECTION { RIGHT,LEFT,SAME,NONE};

 std::tuple<std::function<double(double)>,double,double,double,DIRECTION> jerkFuncForSafeFollowingDist(double a0,double v0,double s0,std::vector<std::tuple<double,double,DIRECTION>> carsAheadData) {
   double tmax = -10;
   // jfn,v1,delta_t,delta_d,Direction,
   std::vector<std::tuple<std::function<double(double)>,double,double,double,DIRECTION>> carsAheadDataNew;
   for(auto carAheadData : carsAheadData) {
     double v1,s1;
     DIRECTION dir;
     std::tie(v1,s1,dir) = carAheadData;
     double const safeTime = 1; //seconds
     double const minDist = lane_width; //meters
     double delta_d = s1-s0-minDist-v1*safeTime;
     std::function<double(double)> jfn;
     double delta_t;
     std::tie(jfn,delta_t) = achieveZeroAcelAndVel(a0,v0-v1,-v1,max_vel-v1,delta_d);
     carsAheadDataNew.push_back(std::make_tuple(jfn,v1,delta_t,delta_d,dir));
     if(tmax<delta_t) {
       tmax = delta_t;
     }
   }
   auto h = [tmax](const std::tuple<std::function<double(double)>,double,double,double,DIRECTION>& x) {
              double v1=std::get<1>(x);
              double delta_d = std::get<3>(x);
              return delta_d+v1*tmax;
            };
   auto cmp = [h,tmax](const std::tuple<std::function<double(double)>,double,double,double,DIRECTION>& a,const std::tuple<std::function<double(double)>,double,double,double,DIRECTION>& b) {
                  return h(a)<h(b);
                };
   return *std::max_element(carsAheadDataNew.begin(),carsAheadDataNew.end(),cmp);
 }



vector<double> map_waypoints_x, map_waypoints_y, map_waypoints_s,map_waypoints_dx, map_waypoints_dy;

void read_map() {
   string map_file_ = "../data/highway_map.csv";
   ifstream in_map_(map_file_.c_str(), ifstream::in);
   string line;
   while (getline(in_map_, line)) {
     istringstream iss(line);
     double x, y;
     float s, d_x, d_y;
     iss >> x >> y >> s >> d_x >> d_y;
     map_waypoints_x.push_back(x);
     map_waypoints_y.push_back(y);
     map_waypoints_s.push_back(s);
     map_waypoints_dx.push_back(d_x);
     map_waypoints_dy.push_back(d_y);
   }
}

 int main() {
   read_map();
   uWS::Hub h;
   h.onMessage([](
       uWS::WebSocket<uWS::SERVER> ws,
       char* data,
       size_t length,
       uWS::OpCode opCode) {
     if (length && length > 2 && data[0] == '4' && data[1] == '2') {
       auto s = hasData(data);
       if (s != "") {
         auto j = json::parse(s);
         string event = j[0].get<string>();
         if (event == "telemetry") {
           double car_x(j[1]["x"]), car_y(j[1]["y"]), car_s(j[1]["s"]), car_d(j[1]["d"]), car_yaw(j[1]["yaw"]),
               car_speed(j[1]["speed"]), end_path_s(j[1]["end_path_s"]), end_path_d = j[1]["end_path_d"];
           car_yaw *= 0.03490658503;
           auto previous_path_x(j[1]["previous_path_x"]), previous_path_y(j[1]["previous_path_y"]);
           auto sensor_fusion = j[1]["sensor_fusion"];
           double car_lane;
           modf(car_d / lane_width, &car_lane);
           double ref_x(car_x), ref_y(car_y), ref_yaw(car_yaw);
           auto prev_size = previous_path_x.size();
          double end_path_lane;
          if (prev_size > 0)
            modf(end_path_d / lane_width, &end_path_lane);
          else
            end_path_lane = car_lane;
          json msgJson;
          vector<double>
              next_x_vals(previous_path_x.begin(), previous_path_x.end()),
              next_y_vals(previous_path_y.begin(), previous_path_y.end());
          if (prev_size < 50) {
            double left_lane_id = end_path_lane - 1;
            double right_lane_id = end_path_lane + 1;
            std::vector<double> laneIds;
            for (auto sf:sensor_fusion) {
              if (sf[6] < 0)
                continue;
              double curLaneId;
              modf(double(sf[6]) / lane_width, &curLaneId);
              laneIds.push_back(curLaneId);
            }
            if (prev_size > 0) {
              car_s = end_path_s;
              car_d = end_path_d;
            }
            double target_vel = max_vel;
            bool left_lane_change_possible(false);
            bool right_lane_change_possible(false);
            bool current_lane_possible(false);
            double left_s, right_s, same_s;
            double left_v, right_v, same_v;
            std::vector<std::tuple<double,double,DIRECTION>> carsAheadData;
            if (left_lane_id > -1) {
              left_lane_change_possible = possibleToChange(left_lane_id,
                                                           sensor_fusion,
                                                           car_s,
                                                           ref_vel_mps,
                                                           prev_size * dt,
                                                           left_s,
                                                           left_v);
              if(left_lane_change_possible) {
                carsAheadData.push_back(std::make_tuple(left_v,left_s,LEFT));
              }
            }
            if (right_lane_id < 3) {
              right_lane_change_possible = possibleToChange(right_lane_id,
                                                            sensor_fusion,
                                                            car_s,
                                                            ref_vel_mps,
                                                            prev_size * dt,
                                                            right_s,
                                                            right_v);
              if(right_lane_change_possible) {
                carsAheadData.push_back(std::make_tuple(right_v,right_s,RIGHT));
              }
            }
            current_lane_possible = possibleToChange(end_path_lane,
                                                     sensor_fusion,
                                                     car_s,
                                                     ref_vel_mps,
                                                     prev_size * dt,
                                                     same_s,
                                                     same_v);
            carsAheadData.push_back(std::make_tuple(same_v,same_s,SAME));
            std::function<double(double)> jfn;

            double delta_t;
            double delta_d;
            DIRECTION dir;
            std::tie(jfn,target_vel,delta_t,delta_d,dir) = jerkFuncForSafeFollowingDist(ref_acc_mps2,ref_vel_mps,car_s,carsAheadData);

            int best_lane_id;
            switch(dir) {
            case LEFT : if(!left_lane_change_possible) throw("left_lane_change_not possible");
              best_lane_id = left_lane_id; break;
            case RIGHT : if(!right_lane_change_possible) throw("right_lane_change_not possible");
              best_lane_id = right_lane_id; break;
            case SAME: best_lane_id = end_path_lane; break;
            default : throw("error");
            }

            {
              for(int i=prev_size;i<num_waypoints;i++) {
              auto nxt_pt = getXY(0,0,map_waypoints_s, map_waypoints_x, map_waypoints_y);
              }
            }
          }
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;
          auto msg = "42[\"control\"," + msgJson.dump() + "]";
          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });
  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse* res, uWS::HttpRequest req, char* data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });
  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });
  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char* message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });
  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
