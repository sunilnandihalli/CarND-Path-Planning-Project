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
const double max_s = 6945.554;
const double mps2mph = 2.23694;
const double dt = 0.02; // time between waypoints
const double time_horizon = 3.0; // time horizon in seconds
const int num_waypoints = time_horizon / dt;
const double lane_width = 4.0;
const double max_jerk = 5; // meters per sec^3
const double max_acc = 5; // meters per sec^2
const double max_vel_mph = 45;
const double max_vel = max_vel_mph / mps2mph; // meters per sec
const double max_acc_change_per_frame = max_jerk * dt;
const double sentinel = 999999999;
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
  /*
  std::cout << " lane_id : " << laneId
            << " nearest_car_back_s : " << nearest_car_back_s << " v : "
            << nearest_car_back_vel
            << " car_s              : " << car_s << " v : " << car_v
            << " nearest_car_front_s : " << nearest_car_front_s << " v : "
            << nearest_car_front_vel << std::endl;
  */
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
double calc_dist(double x1,double y1,double x2,double y2) {
  double dx=x2-x1,dy=y2-y1;
  return sqrt(dx*dx+dy*dy);
}
double calc_vel(double x1,double y1,double x2,double y2) {
  return calc_dist(x1,y1,x2,y2)/dt;
}
double calc_angle(double x1,double y1,double x2 ,double y2 ,double x3,double y3){
  double theta1 = atan2(x2-x1,y2-y1);
  double theta2 = atan2(x3-x2,y3-y2);
  return theta2-theta1;
}
double calc_omega(double x1,double y1,double x2 ,double y2 ,double x3,double y3){
  return calc_angle(x1,y1,x2,y2,x3,y3)/dt;
}
double calc_radius(double x1,double y1,double x2,double y2,double x3,double y3) {
  return calc_dist(x2,y2,x3,y3)/calc_angle(x1,y1,x2,y2,x3,y3);
}
double calc_linear_acc(double x1,double y1,double x2,double y2,double x3,double y3) {
  double v1 = calc_vel(x1,y1,x2,y2);
  double v2 = calc_vel(x2,y2,x3,y3);
  return (v2-v1)/dt;
}
double calc_radial_acc(double x1,double y1,double x2,double y2,double x3,double y3) {
  double v = calc_vel(x2,y2,x3,y3);
  double r = calc_radius(x1,y1,x2,y2,x3,y3);
  return v*v/r;
}
double calc_acc(double x1,double y1,double x2,double y2,double x3,double y3) {
  double al = calc_linear_acc(x1,y1,x2,y2,x3,y3);
  double ar = calc_radial_acc(x1,y1,x2,y2,x3,y3);
  return sqrt(al*al+ar*ar);
}
double calc_jerk(double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4) {
  double a1 = calc_acc(x1,y1,x2,y2,x3,y3);
  double a2 = calc_acc(x2,y2,x3,y3,x4,y4);
  return (a2-a1)/dt;
}
struct l4points {
  int n;
  double x1,x2,x3,x4,y1,y2,y3,y4;
  l4points() {
    n = 0;
    x1=x2=x3=x4=y1=y2=y3=y4=sentinel;
  }
  void addPoint(double x,double y) {
    x1=x2;
    x2=x3;
    x3=x4;
    x4=x;
    y1=y2;
    y2=y3;
    y3=y4;
    y4=y;
    n+=1;
  }

};
double jerk_s(double t, double k, double T) {
  if (t < 0.5 * T)
    return k;
  else if (t <= T)
    return -k;
}

double acc_s(double t, double k, double T) {
  if (t < 0.5 * T)
    return k * t;
  else if (t <= T)
    return k * (T - t);
}

double vel_s(double t, double k, double T, double v0) {
  if (t < 0.5 * T)
    return v0 + k * t * t * 0.5;
  else if (t <= T)
    return v0 + k * (T * t - t * t * 0.5 - T * T * 0.25);
}

double jerk_d(double t, double k, double T) {
  if (t < 0.25 * T)
    return k;
  else if (t < 0.75 * T)
    return -k;
  else if (t < T)
    return k;
}

double acc_d(double t, double k, double T) {
  if (t < 0.25 * T)
    return k * t;
  else if (t < 0.75 * T)
    return k * (0.5 * T - t);
  else if (t < T)
    return k * (t - T);
}

double vel_d(double t, double k, double T) {
  if (t < 0.25 * T)
    return k * t * t * 0.5;
  else if (t < 0.75 * T)
    return k * 0.5 * (T * t - t * t - T * T / 16);
  else if (t < T)
    return k * (t * t * 0.5 + T * T * 0.5 - t * T);
}

double calc_cur_delta_d(double t, double k, double T) {
  if (t < 0.25 * T)
    return k * t * t * t / 6;
  else if (t < 0.75 * T)
    return k * (T * t * t / 4 - t * t * t / 6 - T * T * t / 16
        + T * T * T / (3 * 64));
  else if (t < T)
    return k * (-T * t * t / 2 + t * t * t / 6 + T * T * t / 2
        - 26 * T * T * T / (3 * 64));
}

double calc_kd(double delta_d, double Td) {
  return 32 * delta_d / (Td * Td * Td);
}
bool should_drive_acc_to_zero(double acc,double v,
                       double j/*jerk*/,double delta_v /* final-current*/
    ,double delta_s) {

  if(acc*delta_v<0)
    return false;
  j = (acc>0?-1:1)*fabs(j);
  double delta_t = -acc/j;
  double smallest_delta_v = j * delta_t * delta_t * 0.5 + acc * delta_t;
  double smallest_delta_s =
      j * delta_t * delta_t * delta_t / 6 + acc * delta_t * delta_t * 0.5
          + v * delta_t;
  if(delta_v>0) {
    return delta_v < smallest_delta_v || delta_s < smallest_delta_s;
  } else if (delta_v<0) {
    return delta_v > smallest_delta_v ;
  }
}
bool calc_params(double delta_d,
                 double delta_sdot,
                 double& kd,
                 double& ks,
                 double Td,
                 double Ts) {
  kd = 32 * delta_d / (Td * Td * Td);
  ks = 4 * delta_sdot / (Ts * Ts);
  double cur_max_jerk = sqrt(kd * kd + ks * ks);
  double cur_max_acc_d = kd * Td / 4;
  double cur_max_acc_s = ks * Ts / 2;
  double cur_max_acc =
      sqrt(cur_max_acc_d * cur_max_acc_d + cur_max_acc_s * cur_max_acc_s);
  std::cout << " kd : " << kd << " ks : " << ks << " cur_max_acc_d : "
            << cur_max_acc_d << " cur_max_acc_s " << cur_max_acc_s
            << " cur_max_acc : " << cur_max_acc << " Ts : " << Ts << " Td : "
            << Td << " delta_d : " << delta_d << " delta_sdot : " << std::endl;
  return cur_max_jerk < max_jerk && cur_max_acc < max_acc;

}


int main() {
  uWS::Hub h;
  // Load up map values for waypoint's x,y,s and d normalized normal vectors

  int target_lane = 1;
  double ref_vel_mps = 0.0;
  double ref_acc_mps2 = 0.0;
  vector<double> map_waypoints_x, map_waypoints_y, map_waypoints_s,
      map_waypoints_dx, map_waypoints_dy;
  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0

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
  std::ofstream fout_wp("way_points.csv");
  std::ofstream fout_calls("calls.csv");
  std::ofstream fout_sf("sensor_fusion.csv");
  fout_wp<<"wp_id,x,y,call_id,radius,theta,v,linear_acc,radial_acc,total_acc,j"<<std::endl;
  fout_calls<<"call_id,prev_path_size,next_path_size,car_x,car_y,car_s,car_d,car_yaw,car_speed,end_path_s,end_path_d,end_path_x,end_path_y,num_cars_visible"<<std::endl;
  fout_sf<<"call_id,car_id,x,y,vx,vy,s,d"<<std::endl;
  int wp_id = 0;
  int call_id = 0;
  l4points l4p;
  h.onMessage([&l4p,&fout_wp,&fout_calls,&fout_sf,&wp_id,&call_id,&ref_vel_mps, &ref_acc_mps2, &map_waypoints_x,
               &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](
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
          // j[1] is the data JSON object
          // Main car's localization Data
          double car_x(j[1]["x"]), car_y(j[1]["y"]), car_s(j[1]["s"]),
              car_d(j[1]["d"]), car_yaw(j[1]["yaw"]),
              car_speed(j[1]["speed"]), end_path_s(j[1]["end_path_s"]),
              end_path_d = j[1]["end_path_d"];
          auto previous_path_x(j[1]["previous_path_x"]), previous_path_y
              (j[1]["previous_path_y"]);// Previous path data given to the Planner
          auto sensor_fusion = j[1]["sensor_fusion"];
          double car_lane;
          modf(car_d / lane_width, &car_lane);

          double ref_x(car_x), ref_y(car_y), ref_yaw(car_yaw);
          vector<double> ptsx, ptsy;
          auto prev_size = previous_path_x.size();
          double end_path_lane;
          if (prev_size > 0)
            modf(end_path_d / lane_width, &end_path_lane);
          else
            end_path_lane = car_lane;
          /*
          std::cout << "prev_size : " << prev_size << " car_lane : " << car_lane
                    << " car_d : " << car_d << " end_path_d : " << end_path_d
                    << " end_path_lane : " << end_path_lane << std::endl;
          */
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
                break;
              double curLaneId;
              modf(double(sf[6]) / lane_width, &curLaneId);
              laneIds.push_back(curLaneId);
            }
            if (prev_size > 0)
              car_s = end_path_s;
            double target_vel = max_vel;
            bool left_lane_change_possible(false);
            bool right_lane_change_possible(false);
            bool current_lane_possible(false);
            double left_s, right_s, same_s;
            double left_v, right_v, same_v;
            if (left_lane_id > -1) {
              left_lane_change_possible = possibleToChange(left_lane_id,
                                                           sensor_fusion,
                                                           car_s,
                                                           ref_vel_mps,
                                                           prev_size * dt,
                                                           left_s,
                                                           left_v);
            }
            if (right_lane_id < 3) {
              right_lane_change_possible = possibleToChange(right_lane_id,
                                                            sensor_fusion,
                                                            car_s,
                                                            ref_vel_mps,
                                                            prev_size * dt,
                                                            right_s,
                                                            right_v);
            }
            current_lane_possible = possibleToChange(end_path_lane,
                                                     sensor_fusion,
                                                     car_s,
                                                     ref_vel_mps,
                                                     prev_size * dt,
                                                     same_s,
                                                     same_v);
            double best_lane_id = end_path_lane;
            double best_speed = same_v;
            double best_s = same_s;
            if (left_lane_change_possible && left_v > best_speed) {
              best_lane_id = left_lane_id;
              best_speed = left_v;
              best_s = left_s;
            }
            if (right_lane_change_possible && right_v > best_speed) {
              best_lane_id = right_lane_id;
              best_speed = right_v;
              best_s = right_s;
            }
            /*
            std::cout << "left : " << left_lane_id << " left_s : " << left_s
                      << " left_v : " << left_v << std::endl;
            std::cout << "right : " << right_lane_id << " right_s : " << right_s
                      << " right_v : " << right_v << std::endl;
            std::cout << " best_lane_id : " << best_lane_id
                      << " end_path_lane : "
                      << end_path_lane << " best_speed : " << best_speed
                      << std::endl;
            std::cout << " left_possible : " << left_lane_change_possible
                      << " right_possible : " << right_lane_change_possible
                      << std::endl;
            */

            if (true) {
              if (prev_size < 2) {
                double prev_car_x = car_x - cos(car_yaw);
                double prev_car_y = car_y - sin(car_yaw);
                ptsx.insert(ptsx.end(), {prev_car_x, car_x});
                ptsy.insert(ptsy.end(), {prev_car_y, car_y});
              } else {
                ref_x = previous_path_x[prev_size - 1];
                ref_y = previous_path_y[prev_size - 1];
                double ref_x_prev = previous_path_x[prev_size - 2];
                double ref_y_prev = previous_path_y[prev_size - 2];
                ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);
                ptsx.insert(ptsx.end(), {ref_x_prev, ref_x});
                ptsy.insert(ptsy.end(), {ref_y_prev, ref_y});
              }


              double final_d = lane_width * (best_lane_id + 0.5);
              double final_s =
                  car_s + (ref_vel_mps + best_speed) * 0.5 * (time_horizon-dt*prev_size);

              for (double alpha:{0.9, 1.0}) {
                auto nxt_pt = getXY(car_s * (1 - alpha) + final_s * alpha,
                                    car_d * (1 - alpha) + final_d * alpha,
                                    map_waypoints_s,
                                    map_waypoints_x,
                                    map_waypoints_y);
                ptsx.push_back(nxt_pt[0]);
                ptsy.push_back(nxt_pt[1]);
              }
              for (int i = 0; i < ptsx.size(); i++) {
                double shift_x = ptsx[i] - ref_x;
                double shift_y = ptsy[i] - ref_y;
                ptsx[i] =
                    (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
                ptsy[i] =
                    (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
              }
              tk::spline s;
              s.set_points(ptsx, ptsy);
              double target_x = final_s;
              double target_y = s(target_x);
              double x_add_on = 0;
              double time_left = time_horizon-dt*prev_size;
              target_vel = best_speed;
              for (int i = 1; i <= num_waypoints - previous_path_x.size();
                   i++) {
                double desired_acc = (target_vel - ref_vel_mps) / time_left;
                if (desired_acc > ref_acc_mps2 + max_acc_change_per_frame)
                  desired_acc = ref_acc_mps2 + max_acc_change_per_frame;
                else if (desired_acc < ref_acc_mps2 - max_acc_change_per_frame)
                  desired_acc = ref_acc_mps2 - max_acc_change_per_frame;
                if (desired_acc > max_acc)
                  desired_acc = max_acc;
                else if (desired_acc < -max_acc)
                  desired_acc = -max_acc;
                ref_acc_mps2 = desired_acc;
                ref_vel_mps += ref_acc_mps2 * dt;
                time_left -= dt;
                double x_point = x_add_on + ref_vel_mps * dt;
                double y_point = s(x_point);
                x_add_on = x_point;
                double x_ref = x_point;
                double y_ref = y_point;
                x_point = x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw);
                y_point = x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw);
                x_point += ref_x;
                y_point += ref_y;
                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
              }
            } else if (false) {
              double delta_d = lane_width * (best_lane_id - end_path_lane);
              double delta_sdot = best_speed - ref_vel_mps;
              double Td = time_horizon;
              double kd;
              double ks;
              double Ts = time_horizon;
              double v0 = ref_vel_mps;
              double cur_s = end_path_s;
              double prev_sdot = v0;
              while (!calc_params(delta_d, delta_sdot, kd, ks, Td, Ts))
                Ts *= 1.1;
              std::cout << " best_speed : " << best_speed << " ref_vel_mps : "
                        << ref_vel_mps << std::endl;
              std::cout << " delta_d : " << delta_d << " delta_sdot : "
                        << delta_sdot << " kd : " << kd << " ks : " << ks
                        << " Td : " << Td << " Ts : " << Ts << std::endl;
              for (double cur_dt = dt; cur_dt <= Td; cur_dt += dt) {
                double cur_delta_d = calc_cur_delta_d(cur_dt, kd, Td);
                double cur_sdot = vel_s(cur_dt, ks, Ts, v0);
                cur_s += (prev_sdot + cur_sdot) * 0.5 * dt;
                std::cout << " s : " << cur_s << " d : " << cur_delta_d
                          << " cur_dt : " << cur_dt << " cur_sdot : "
                          << cur_sdot << std::endl;
                auto next_pt =
                    getXY(cur_s, end_path_d + cur_delta_d, map_waypoints_s,
                          map_waypoints_x, map_waypoints_y);
                next_x_vals.push_back(next_pt[0]);
                next_y_vals.push_back(next_pt[1]);
                ref_vel_mps = cur_sdot;
              }
            } else if (false) {
              double delta_d = lane_width * (best_lane_id - end_path_lane);
              double Td = time_horizon;
              double kd = calc_kd(delta_d, Td);
              double ks = sqrt(max_jerk * max_jerk - kd * kd);
              double max_acc_change = ks * dt;

              for (double cur_dt = dt; cur_dt <= Td; cur_dt += dt) {
                double cur_delta_d = calc_cur_delta_d(cur_dt, kd, Td);
                double delta_sdot = best_speed - ref_vel_mps;
                bool reduce_abs_acc = should_drive_acc_to_zero(ref_acc_mps2,delta_sdot,ks,delta_sdot,best_s-car_s);
                if (reduce_abs_acc) {
                  if(ref_acc_mps2>0) {
                    ref_acc_mps2=max(0.0,ref_acc_mps2-max_acc_change);
                  } else if (ref_acc_mps2<0) {
                    ref_acc_mps2 = min(0.0,ref_acc_mps2+max_acc_change);
                  }
                } else if(delta_sdot > 0) {
                  ref_acc_mps2 = std::min(ref_acc_mps2+max_acc_change,max_acc);
                } else if(delta_sdot <0) {
                  ref_acc_mps2 = std::max(ref_acc_mps2-max_acc_change,-max_acc);
                }
                ref_vel_mps = std::min(max_vel,ref_vel_mps+ref_acc_mps2 * dt);
                car_s = advance_s(dt,car_s,ref_vel_mps);
                auto next_pt = getXY(car_s,
                                     car_d + cur_delta_d,
                                     map_waypoints_s,
                                     map_waypoints_x,
                                     map_waypoints_y);
                next_x_vals.push_back(next_pt[0]);
                next_y_vals.push_back(next_pt[1]);
              }
            }
          }
          {
            for(int nxt_id =  prev_size;nxt_id<next_x_vals.size();nxt_id++) {
              l4p.addPoint(next_x_vals[nxt_id],next_y_vals[nxt_id]);
              double a(sentinel),j(sentinel),v(sentinel),radius(sentinel),theta(sentinel),
                  linear_acc(sentinel),radial_acc(sentinel);
              if(l4p.n>=4){
                radius = calc_radius(l4p.x2,l4p.y2,l4p.x3,l4p.y3,l4p.x4,l4p.y4);
                theta = calc_angle(l4p.x2,l4p.y2,l4p.x3,l4p.y3,l4p.x4,l4p.y4);
                linear_acc = calc_linear_acc(l4p.x2,l4p.y2,l4p.x3,l4p.y3,l4p.x4,l4p.y4);
                radial_acc = calc_radial_acc(l4p.x2,l4p.y2,l4p.x3,l4p.y3,l4p.x4,l4p.y4);
                v = calc_vel(l4p.x3,l4p.y3,l4p.x4,l4p.y4);
                a = calc_acc(l4p.x2,l4p.y2,l4p.x3,l4p.y3,l4p.x4,l4p.y4);
                j = calc_jerk(l4p.x1,l4p.y1,l4p.x2,l4p.y2,l4p.x3,l4p.y3,l4p.x4,l4p.y4);
              }
              fout_wp<<wp_id<<","<<next_x_vals[nxt_id]<<","<<next_y_vals[nxt_id]<<","
                     <<call_id<<","<<radius<<","<<theta<<","<<v<<","<<linear_acc<<","
                     <<radial_acc <<","<<a<<","<<j<<std::endl;
              wp_id+=1;
            }
            fout_wp<<std::flush;
            int num_cars_visible = 0;
            for(auto sf : sensor_fusion) {
              double d(sf[6]);
              if(d>=0.0 && d<=lane_width*3) {
                fout_sf<<call_id<<","<<sf[0]<<","<<sf[1]<<","<<sf[2]<<","<<sf[3]<<","<<sf[4]<<","<<sf[5]<<","<<sf[6]<<std::endl;
                num_cars_visible+=1;
              }
            }
            auto end_pt = getXY(j[1]["end_path_s"],j[1]["end_path_d"],map_waypoints_s,map_waypoints_x,map_waypoints_y);
            fout_calls<<call_id<<","<<prev_size<<","<<next_x_vals.size()<<","<<j[1]["x"]<<","<<j[1]["y"]<<","<<j[1]["s"]<<","
                      <<j[1]["d"]<<","<<j[1]["yaw"]<<","<<j[1]["speed"]<<","<<j[1]["end_path_s"]<<","
                      <<j[1]["end_path_d"]<<","<<end_pt[0]<<","<<end_pt[1]<<","<<num_cars_visible<<std::endl<<std::flush;
            std::cout<<"call_id : "<<call_id<<" prev_size : "<<prev_size<<" next_size : "<<next_x_vals.size()<<std::endl;
            call_id+=1;
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
