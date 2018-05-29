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
const double time_horizon = 2; // time horizon in seconds
const int num_waypoints = time_horizon / dt;
const double lane_width = 4.0;
const double max_jerk = 5; // meters per sec^3
const double max_acc = 5; // meters per sec^2
const double max_vel_mph = 49;
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
  if(abs(ret)>thr) {
    std::cout<<"s1 : "<<s1<<" s2 : "<<s2<<" ret : "<<ret<<" thr : "<<thr<<std::endl;
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
  std::cout<<" lane_id : "<<laneId
           <<" nearest_car_back_s : "<<nearest_car_back_s<<" v : "<<nearest_car_back_vel
           <<" car_s              : "<<car_s <<" v : "<<car_v
           <<" nearest_car_front_s : "<<nearest_car_front_s<<" v : "<<nearest_car_front_vel<<std::endl;
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
double jerk_s(double t,double k,double T) {
  if(t<0.5*T)
    return k;
  else if(t<=T)
    return -k;
}

double acc_s(double t,double k,double T) {
  if(t<0.5*T)
    return k*t;
  else if (t<=T)
    return k*(T-t);
}

double vel_s(double t,double k,double T,double v0) {
  if(t<0.5*T)
    return v0+k*t*t*0.5;
  else if (t<=T)
    return v0+k*(T*t-t*t*0.5-T*T*0.25);
}

double jerk_d(double t,double k,double T) {
  if(t<0.25*T)
    return k;
  else if (t<0.75*T)
    return -k;
  else if (t<T)
    return k;
}

double acc_d(double t,double k,double T) {
  if(t<0.25*T)
    return k*t;
  else if(t<0.75*T)
    return k*(0.5*T-t);
  else if(t<T)
    return k*(t-T);
}

double vel_d(double t,double k,double T) {
  if(t<0.25*T)
    return k*t*t*0.5;
  else if(t<0.75*T)
    return k*0.5*(T*t-t*t-T*T/16);
  else if(t<T)
    return k*(t*t*0.5+T*T*0.5-t*T);
}

double delta_d(double t,double k,double T) {
  if(t<0.25*T)
    return k*t*t*t/6;
  else if(t<0.75*T)
    return k*(T*t*t/4 - t*t*t/6 - T*T*t/16 + T*T*T/(3*64));
  else if(t<T)
    return k*(-T*t*t/2+t*t*t/6+T*T*t/2-26*T*T*T/(3*64));
}

void calc_params(double delta_d,double delta_sdot,double& kd,double& ks,double& Td,double& Ts) {
  kd = 32*delta_d/(Td*Td*Td);
  ks = 4*delta_sdot/(Ts*Ts);
}

int main() {
  uWS::Hub h;
  // Load up map values for waypoint's x,y,s and d normalized normal vectors

  int target_lane = 1;
  double ref_vel_mps = 0.0;
  double ref_acc_mps2 = max_acc;
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
  h.onMessage([&ref_vel_mps, &ref_acc_mps2, &map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](
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
          std::cout << "prev_size : " << prev_size << " car_lane : " << car_lane
                    << " car_d : " << car_d << " end_path_d : " << end_path_d
                    << " end_path_lane : " << end_path_lane << std::endl;
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
          if(left_lane_change_possible && left_v>best_speed) {
            best_lane_id = left_lane_id;
            best_speed = left_v;
          }
          if(right_lane_change_possible && right_v>best_speed) {
            best_lane_id = right_lane_id;
            best_speed = right_v;
          }
          std::cout<<"left : "<<left_lane_id<<" left_s : "<<left_s<<" left_v : "<<left_v<<std::endl;
          std::cout<<"right : "<<right_lane_id<<" right_s : "<<right_s<<" right_v : "<<right_v<<std::endl;
          std::cout<<" best_lane_id : "<<best_lane_id<<" end_path_lane : "<<end_path_lane<<" best_speed : "<<best_speed<<std::endl;
          std::cout<<" left_possible : "<<left_lane_change_possible<<" right_possible : "<<right_lane_change_possible<<std::endl;
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


          double final_d = lane_width*(best_lane_id+0.5);
          double final_s = car_s + (ref_vel_mps+best_speed)*0.5*time_horizon;

          for(double alpha:{0.9,1.0}) {
            auto nxt_pt = getXY(car_s*(1-alpha)+final_s*alpha,
            car_d*(1-alpha)+final_d*alpha,
            map_waypoints_s,map_waypoints_x,map_waypoints_y);
            ptsx.push_back(nxt_pt[0]);
            ptsy.push_back(nxt_pt[1]);
          }
          for (int i = 0; i < ptsx.size(); i++) {
            double shift_x = ptsx[i] - ref_x;
            double shift_y = ptsy[i] - ref_y;
            ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
            ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
          }
          tk::spline s;
          s.set_points(ptsx, ptsy);
          json msgJson;
          vector<double>
              next_x_vals(previous_path_x.begin(), previous_path_x.end()),
              next_y_vals(previous_path_y.begin(), previous_path_y.end());
          double target_x = final_s;
          double target_y = s(target_x);
          double x_add_on = 0;
          double time_left = time_horizon;
          target_vel = best_speed;
          for (int i = 1; i <= num_waypoints - previous_path_x.size(); i++) {
            double desired_acc = (target_vel - ref_vel_mps) / time_left;
            if (desired_acc > ref_acc_mps2 + max_acc_change_per_frame)
              desired_acc = ref_acc_mps2 + max_acc_change_per_frame;
            else if (desired_acc < ref_acc_mps2 - max_acc_change_per_frame)
              desired_acc = ref_acc_mps2 - max_acc_change_per_frame;
            if(desired_acc>max_acc)
              desired_acc = max_acc;
            else if (desired_acc<-max_acc)
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
