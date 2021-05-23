/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;

static default_random_engine rand_gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
 
    particles.clear();
    weights.clear();
    num_particles = 300;  // TODO: Set the number of particles
    normal_distribution<double> rand_n_x{x, std[0]}, rand_n_y{y, std[1]}, rand_n_theta{theta, std[2]};

    for (int i = 0; i < num_particles; ++i){
      double p_x{rand_n_x(rand_gen)}, p_y{rand_n_y(rand_gen)}, p_theta{rand_n_theta(rand_gen)}, p_weight{1.0};
      Particle p;
      p.id = i;
      p.x = p_x;
      p.y = p_y;
      p.theta = p_theta;
      p.weight = p_weight;
      
      particles.push_back(p);
    }
	
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  
  normal_distribution<double> rand_n_x{0, std_pos[0]}, rand_n_y{0, std_pos[1]}, rand_n_theta{0, std_pos[2]};
  const double eps = 0.0001, vel_d = velocity * delta_t, yaw_d = yaw_rate * delta_t;
  
  double p_x_mean, p_y_mean, p_theta_mean;
  
  for (auto &p : particles){
    
    if (fabs(yaw_rate) > eps ){
      const double vel_yaw = velocity / yaw_rate;
      p_x_mean = p.x + vel_yaw * (sin(p.theta + yaw_d) - sin(p.theta));
      p_y_mean = p.y + vel_yaw * (cos(p.theta) - cos(p.theta + yaw_d));
	  p_theta_mean = p.theta + yaw_d;
    }
    else {
      p_x_mean = p.x + vel_d * cos(p.theta);
      p_y_mean = p.y + vel_d * sin(p.theta);
      p_theta_mean = p.theta;   
    }
      
    p.x = p_x_mean + rand_n_x(rand_gen);
    p.y = p_y_mean + rand_n_y(rand_gen);
    p.theta = p_theta_mean + rand_n_theta(rand_gen);
  }
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
 
  
  
  double min_dist, dist_value;
  //int lmk_id;

  for (auto &obs : observations) {
    min_dist = numeric_limits<double>::max();
    
    for (const auto &lmk : predicted){
     dist_value = dist(obs.x, obs.y, lmk.x, lmk.y);
     if (dist_value < min_dist) {
       min_dist = dist_value;
       //lmk_id = lmk.id;
       obs.id = lmk.id;
     }
    }
  }
  
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  
  // Define consts
  const double norm_2d_scaling_factor = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
  const double var_sens_x = pow(std_landmark[0], -2.0), var_sens_y = pow(std_landmark[1], -2.0); 
  
  for (auto &p : particles){
    
    // Find landmarks within sensor_range distance from particle p
    vector<LandmarkObs> predicted_obs;
    
    for (auto lmk : map_landmarks.landmark_list){
      if (dist(lmk.x_f, lmk.y_f, p.x, p.y) <  sensor_range){
        predicted_obs.push_back(LandmarkObs {lmk.id_i, lmk.x_f, lmk.y_f});
      }
    }
    
    // Transform observations to map coordinates
    vector<LandmarkObs> transformed_obs;
 	double x_m, y_m;
    
    for (auto obs : observations){
      x_m = cos(p.theta) * obs.x - sin(p.theta) * obs.y + p.x;
      y_m = sin(p.theta) * obs.x + cos(p.theta) * obs.y + p.y;
      
      transformed_obs.push_back(LandmarkObs {obs.id, x_m, y_m});  
      
    }
    
    // Data Association
    dataAssociation(predicted_obs, transformed_obs);
    
    // Update Weights
    // reset weight
    p.weight = 1.0;
    
    for (const auto &obs : transformed_obs){
      for (auto lmk : map_landmarks.landmark_list){
        if (obs.id == lmk.id_i){
          p.weight *= norm_2d_scaling_factor * exp(-0.5 * (pow(lmk.x_f - obs.x, 2.0) * var_sens_x +  pow(lmk.y_f - obs.y, 2.0) * var_sens_y));
        }   
      }     
    }   
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
  
  vector<double> weights;
  
  for (const auto &p : particles){
    weights.push_back(p.weight);
  }
  
  discrete_distribution<int> rand_d(weights.begin(), weights.end());
  
  vector<Particle> posterior_particles;
  
  for (int i = 0; i < num_particles; ++i){
    int new_idx = rand_d(rand_gen);
    posterior_particles.push_back(particles[new_idx]);
  }
  
  particles = posterior_particles; 
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}