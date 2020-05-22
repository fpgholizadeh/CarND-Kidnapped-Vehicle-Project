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
using std::normal_distribution;
std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  std::default_random_engine gen;
  num_particles = 100;  // Set the number of particles
  //normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std[0]);
  // normal distributions for y and theta
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++) {
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1;
    particles.push_back(particle);
    weights.push_back(particle.weight);
  }
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

     /// creating a noise with zero mean and given stadard deviations
   normal_distribution<double> dist_x(0, std_pos[0]);
   normal_distribution<double> dist_y(0, std_pos[1]);
   normal_distribution<double> dist_theta(0, std_pos[2]);
  
   for (unsigned int i=0; i < particles.size(); i++){
  	double x_p = particles[i].x;;
  	double y_p = particles[i].y; 
  	double theta_p = particles[i].theta;
   	if (fabs(yaw_rate) >= 0.00001)
   {
    x_p +=  (velocity/yaw_rate)*(sin(theta_p+yaw_rate*delta_t)-sin(theta_p));
    y_p += (velocity/yaw_rate)*(cos(theta_p)-cos(theta_p+yaw_rate*delta_t));
    theta_p += yaw_rate*delta_t;
   } 
   else
   {
    x_p += velocity*delta_t*cos(theta_p);
    y_p += velocity*delta_t*sin(theta_p);
     
   }

   /// Adding the noise to the measured data     
   particles[i].x = x_p +  dist_x(gen);
   particles[i].y = y_p + dist_y(gen);
   particles[i].theta = theta_p + dist_theta(gen);
   
   }
}


void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
   
   for (unsigned int i = 0; i < observations.size(); i++)
   { 
     /// initialize the nearest_dist and id to some unrelaestic values
     double nearest_dist = 1000000;
     int nearest_id = -1;
     
     for (unsigned int j = 0; j < predicted.size(); j++)
     {
     double x_dis = observations[i].x - predicted[j].x;
     double y_dis = observations[i].y - predicted[j].y;
     double distance = sqrt(x_dis*x_dis + y_dis*y_dis); // calculate the distnace
     if (distance < nearest_dist)
     {                     
       //update the nearest distance and id
       nearest_dist = distance;
       nearest_id = predicted[j].id;
     }
     
     }
     observations[i].id = nearest_id;  //// assign the nearest landmark to the obseravation
   }
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * Update the weights of each particle using a mult-variate Gaussian 
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

    double sigma_x = std_landmark[0];
    double sigma_y = std_landmark[1];
    for (unsigned int i = 0 ; i < particles.size(); i++)
    {
        /// for each particle, the predicted measurement is calculated with respect to each landmark
        vector<LandmarkObs> predicted_measurements;
        double x_ = particles[i].x;
        double y_ = particles[i].y;
        double theta_ = particles[i].theta;
       
       for (unsigned int k =0 ; k < map_landmarks.landmark_list.size(); k++)
      {
      	double ld_list_x = map_landmarks.landmark_list[k].x_f;
       	double ld_list_y = map_landmarks.landmark_list[k].y_f;
       	double ld_list_id = map_landmarks.landmark_list[k].id_i;
         
       	double x_dis = x_ - ld_list_x ;
       	double y_dis = y_ - ld_list_y ;
       	double predict_dis = sqrt(x_dis*x_dis + y_dis*y_dis);
         
       	if (predict_dis <= sensor_range)
       	{
          predicted_measurements.push_back(LandmarkObs{ld_list_id ,ld_list_x , ld_list_y});
       	}
      }
      
    /// for each particle, the observed measurement is transfered from the vehicle's coordinate to the map's coordinate system
    vector <LandmarkObs> observations_map_cord;
    for (unsigned int j = 0; j < observations.size(); j++)
    {
     	double x_map;
     	double y_map;
     	x_map = x_ + (cos(theta_) * observations[j].x) - (sin(theta_) * observations[j].y);
     	y_map = y_ + (sin(theta_) * observations[j].x) + (cos(theta_) * observations[j].y);
    
    	observations_map_cord.push_back(LandmarkObs{observations[j].id,x_map, y_map});
    } 
      
    ///// using nearest neighbor technique to find the closed landmark to the observed measurements 
    dataAssociation(predicted_measurements,observations_map_cord);
    
    //// calculating posterior probablities
    double updated_weight = 1.0;

    for (unsigned int j = 0; j < observations_map_cord.size(); j++)
    {
     double xmap_obs = observations_map_cord[j].x;
     double ymap_obs = observations_map_cord[j].y;
     double id_obs = observations_map_cord[j].id;

     double mu_x;
     double mu_y;
    for (unsigned int k =0; k < predicted_measurements.size(); k++)
    {
     if(id_obs == predicted_measurements[k].id)
     {
     mu_x = predicted_measurements[k].x;
     mu_y = predicted_measurements[k].y;
     }
    }
    
    double normalizer;
    normalizer = (1/(2 * M_PI * sigma_x * sigma_y));
    // calculating the exponent term
   	double exponent;
   	exponent = (pow(xmap_obs - mu_x, 2) / (2 * pow(sigma_x, 2)))
               + (pow(ymap_obs - mu_y, 2) / (2 * pow(sigma_y, 2)));
    
  // calculate weight by multiplying normalization and exponent terms
  	double weight_;
  	weight_ = normalizer * exp(-exponent);                                          
  	updated_weight *= weight_;
  	
    }
    particles[i].weight = updated_weight;  
  }
}

void ParticleFilter::resample() {
  /**
   * Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
 // Get weights 
  vector<double> weights;
  for(int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }
  // get max weight
  double max_weight = *max_element(weights.begin(), weights.end());

  // Creating distributions
  std::uniform_real_distribution<float> dist_double(0.0, max_weight);
  std::uniform_int_distribution<int> dist_int(0, num_particles - 1);

  // creat index
  int index = dist_int(gen);

  double beta = 0.0;

  // wheel ( refer to the udacity course)
  vector<Particle> resampled_particles;
  for(int i = 0; i < num_particles; i++) {
    beta += dist_double(gen) * 2.0;
    while( beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampled_particles.push_back(particles[index]);
  }

  particles = resampled_particles;
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