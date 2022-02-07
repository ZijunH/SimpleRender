#ifndef CAMERA_HPP
#define CAMERA_HPP

#include "rt.hpp"

class camera {
   public:
    camera(point3 lookfrom, point3 lookat, vec3 vup, double vfov, double aspect_ratio, double aperture, double focus_dist, double o_time = 0.0, double c_time = 0.0) {
        const double theta = degrees_to_radians(vfov);
        const double h = tan(theta/2);
        const double viewport_height = 2.0 * h;
        const double viewport_width = aspect_ratio * viewport_height;

        w = unit(lookfrom - lookat);
        u = unit(cross(vup, w));
        v = cross(w, u);

        origin = lookfrom;
        horizontal = focus_dist * viewport_width * u;
        vertical = focus_dist * viewport_height * v;
        lower_left = origin - horizontal / 2 - vertical / 2 - focus_dist * w;

        lens_radius = aperture / 2;

        open_time = o_time;
        close_time = c_time;
    }

    ray get_ray(double s, double t) const {
        vec3 rd = lens_radius * random_in_unit_disk();
        vec3 offset = u * rd.x() + v * rd.y();

        return ray(origin + offset, lower_left + s * horizontal + t * vertical - origin - offset, random_double(open_time, close_time));
    }

   private:
    point3 origin;
    point3 lower_left;
    vec3 horizontal;
    vec3 vertical;
    vec3 u;
    vec3 v;
    vec3 w;
    double lens_radius;
    double open_time;
    double close_time;
};

#endif