#include <cglm/cglm.h>   /* for inline */
#include <cglm/call.h>   /* for library call (this also includes cglm.h) */
#include <stdio.h>
#include "gnuplot_i.h"


double max(double a, double b) {
    return a < b ? b : a;
}

double min(double a, double b) {
    return a < b ? a : b;
}

double sat(double a, double min, double max) {
    if (a > max) return max;
    if (a < min) return min;
    return a;
}
typedef struct QuadConstants {
    mat3 inv_inertia_tensor;
} QuadConstants;


typedef struct Forces {
    vec3 moment;
} Forces;


typedef struct QuadState {
    //vec3 pos; // position
    //vec3 vel; // velocity
    vec3 rot; // rotation
    versor ori; // orientation
    //mat3 transform;
} QuadState;

void qs_create(QuadState *qs) {

    mat4  rot;
    glm_euler((vec3){0.10f, 0.0f, 0.0f}, rot);
    glm_mat4_quat(rot, qs->ori);
    //    glm_quatv(qs->ori, 10, (vec3) {0.0f, 1.0f, 1.0f});
    //glm_vec3((vec3) {0.0f, 0.0f, 1.0f}, qs->pos);
    //glm_vec3((vec3) {0.0f, 0.0f, 0.0f}, qs->vel);
    glm_vec3((vec3) {0.0f, 0.0f, 0.0f}, qs->rot);
    glm_vec3((vec3) {0.0f, 0.0f, 0.0f}, qs->rot);
    glm_vec3((vec3) {0.0f, 0.0f, 0.0f}, qs->rot);
}


void qc_create(QuadConstants *qc) {
    mat3 temp_inertia;
    glm_mat3_copy((mat3) {
            1200365.8235*10e-9, -10487.8454*10e-9, -44.8982*10e-9,
            -10487.8454*10e-9, 1200365.8235*10e-9,-44.8982*10e-9,
            -44.8982*10e-9,  -44.8982*10e-9, 2391285.2518*10e-9
            },
            temp_inertia);
    glm_mat3_inv(temp_inertia, qc->inv_inertia_tensor);

}


QuadState qs_iterate(QuadState s, Forces f, QuadConstants c, double dt) {
    QuadState nqs;
    // Calculating derivative of rotation
    // theta'' = I^-1*Moments in the local basis
    vec3 d_dt_rotation_local_coords;
    glm_mat3_mulv(c.inv_inertia_tensor, f.moment, d_dt_rotation_local_coords );

    vec3 d_dt_rotation_world_coords;
    glm_quat_rotatev(s.ori, d_dt_rotation_local_coords, d_dt_rotation_world_coords);

    // Updating angular speed rot from forces
    // Assuming forces are in world coordinates
    for(unsigned int i = 0; i < 3 ; i++) {
        nqs.rot[i] = s.rot[i] + d_dt_rotation_world_coords[i]*dt;
    }


    // Updating orientation quaternion from angular speed rot
    versor omega;
    glm_quat_init(omega, nqs.rot[0], nqs.rot[1], nqs.rot[2], 0);
    versor omega_times_ori;
    glm_quat_mul(omega, s.ori, omega_times_ori);
    for(unsigned int i = 0; i < 4; i++) {
        nqs.ori[i] = s.ori[i] + dt/2*omega_times_ori[i];
    }
    return nqs;
}

#define K_P 0.4f
#define K_D 0.04f
void compute_commands(QuadState *nqs, QuadState *oqs, double dt, vec3 out) {
    mat4 orientation;
    vec3 euler_nqs;
    vec3 euler_oqs;
    glm_quat_mat4(nqs->ori, orientation);
    glm_euler_angles(orientation, euler_nqs);

    glm_quat_mat4(oqs->ori, orientation);
    glm_euler_angles(orientation, euler_oqs);
    for(unsigned int i = 0; i < 3; i++) {
        out[i] =- K_P * euler_nqs[i] + K_D*(euler_oqs[i]-euler_nqs[i])/dt;

    }
}

void compute_moments(vec3 commands, vec3 moments) {

    //    for (unsigned int i = 0; i < 3; i++)
    //        commands[i] = sat(commands[i], -10, 10);
    vec4 motor_value;
    float command_h = 0.0f;
    motor_value[1] = ( -1 * commands[1]) + (1 * commands[0]) + (-1 * commands[2]) + command_h;
    motor_value[2] = ( -1 * commands[1]) + (-1 * commands[0]) + ( 1 * commands[2]) + command_h;
    motor_value[3] = (1 * commands[1]) + ( -1 * commands[0]) + (-1 * commands[2]) + command_h;
    motor_value[0] = (1 * commands[1]) + (1 * commands[0]) + ( 1 * commands[2]) + command_h;
    vec4 motor_force;
    const float proj = 0.70710678118f*0.25;
    vec3 origin_to_motor[4][3] = {
        {-proj, proj, 0.0f},
        {proj, proj, 0.0f },
        {proj, -proj, 0.0f},
        {-proj, -proj, 0.0f},
    };
    for(unsigned int i = 0; i < 4; i++) {
        motor_value[i] = sat(motor_value[i], -0.5f, 0.5f);
        //printf("%f\n", motor_value[i]);
        motor_force[i] = motor_value[i]*0.4500f*9.81f; // (450g max)
        vec3 force = {0.0f, 0.0f, 0.0f};
        force[2] = motor_force[i];
        vec3 moment;
        //glmc_vec3_print(force, stderr);
        //glmc_vec3_print(origin_to_motor[i], stderr);
        glm_vec3_cross(origin_to_motor[i], force, moment);
        for(unsigned int j = 0; j < 3; j++)
            moments[j] += moment[j];

        //glmc_vec3_print(moment, stderr);
        // printf("---------\n");
    }

}


#define LENGTH 1
#define STEPS 10000
#define SIM_LENGTH  (LENGTH*STEPS)
int main() {
    vec3 moments;
    /*
       moments[0] = 0;
       moments[1] = 0;
       moments[2] = 0;
       compute_moments((vec3) {-1.0f, -0.0f, 0.0f}, moments);
       glmc_vec3_print(moments, stderr);
       */
    QuadConstants c; // 0,2÷12×(0,15²+0,04²)
    QuadState states[SIM_LENGTH];
    double dt = 1.0f/STEPS;

    double X_s[SIM_LENGTH];
    double Y1_s[SIM_LENGTH];
    double Y2_s[SIM_LENGTH];

    double Y3_s[SIM_LENGTH];
    double Y4_s[SIM_LENGTH];

    qs_create(&states[0]);
    qs_create(&states[1]);
    qc_create(&c);

    X_s[0] = 0;
    Y1_s[0] = 0;
    Y2_s[0] = 0;
    Y3_s[0] = 0;
    Y4_s[0] = 0;


    for(unsigned int i = 2; i < SIM_LENGTH ; i++) {
        int s = i/STEPS;
        printf("%i\n", s);

        //glm_quat_normalize(states[i-1].ori);
        vec3 command;
        compute_commands(&states[i-1], &states[i-2], dt, command);

        vec3 moments;
        compute_moments(command, moments);

        Forces f;

        glm_vec3_copy(moments, f.moment);

        states[i] = qs_iterate(states[i-1], f, c, dt);
        float val = glm_quat_norm(states[i].ori);
        float angle = glm_quat_angle(states[i].ori);
        vec3 axis;
        glm_quat_axis(states[i].ori, axis);


        mat4 orientation;
        glm_quat_mat4(states[i].ori, orientation);
        vec3 euler;
        glm_euler_angles(orientation, euler);

        printf("[Quat] Orientation (%f, %f, %f, %f) || ||=%f\n", states[i].ori[0], states[i].ori[1], states[i].ori[2], states[i].ori[3], val);

        printf("[Vec3] Axe, angle (%f, %f, %f) angle=%f\n", axis[0], axis[1], axis[2], angle  );
        printf("[Euler] Angles (%f, %f, %f)\n", euler[0], euler[1], euler[2]);
        printf("[Vec3] Rotation (%f, %f, %f)\n", states[i].rot[0], states[i].rot[1], states[i].rot[2]  );



        Y1_s[i] = euler[0]*500;
        Y2_s[i] = euler[1]*500;
        //Y3_s[i] = sat(command[0], -10, 10);
        Y3_s[i] = moments[0];
        Y4_s[i] = moments[1];
        X_s[i] = i*dt;
    }


    gnuplot_ctrl * h ;
    h = gnuplot_init() ;
    //    gnuplot_cmd(h, "set multiplot layout 2,1 rowsfirst");
    gnuplot_plot_xy(h, X_s, Y1_s, SIM_LENGTH, "euler X") ;
    gnuplot_plot_xy(h, X_s, Y2_s, SIM_LENGTH, "euler Y") ;
    gnuplot_plot_xy(h, X_s, Y3_s, SIM_LENGTH, "cmd X") ;
    gnuplot_plot_xy(h, X_s, Y4_s, SIM_LENGTH, "cmd Y") ;
    sleep(1000);

    gnuplot_close(h) ;
    return 0;
}
