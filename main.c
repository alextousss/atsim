#include <cglm/cglm.h>   /* for inline */
#include <cglm/call.h>   /* for library call (this also includes cglm.h) */
#include <stdio.h>
#include "gnuplot_i.h"


double get_random() { return (double)rand() / (double)RAND_MAX; }
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
    vec3 rot_noised; // rotation
    vec3 rot_filtered; // rotation
    versor ori; // orientation
    vec4 motor_speed;
    //mat3 transform;
} QuadState;

void qs_create(QuadState *qs) {

    glm_quatv(qs->ori, 0.10f, (vec3) {0.70f, 0.70f, 0.0f});
    glm_vec3((vec3) {0.0f, 0.0f, 0.0f}, qs->rot);
    glm_vec3((vec3) {0.0f, 0.0f, 0.0f}, qs->rot_noised);
    glm_vec3((vec4) {0.0f, 0.0f, 0.0f, 0.0f}, qs->motor_speed);
}


void qc_create(QuadConstants *qc) {
    mat3 temp_inertia;
    glm_mat3_copy((mat3) {
            1200365.8235e-9, -10487.8454e-9, -44.8982e-9,
            -10487.8454e-9, 1200365.8235e-9,-44.8982e-9,
            -44.8982e-9,  -44.8982e-9, 2391285.2518e-9
            },
            temp_inertia);
    glm_mat3_inv(temp_inertia, qc->inv_inertia_tensor);

}


QuadState qs_iterate(QuadState s, Forces f, QuadConstants c, unsigned int idx, double dt) {
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
    if(idx%10==0) {
        const float beta = exp(-2.0f*3.14f*20.0f*dt); // filtrage passe-bas 200hz
        for(int i = 0; i < 3 ; i++) {
            nqs.rot_noised[i] = nqs.rot[i] + 0.06f*(get_random()-0.5f);
            nqs.rot_filtered[i] = beta*s.rot_filtered[i] + (1.0f-beta)*nqs.rot_noised[i];
        }
    } else {
        for(int i = 0; i < 3 ; i++) {
            nqs.rot_noised[i] = s.rot_noised[i];
            nqs.rot_filtered[i] = s.rot_filtered[i];
        }
    }

    return nqs;
}


float X_K_P = 0.7f;
float X_K_D = 1.0f;
float K_P = 0.22f;
float K_D = 0.0012f;
#define LATENCY 15 // 30ms
void compute_commands(QuadState *nqs, QuadState *oqs, double dt, vec3 out) {
    for(int i = 0; i < 3; i++) {
        out[i] = 0;
    }
    vec3 euler_nqs;
    mat4 orientation;
    glm_quat_mat4(nqs->ori, orientation);
    glm_euler_angles(orientation, euler_nqs);

    vec3 euler_oqs;
    glm_quat_mat4(oqs->ori, orientation);
    glm_euler_angles(orientation, euler_oqs);



    vec3 ocsigne;
    vec3 ncsigne;




    for(unsigned int i = 0; i < 3; i++) {
        ncsigne[i] = X_K_P*-euler_nqs[i] + X_K_D*nqs->rot_filtered[i];
        ocsigne[i] = X_K_P*-euler_oqs[i] + X_K_D*oqs->rot_filtered[i];
        out[i] =   K_P * (ncsigne[i]-nqs->rot_filtered[i]) + K_D*((ncsigne[i]-nqs->rot_filtered[i])-(ocsigne[i]-oqs->rot_filtered[i]))/dt;
    }
}

void compute_moments(vec3 commands, vec3 moments, vec4 old_motor_value, double dt, vec4 new_motor_value) {
    for(int i = 0; i  < 3 ; i++)
        moments[i] = 0.0f;

    vec4 motor_value;
    //    for (unsigned int i = 0; i < 3; i++)
    //        commands[i] = sat(commands[i], -10, 10);
    float command_h = 0.35f;

    motor_value[0] = (1 * commands[1]) + (1 * commands[0]) + ( 1 * commands[2]) + command_h;
    motor_value[1] = ( -1 * commands[1]) + (1 * commands[0]) + (-1 * commands[2]) + command_h;
    motor_value[2] = ( -1 * commands[1]) + (-1 * commands[0]) + ( 1 * commands[2]) + command_h;
    motor_value[3] = (1 * commands[1]) + ( -1 * commands[0]) + (-1 * commands[2]) + command_h;



    vec4 motor_force;
    const float proj = 0.70710678118f*0.25;
    vec3 origin_to_motor[4][3] = {
        {-proj, proj, 0.0f},
        {proj, proj, 0.0f },
        {proj, -proj, 0.0f},
        {-proj, -proj, 0.0f},
    };
    vec4 temp;
    for(unsigned int i = 0; i < 4; i++) {
        motor_value[i] = sat(motor_value[i], 0.0f, 1.0f);
        const float beta = exp(-2.0f*3.14f*25.0f*dt);

        temp[i] = motor_value[i];
        motor_value[i] = beta*old_motor_value[i] + (1.0f-beta)*motor_value[i];
        old_motor_value[i] = temp[i];
        new_motor_value[i] = motor_value[i];

        motor_force[i] = motor_value[i]*0.4500f*9.81f; // (450g max)
//        printf("motor %d force is %f\n", i, motor_force[i]);
        vec3 force = {0.0f, 0.0f, 0.0f};
        force[2] = motor_force[i];
        vec3 moment;
        //glmc_vec3_print(force, stderr);
        //glmc_vec3_print(origin_to_motor[i], stderr);
        glm_vec3_cross(origin_to_motor[i], force, moment);
        for(unsigned int j = 0; j < 3; j++)
            moments[j] += moment[j];
        moments[0] += 10e-3f*9.81f*0.25f; // 10g à 25cm de perturbation (10e-3*9,81*0,25)
        //glmc_vec3_print(moment, stderr);
        // printf("---------\n");
    }

}


#define LENGTH 1
#define STEPS 5000
#define SIM_LENGTH  (LENGTH*STEPS)
#define N_SIMS 10e8
typedef struct SimInfo {
    float X_K_D;
    float X_K_P;
    float K_P;
    float K_D;
} SimInfo;
int main() {
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
    double Y5_s[SIM_LENGTH];
    double Y6_s[SIM_LENGTH];
    double Y7_s[SIM_LENGTH];
    qc_create(&c);
    float best_rms = 10000000.0f;
    SimInfo best_coeffs;
    for(unsigned long m_i = 0; m_i < N_SIMS; m_i++) {
        if(m_i%100 == 0) {
            X_K_D *= 1.01f;
        } else if (m_i%100000 == 0) {
            X_K_P *= 1.01f;
        } else if (m_i%10000000 == 0) {
            K_P *= 1.01f;
        } else if (m_i%1000000000 == 0) {
            K_D *= 1.01f;
        }

        for(unsigned int i =0; i <SIM_LENGTH;i++) {
            qs_create(&states[i]);

            X_s[i] = 0;
            Y1_s[i] = 0;
            Y2_s[i] = 0;
            Y3_s[i] = 0;
            Y4_s[i] = 0;
            Y5_s[i] = 0;
            Y6_s[i] = 0;
            Y7_s[i] = 0;
        }


        vec3 command;
        double rms = 0.0f;
        for(unsigned int i = 1000; i < SIM_LENGTH ; i++) {
            int s = i/STEPS;
          //  printf("%i\n", s);

            //glm_quat_normalize(states[i-1].ori);
            if(i%10==0) {
                compute_commands(&states[i-LATENCY*10], &states[i-(LATENCY*10+11)], dt, command);
            }

            vec3 moments;
            vec4 new_motor_speed;
            compute_moments(command, moments, states[i-1].motor_speed, dt, new_motor_speed);

            Forces f;

            glm_vec3_copy(moments, f.moment);

            states[i] = qs_iterate(states[i-1], f, c, i, dt);
            for(int j = 0; j < 4 ; j++)
                states[i].motor_speed[j] = new_motor_speed[j];
            float val = glm_quat_norm(states[i].ori);
            float angle = glm_quat_angle(states[i].ori);
            vec3 axis;
            glm_quat_axis(states[i].ori, axis);


            mat4 orientation;
            glm_quat_mat4(states[i].ori, orientation);
            vec3 euler;
            glm_euler_angles(orientation, euler);

            rms += euler[0]*euler[0];


            //printf("[Quat] Orientation (%f, %f, %f, %f) || ||=%f\n", states[i].ori[0], states[i].ori[1], states[i].ori[2], states[i].ori[3], val);

            //printf("[Vec3] Axe, angle (%f, %f, %f) angle=%f\n", axis[0], axis[1], axis[2], angle  );
            //printf("[Euler] Angles (%f, %f, %f)\n", euler[0], euler[1], euler[2]);
            //printf("[Vec3] Rotation (%f, %f, %f)\n", states[i].rot[0], states[i].rot[1], states[i].rot[2]  );



            /*   Y1_s[i] = -X_K_P*euler[0];
                 Y2_s[i] = states[i].rot_noised[0];
                 Y3_s[i] = states[i].rot_filtered[0]; //-X_K_P*euler[0];
                 Y4_s[i] = states[i].motor_speed[0];
                 Y5_s[i] = states[i-1].motor_speed[0];
                 Y6_s[i] = states[i].motor_speed[2];
                 Y7_s[i] = states[i-1].motor_speed[2];
                 X_s[i] = i*dt;
                 */
        }
        rms /= SIM_LENGTH;
        if(rms < best_rms) {
            best_rms = rms;
            best_coeffs = (SimInfo) {
                X_K_P,
                X_K_D,
                K_P,
                K_D
            };
        }
        printf("%d\t%lf\n", m_i, rms);

        /*
           gnuplot_ctrl * h ;
           h = gnuplot_init() ;
        //    gnuplot_cmd(h, "set multiplot layout 2,1 rowsfirst");
        gnuplot_plot_xy(h, X_s, Y1_s, SIM_LENGTH, "euler X") ;
        gnuplot_plot_xy(h, X_s, Y2_s, SIM_LENGTH, "rot noised X") ;
        gnuplot_plot_xy(h, X_s, Y3_s, SIM_LENGTH, "rot filtered X ") ;
        gnuplot_plot_xy(h, X_s, Y4_s, SIM_LENGTH, "motor 0") ;
        gnuplot_plot_xy(h, X_s, Y5_s, SIM_LENGTH, "motor 0 noised") ;
        gnuplot_plot_xy(h, X_s, Y6_s, SIM_LENGTH, "motor 2") ;
        gnuplot_plot_xy(h, X_s, Y7_s, SIM_LENGTH, "motor 2 noised") ;
        sleep(1000);

        gnuplot_close(h) ;
        */
    }
    printf("%f\t%f\t%f\t%f\n", best_coeffs.X_K_P, best_coeffs.X_K_D, best_coeffs.K_P, best_coeffs.K_D);
    return 0;
}
