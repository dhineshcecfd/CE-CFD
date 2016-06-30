#include <iostream>
#include <math.h>
#include <vector>
#include<fstream>
#include <stdlib.h>
#include<sys/time.h>
#include<omp.h>

void boundary(std::vector <double> &u, double N){
    double h = 1.0/(N-1);

    //bottom left
    double x1 = -1.0;
    double y1 = -1.0;
        for (int j1=0; j1<(N/2);j1++){
//            std::cout << x1 << std::endl;
//            std::cout << "J : " << j1 << std::endl;
            u[j1] = pow((x1*x1 + y1*y1), 0.25) * sin (0.5 * (atan(y1/x1)+M_PI));
            x1 += (2*h);
        }
        //bottom right
            for (int j1=(N/2)+1; j1<N;j1++){
//                std::cout << x1 << std::endl;
//                std::cout << "J : " << j1 << std::endl;
                u[j1] = pow ((x1*x1 + y1*y1), 0.25) * sin (0.5 * (atan(y1/x1)+(2*M_PI)));
                x1 += (2*h);
            }

        //right bottom
        double x2 = 1.0;
        double y2 = -1.0;
            for (int j2=1; j2<(N/2);j2++){
//                std::cout << y2 << std::endl;
//                std::cout << "J : " << ((j2)*N)-1 << std::endl;
                u[(j2*N)-1] = pow((x2*x2 + y2*y2),0.25) * sin (0.5 * (atan(y2/x2) +(2*M_PI)));
                y2 += 2*h;
            }
            //right top
                for (int j2=(N/2)+1; j2<=N;j2++){
//                    std::cout << y2 << std::endl;
//                    std::cout << "J : " << ((j2)*N)-1 << std::endl;
                    u[(j2*N)-1] = pow((x2*x2 + y2*y2), 0.25) * sin (0.5 * atan(y2/x2));
                    y2 += 2*h;
                }

            //left
            double x = -1.0;
            double y = -1.0;
                for (int j=1; j<N;j++){
                    y += 2*h;
//                    std::cout << y << std::endl;
//                    std::cout << "J : " << (j*N) << std::endl;
                    u[(j*N)] = pow ((x*x + y*y), 0.25) * sin (0.5 * (atan(y/x)+M_PI));

                }

                //top left
                double y3 = 1.0;
                double x3 = -1.0;
                    for (int j3=1; j3<(N/2)-1;j3++){
                        x3 += 2*h;
//                        std::cout << x3 << std::endl;
//                        std::cout << "J : " << (N*(N-1))+j3 << std::endl;
                        u[(N*(N-1))+j3] = pow ((x3*x3 + y3*y3), 0.25) * sin (0.5 * (atan(y3/x3)+M_PI));

                    }
                    //top right
                    for (int j3=N/2; j3<N-1;j3++){
                        x3 += 2*h;
//                        std::cout << x3 << std::endl;
//                        std::cout << "J : " << (N*(N-1))+j3 << std::endl;
                        u[(N*(N-1))+j3] = pow ((x3*x3 + y3*y3), 0.25) * sin (0.5 * atan(y3/x3));

                    }

                    //slit implementaion
                    double x0 = 0.0;
                    double y0 = 0.0;
                    for (int j0=1; j0<(N/2); j0++){
                        x0 += 2*h;
                        u[(N*N-1)/2+j0] = pow((x0*x0 + y0*y0),0.25) * sin (0.5 * atan(y0/x0));
                    }
}

void smoothing_rbgs(std::vector<double> &u, std::vector<double> &func, double l){

    int N = pow (2, l)+1;
    double N_h = pow (2.0, l);
    double h = 2.0 / N_h;
    double left_right = 1.0 / (h*h);
    double top_bottom = 1.0 / (h*h);
    double centre = 4.0 / (h*h);
    int temp, j;
#pragma omp parallel
{
#pragma omp for private(temp, j)
    //red
    for (int i=1; i<N-1; i++){
        temp = 1+(i%2);
        for (j=temp; j<N-1; j+=2){
            if ((i*N+j) >= (N*N-1)/2 && (i*N+j) <= ((N*N-1)/2) + (N/2)){
                //leave points
//                std::cout<<"Red : " << i*N+j <<std::endl;
            }else{
            u[i*N + j] = (func [i*N + j] +
                                      left_right * (u [i*N + j - 1] + u [i*N + j + 1]) +
                                            top_bottom * (u [(i - 1)*N + j] + u [(i + 1)*N + j])) / centre;
    }
        }
    }
#pragma omp barrier

#pragma omp for private(temp, j)
    //black
    for (int i = 1; i < N-1; i++) {
        temp = 2 - (i % 2);
            for (j = temp; j < N-1; j += 2) {
                if ((i*N+j) >= (N*N-1)/2 && (i*N+j) <= ((N*N-1)/2) + (N/2)){
                    //leave poits
//                    std::cout<<"Black : " << i*N+j <<std::endl;
                }else{
                u[i*N + j] = (func [i*N + j] +
                                          left_right * (u [i*N + j - 1] + u [i*N + j + 1]) +
                                                top_bottom * (u [(i - 1)*N + j] + u [(i + 1)*N + j])) / centre;

                }
            }
    }
#pragma omp barrier
}
}

void residual(std::vector<double>&u, std::vector<double>&func, std::vector<double>&res, double l){
    int N = pow (2, l)+1;
    double N_h = pow (2.0, l);
    double h = 2.0 / N_h;
    double left_right = 1.0 / (h*h);
    double top_bottom = 1.0 / (h*h);
    double centre = -4.0 / (h*h);
    int j;

#pragma omp parallel for private(j)
    //calculate residual
    for (int i=1; i<N-1; i++){
        for (j=1; j<N-1;j++){
            if((i*N+j) >= (N*N-1)/2 && (i*N+j) <= ((N*N-1)/2) + (N/2)){
                //leave points
            }else{
            res[i*N+j] = func[i*N+j] + (centre * u[i*N+j]) +
                                            (left_right * (u[i*N+j-1] + u[i*N+j+1])) +
                                                (top_bottom * (u[(i-1)*N+j] + u[(i+1)*N+j]));
            }
        }
    }
}

//double L2norm_v_cycle(std::vector<double>&u, std::vector<double>&func, std::vector<double>&res, double l){
//    int N = pow (2, l)+1;
//    double N_h = pow (2.0, l);
//    double h = 1.0 / N_h;
//    double left_right = 1.0 / (h*h);
//    double top_bottom = 1.0 / (h*h);
//    double centre = -4.0 / (h*h);
//    double temp = 0.0;
//    double L2norm = 0.0;

////#pragma omp parallel sections private(j) reduction(+:temp)
//    //calculate residual
//    for (int i=1; i<N-1; i++){
//        for (int j=1; j<N-1;j++){
//            if((i*N+j) >= (N*N-1)/2 && (i*N+j) <= ((N*N-1)/2) + (N/2)){
//                //leave points
//            }else{
//            res[i*N+j] = func[i*N+j] + (centre * u[i*N+j]) +
//                                            (left_right * (u[i*N+j-1] + u[i*N+j+1])) +
//                                                (top_bottom * (u[(i-1)*N+j] + u[(i+1)*N+j]));
//            temp += pow (res[i*N+j], 2);
//            }
//        }
//    }
//    L2norm = sqrt(temp/(N-2)*(N-2));
////    std::cout << "L2norm : "<< L2norm <<std::endl;
//    return L2norm;

//}

void restriction_residual(std::vector<double>&res, std::vector<double>&func, int l){
    int N = pow (2.0, l)+1;
    int New_N = (N/2)+1;
    int J;
#pragma omp parallel for schedule(static)
    for (int I=1; I<New_N-1; I++){
        int i = 2*I;
//        std::cout<< "I : " << I << "    "<< "i : " << i <<std::endl;
        for (J=1; J<New_N-1; J++){
            int j = 2*J;
//            std::cout<< "J : " << J << "j : " << j <<std::endl;
            func[I*New_N+J] = 0.0625 * (4*res[i*N+j] +  //center
                                            2*res[i*N+j-1] + 2*res[i*N+j+1] +   //left_right
                                                2*res[(i+1)*N+j] + 2*res[(i-1)*N+j] +   //top_bottom
                                                    res[(i+1)*N+j-1] + res[(i+1)*N+j+1] +   //topleft_topright
                                                        res[(i-1)*N+j-1] + res[(i-1)*N+j+1]);   //bottomleft_bottomright
//            std::cout << "I*New_N+J : " << I*New_N+J << "   " << "i*N+j : " << i*N+j << std::endl;
//            std::cout << "func : " << func[I*New_N+J] << std::endl;
        }
    }
}

void interpolate (std::vector<double>&u_old, std::vector<double>&u_new, int l){
    int N = pow (2.0, l)+1;
    int New_N = (2*N)-1;
    int j;

//#pragma omp parallel for schedule(static)
    for(int i=1; i<N-1; i++){
        int I = 2*i;
        for (j=1; j<N-1; j++){
            int J = 2*j;
//            std::cout<<"old"<<"\t"<<i*N+j<<"\t"<<"new"<<"\t"<<I*New_N+J<< std::endl;
//            std::cout<<"old"<<"\t"<<i*N+j<<"\t"<<"new"<<"\t"<<I*New_N+J-1<< std::endl;
//            std::cout<<"old"<<"\t"<<i*N+j<<"\t"<<"new"<<"\t"<<I*New_N+J+1<< std::endl;
//            std::cout<<"old"<<"\t"<<i*N+j<<"\t"<<"new"<<"\t"<<(I-1)*New_N+J<< std::endl;
//            std::cout<<"old"<<"\t"<<i*N+j<<"\t"<<"new"<<"\t"<<(I+1)*New_N+J<< std::endl;
//            std::cout<<"old"<<"\t"<<i*N+j<<"\t"<<"new"<<"\t"<<(I-1)*New_N+(J-1)<< std::endl;
//            std::cout<<"old"<<"\t"<<i*N+j<<"\t"<<"new"<<"\t"<<(I-1)*New_N+(J+1)<< std::endl;
//            std::cout<<"old"<<"\t"<<i*N+j<<"\t"<<"new"<<"\t"<<(I+1)*New_N+(J-1)<< std::endl;
//            std::cout<<"old"<<"\t"<<i*N+j<<"\t"<<"new"<<"\t"<<(I+1)*New_N+(J+1)<< std::endl;
            u_new[I*New_N+J] += u_old[i*N+j]; //centre
            u_new[I*New_N+J-1] += u_old[i*N+j]/2.0; //left
            u_new[I*New_N+J+1] += u_old[i*N+j]/2.0; //right
            u_new[(I-1)*New_N+J] += u_old[i*N+j]/2.0; //bottom
            u_new[(I+1)*New_N+J] += u_old[i*N+j]/2.0; //top
            u_new[(I-1)*New_N+(J-1)] += u_old[i*N+j]/4.0; //bottom_left
            u_new[(I-1)*New_N+(J+1)] += u_old[i*N+j]/4.0; //bottom_right
            u_new[(I+1)*New_N+(J-1)] += u_old[i*N+j]/4.0; //top_left
            u_new[(I+1)*New_N+(J+1)] += u_old[i*N+j]/4.0; //top_right

        }
    }
}

void u_exact_solution(std::vector<double> &u_exact, double n){
    double x, r , y, theta;
        for(double i = 0.0; i < n; i++){
            y = (2*i/(n-1)) - 1;
            if( y < 0 ){//bottom
                for(double j = 0.0; j < n; j++){
                    x = (2*j/(n-1)) - 1;
                    if(x <= 0 ){ //bottom left
                    r = sqrt(sqrt(pow(x,2.0)+pow(y,2.0)));
                    theta =  0.5 * (atan(y/x) + M_PI);
                    u_exact[i*n + j] = r * sin(theta);
                    }
                    else{ //bottom right
                        r = sqrt(sqrt(pow(x,2.0)+pow(y,2.0)));
                        theta =  0.5 * (atan(y/x) + 2*M_PI);
                        u_exact[i*n + j] = r * sin(theta);
                    }
                }
            }
            if( y >= 0 ){//top
                for(double j = 0.0; j < n; j++){
                    x = (2*j/(n-1)) - 1;
                    if(x < 0 ){//top left
                    r = sqrt(sqrt(pow(x,2.0)+pow(y,2.0)));
                    theta =  0.5 * (atan(y/x) + M_PI);
                    u_exact[i*n + j] = r * sin(theta);
                    }
                    else{//top right
                        r = sqrt(sqrt(pow(x,2.0)+pow(y,2.0)));
                        theta =  0.5 * atan(y/x);
                        u_exact[i*n + j] = r * sin(theta);
                    }
                }
            }
        }

        // IMPLEMENTATION OF THE SLIT IN BC

        int temp;
        temp = (n*n - 1)/2;
        u_exact[temp] = 0.0;
        for(double sl = 1.0; sl<((n+1)/2); sl++){
        y = 0.0;
        x = (2*sl/(n-1));
        r = sqrt(sqrt(pow(x,2.0)+pow(y,2.0)));
        theta =  0.5 * atan2(y,x)*(180/M_PI);
        u_exact[temp + sl] = r * sin(theta);
        }

}

double error_norm (std::vector<double>&u_exact, std::vector<double>&u, double N){
    double err = 0.0;
    double err_norm = 0.0;

    for(int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            err += pow((u_exact[i*N+j]- u[i*N+j]), 2.0);
        }
    }
    err_norm = sqrt(err/(N*N));
    return err_norm;
}


int main(int argc , char *argv[]){

    if (argc <= 1)
    {
        std::cout << "Not enough arguments" << std::endl;
    }
    int l = atoi(argv[1]);
    int N = pow(2, l)+1;
    double h = 2.0/(N-1);

    //declar arrays
    std::vector <std::vector <double> > u(l);
    std::vector <std::vector <double> > res(l);
    std::vector <std::vector <double> > func(l);
    std::vector<double> u_exact(N*N);

    //initialize arrays to zero
    for (int i=l-1; i>=0; i--){
        int local_N = pow (2, i+1)+1;
        u[i] = std::vector <double> (local_N*local_N);
        res[i] = std::vector <double> (local_N*local_N);
        func[i] = std::vector <double> (local_N*local_N);
//        std::cout << i+1 << std::endl;
    }

    //building boundary
    boundary(u[l-1], N);


//    //write the init.dat
//    std::ofstream solution2("init.dat", std::ofstream::out);
//    for (int i=0; i < N; ++i){
//        for(int j=0; j < N; ++j){
//            solution2 << (2*j*h)-1 << "\t";
//            solution2 << (2*i*h)-1 << "\t";
//            solution2 << u[l-1][i*N+j] << std::endl;
//        }
//    //    solution2 << std::endl;
//    }
//    solution2.close();

    //exact solution
    u_exact_solution(u_exact, N);

    //comment out for convergence
//    double L2temp = 0.0;
//    double conver = 0.0;

    std::cout<<"Your Alias: "<<"ADG"<<std::endl;
    struct timeval t0, t;
    gettimeofday(&t0, NULL);
//    struct timeval tim;
//        gettimeofday (&tim, NULL);
//        double t1 = tim.tv_sec + (tim.tv_usec/1000000.0);

    for (int v=0; v<13; v++){
        //complete V-cycle
        int f = l;
        for (int r=1; r<l;r++){

//        std::cout << "F : " <<f << std::endl;
        //smoothing
        for (int i=0; i<2; i++)
            smoothing_rbgs(u[f-1], func[f-1], f);

        //residual, r=f-Au
        residual(u[f-1], func[f-1], res[f-1], f);

        //restric residual
        restriction_residual(res[f-1], func[f-2], f);

        //reducing one level downwards
        f=f-1;

    }

    for (int p=1; p<l; p++){

//        std::cout << "L : " <<l << "P : " << p << std::endl;
        //post-smoothing
        smoothing_rbgs(u[p-1], func[p-1], p);

        //interpolatation
        interpolate(u[p-1], u[p], p);

    }

    //initialize sub-arrays to zero - "U"
    for (int i=l-2; i>=0; i--){
        int local_N = pow (2, i+1)+1;
        for (int j=0; j<local_N*local_N; j++){
            u[i][j] = 0;
//            std::cout << "u[" << i << "][" << j << "] : " << u[i][j] << std::endl;
        }
    }

    //comment out for convergence
//    //find the L2norm after each V-cycle
//        double L2norm = L2norm_v_cycle(u[l-1], func[l-1], res[l-1], l);
//        std::cout << "L2norm of : "<< v+1 << " cycle is : " << L2norm <<std::endl;


//    //convergence after each V-cycle
//        if( v != 0 ){
//        conver = L2norm / L2temp;
//        std::cout << "Convergence rate after : " << v+1 << " cycle(s) : " << conver << std::endl;
//        }
//        L2temp = L2norm;

//        //initialize sub-arrays to zero - "res"
//        for (int i=l-1; i>=0; i--){
//            int local_N = pow (2, i+1)+1;
//            for (int j=0; j<local_N*local_N; j++){
//                res[i][j] = 0;
////                std::cout << "res[" << i << "][" << j << "] : " << res[i][j] << std::endl;
//            }
//        }

//        //error norm
//         std::cout << "Error norm : " << error_norm(u_exact, u[l-1], N) << std::endl;

  }

//    gettimeofday (&tim, NULL);
//    double t2 = tim.tv_sec + (tim.tv_usec/1000000.0);
//    double time = t2 - t1;

    gettimeofday(&t, NULL);
    std::cout << "Wall clock time of MG execution: " << ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t)1000000 + (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3
    << " ms" << std::endl;

//    //print the running time
//    std::cout << "Time to solve the problem : " << time << std::endl;

    //error norm
     std::cout << "Error norm : " << error_norm(u_exact, u[l-1], N) << std::endl;




//Write solution.txt
     std::cout << "Writing solution.txt ......" << std::endl;
//std::ofstream solution("solution.dat", std::ofstream::out);

//// for loop to write the output
//for (int i=0; i < N; ++i){
//    for(int j=0; j < N; ++j){
//        solution << (2*j*h)-1 << "\t";
//        solution << (2*i*h)-1 << "\t";
//        solution << u[l-1][i*N+j] << std::endl;
//    }
////    solution << std::endl;
//}
//solution.close();

////Write residual.txt
//std::ofstream solution1("residual.txt", std::ofstream::out);

//// for loop to write the output
//for (int i=0; i < N; ++i){
//    for(int j=0; j < N; ++j){
//        solution1 << i*h << "\t";
//        solution1 << j*h << "\t";
//        solution1 << res[l-1][i*N+j] << std::endl;
//    }
//    solution1 << std::endl;
//}
//solution1.close();

////Write u_exact.txt
//std::ofstream solution2("u_exact.txt", std::ofstream::out);

//// for loop to write the output
//for (int i=0; i < N; ++i){
//    for(int j=0; j < N; ++j){
//        solution2 << (2*i*h)-1 << "\t";
//        solution2 << (2*j*h)-1 << "\t";
//        solution2 << u_exact[i*N+j] << std::endl;
//    }
////    solution2 << std::endl;
//}
//solution2.close();
}



    // x --> -1->1
    // y --> -1->1

