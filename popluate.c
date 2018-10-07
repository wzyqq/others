#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_erf.h>
#include <math.h>
#include <time.h>
#define nhalo         102885216
//#define nhalo         43
#define log_Mmin      12.78
#define sigma_logM    0.68
#define M0        	  pow(10.0,12.71)
#define M1_prime      pow(10.0,13.76)
#define alpha 		  1.15
#define H0            100*1000/(3.08568E22)
#define pi            3.141592653589793
#define G             6.67259E-11
#define delta_c       1.686
#define triangle      200.0
#define eta_num       1000

struct halo_property
{
	double x, y, z; 
	double M;
};
struct galaxy_property
{
	double x,y,z;
	struct galaxy_property *next;
};
int find(double number, double max, double *array)
{
	int loc;
	int i;
	int left = 0, right = max-1;
	int mid;
	while(right>=left)
	{
		mid=(left+right)/2;
		if(array[mid]==number)
		{
			loc = mid;
			break;
		}
		else if(array[mid]>number)
		right=mid-1;
		else
		left=mid+1;
		if(right<=(left+1))
			loc =left;
	}
	return loc;
}
double NFW_random(double random, double mass)
{
	int i;
	int random_loc;
	double random_res;
	double omega_m0 = 0.268;
	double omega_lamba0 = 0.732;
	double H;
	double concentration;
	double eta[eta_num];
	double xi[eta_num];
	double rho,rho_h,r_h,rs;
	double m_nl = 3.3E12;
	rho = 3*H0*H0/(8*pi*G)*omega_m0/1.9891E30*pow((3.08568E22),3);
//	printf("rho=%le\n", rho);
	concentration = 11*pow((mass/m_nl),(-0.13))/1;
//	printf("concentration=%le\n", concentration);
    rho_h = triangle*rho;
//	print('rho_h=',rho_h);
    r_h = pow((3*mass/4/pi/rho_h),(1.0/3));
//	printf("r_h=%lf\n", r_h);
    rs = r_h/concentration;
//	printf("rs=%lf\n",rs);
//	rhos = rh*triangle/3*concentration**3/(log(1+concentration)-concentration/(1+concentration));
    for(i=0; i<eta_num; i++)
    {	
    	eta[i] = concentration*rs/eta_num*(i+0.5);
	    xi[i] = (log(1+eta[i]/rs)-eta[i]/(rs+eta[i]))/(log(1+concentration)-concentration/(1+concentration));
    }
    random_loc = find(random,eta_num,xi);
    random_res = eta[random_loc];
    return random_res;
}
int main()
{
	FILE *fp/*,*fq*/,*fr;
	FILE *TEST1,*TEST2,*TEST3;
	int i,j;
	double temp;
	double N_cen_average, N_sat_average;
	double N_cen,N_sat;
	double N_cen_res,N_sat_res;
	double random1,random2,random3;
	double r_random,theta_random,phi_random;
	double cen_random,sat_random;
	struct halo_property *data = (struct halo_property*)malloc( nhalo*sizeof(struct halo_property));
	struct galaxy_property *head = NULL;
	struct galaxy_property *galaxy = NULL;
	fp = fopen("halo_cat","r");
//	fp = fopen("test.csv","r");
	fr = fopen("/data/s3/zywang/galaxy_cat_zehavi11","w");
	TEST1 = fopen("test1","w");
	TEST2 = fopen("test2","w");
	TEST3 = fopen("test3","w");
	srand((unsigned) time(NULL));
	for (i = 0; i< nhalo; i++)
	{
		fscanf(fp, "%lf%lf%lf%lf\n", &data[i].x, &data[i].y, &data[i].z, &data[i].M);
		data[i].M = pow(10,data[i].M);
//        printf("%lf\t%lf\t%lf\t%lf\n", data[i].x, data[i].y, data[i].z, data[i].M);
		galaxy = (struct galaxy_property*)malloc( sizeof(struct galaxy_property));
		temp = (log10(data[i].M) - log_Mmin) / sigma_logM;
		N_cen_average = (1 + gsl_sf_erf(temp)) * 0.5;
		if((data[i].M - M0)<=0)
		{
			N_sat_average=0;
		}
		else        
		{
			N_sat_average = N_cen_average*pow((data[i].M - M0)/M1_prime,alpha); 		
		}   
//		printf("%le\t%le\t%le\t%le\n", data[i].M, M0, M1_prime, alpha);
//		printf("N_cen_average=%lf\tN_sat_average=%lf\n", N_cen_average, N_sat_average);
		cen_random = rand()/(RAND_MAX+1.0);
//		printf("cen_random=%lf\n", cen_random);	
			N_cen_res = N_cen_average;
//			printf("N_cen_res=%lf\n", N_cen_res);
			if(cen_random < N_cen_res) 
			{
				N_cen = ceil(N_cen_average);
				galaxy->x = data[i].x;
				galaxy->y = data[i].y;
				galaxy->z = data[i].z;
				galaxy->next = head;
				head = galaxy;
				fprintf(fr, "%lf\t%lf\t%lf\t\n", galaxy->x, galaxy->y, galaxy->z);
//				printf("position=%lf\t%lf\t%lf\t\n", galaxy->x, galaxy->y, galaxy->z);
			}
			else N_cen = 0;
//			printf("N_cen=%lf\n", N_cen);
			
			if(N_sat_average>=1)
			{
				for(j=0;j<floor(N_sat_average);j++)
				{
					random1 = rand()/(RAND_MAX+1.0);
					r_random = NFW_random(random1,data[i].M);
					random2 = rand()/(RAND_MAX+1.0);
					phi_random = random2*180.0;
					random3 = rand()/(RAND_MAX+1.0);
					theta_random = asin(random3)*180.0/pi;
					galaxy->x = data[i].x + r_random*cos(theta_random)*cos(phi_random);
					galaxy->y = data[i].y + r_random*cos(theta_random)*sin(phi_random);
					galaxy->z = data[i].z + r_random*sin(theta_random);
					galaxy->next = head;
					head = galaxy;
					fprintf(TEST1, "%lf\n", r_random);
					fprintf(TEST2, "%lf\n", phi_random);
					fprintf(TEST3, "%lf\n", theta_random);
					fprintf(fr, "%lf\t%lf\t%lf\t\n", galaxy->x, galaxy->y, galaxy->z);
//					printf("r_random=%lf\t\n", r_random);
//					printf("%lf\t%lf\t%lf\t\n", galaxy->x, galaxy->y, galaxy->z);
				}
			}
				sat_random = rand()/(RAND_MAX+1.0);
//				printf("sat_random=%lf\n",sat_random);
				N_sat_res = N_sat_average-floor(N_sat_average);
//				printf("N_sat_res=%lf\n", N_sat_res);
				if(sat_random < N_sat_res)
				{
					N_sat = ceil(N_sat_average);
					random1 = rand()/(RAND_MAX+1.0);
					r_random = NFW_random(random1,data[i].M);
					random2 = rand()/(RAND_MAX+1.0);
					phi_random = random2*180.0;
					random3 = rand()/(RAND_MAX+1.0);
					theta_random = asin(random3)*180.0/pi;
					galaxy->x = data[i].x + r_random*cos(theta_random)*cos(phi_random);
					galaxy->y = data[i].y + r_random*cos(theta_random)*sin(phi_random);
					galaxy->z = data[i].z + r_random*sin(theta_random);
					galaxy->next = head;
					head = galaxy;
					fprintf(fr, "%lf\t%lf\t%lf\t\n", galaxy->x, galaxy->y, galaxy->z);
//					printf("%lf\t%lf\t%lf\t\n", galaxy->x, galaxy->y, galaxy->z);
				}
				else  N_sat = floor(N_sat_average);
//				printf("N_sat=%lf\n", N_sat);
	}
//	printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", data[0].row_id, data[0].x, data[0].y, data[0].z, data[0].Mvir, data[0].Mtot, data[0].Rvir, data[0].conc);
//	free(head);
//	free(data);
//	free(galaxy);
//	fclose(fp);
//	fclose(fq);
//	fclose(fr);
	return 0;
}
