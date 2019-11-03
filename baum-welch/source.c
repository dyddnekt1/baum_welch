#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <io.h>
#include <conio.h>
#include <math.h>

#include "hmm.h"
#define HEADER_NEXT "hmm1.h"

#define PHONES_NUM 21
#define TEST_MAX_PRO 40
#define MAX_TIME 600
#define CUT 5
// i don't know why but if i don't cut back and forth all value become #INDO
#define MY_K -20

struct _finddata_t fd;

hmmType *pro(char *name);
int pro_index(char *name);
float log_N_dis(float input[][N_DIMENSION], int t, float mean[], float var[]);

hmmType new_phones[PHONES_NUM];

int main()
{	
	int tau = 0,T = 0, t = 0, d = 0, i = 0, j = 0, k = 0;
	for(i = 0; i < PHONES_NUM; i++){
		new_phones[i].name = phones[i].name;
	}

	//baum-welch alg
	FILE* trn_mono = fopen("trn_mono.txt","r");
	char line[_MAX_PATH];
	fgets(line, _MAX_PATH, trn_mono); //drop first line
	printf(line);
	float gamma_sum[PHONES_NUM * N_STATE] = {};
	float gamma_k_sum[PHONES_NUM * N_STATE][N_PDF] = {};
	float xi_sum[PHONES_NUM * N_STATE + 2][PHONES_NUM * N_STATE + 2] = {};
	float sum[PHONES_NUM * N_STATE][N_PDF][N_DIMENSION] = {};
	float minus_sum[PHONES_NUM * N_STATE][N_PDF][N_DIMENSION] = {};
	float squared_sum[PHONES_NUM * N_STATE][N_PDF][N_DIMENSION] = {};
	float temp;
			
	//for each observation sequence
	while(fgets(line, _MAX_PATH, trn_mono)){
		tau++;
		//get test file route
		line[strlen(line) - 1] = '\0';
		if(!strcmp(line,"")) break;
		line[strlen(line) - 4] = '\0';
		char test_route[_MAX_PATH];
		int index = 0;
		for(index = 0; index < strlen(line); index++){
			test_route[index] = line[index + 1];
		}
		strcat(test_route,"txt");
		//make hmm of test file
		float tp[TEST_MAX_PRO * N_STATE + 2][TEST_MAX_PRO * N_STATE + 2] = {};
		stateType *state[TEST_MAX_PRO * N_STATE] = {};
		int Nth_state[TEST_MAX_PRO * N_STATE] = {};

		int test_state_num = 0;
		
		while(1){
			fgets(line, _MAX_PATH, trn_mono);
			line[strlen(line) - 1] = '\0';
			if(!strcmp(line,".")) break;
			else{
				hmmType *add = pro(line);

				if(test_state_num == 0){
					int n = 0, m = 0;
					for(n = 0; n < N_STATE + 2; n++){
						for(m = 0; m < N_STATE + 2; m++){
							tp[test_state_num + n][test_state_num + m] = add->tp[n][m];
						}
					}
					for(n = 0; n < N_STATE ; n++){
						state[test_state_num + n] = &(add->state[n]);
						Nth_state[test_state_num + n] = N_STATE * pro_index(line) + n;
					}
					test_state_num = N_STATE;
				}
				
				else{
					int n = 0, m = 0;
					for(n = 1; n < N_STATE + 2; n++){
						state[test_state_num + n - 1] = &(add->state[n - 1]);
						Nth_state[test_state_num + n - 1] = N_STATE * pro_index(line) + n - 1;
						for(m = 1; m < N_STATE + 2; m++){
							tp[test_state_num + n][test_state_num + m] = add->tp[n][m];
						}
					}
					if(!strcmp(line,"sp")){
						tp[test_state_num][test_state_num + 2] = tp[test_state_num][test_state_num + 1] * add->tp[0][1];
						tp[test_state_num][test_state_num + 1] = tp[test_state_num][test_state_num + 1] * add->tp[0][2];
						test_state_num = test_state_num - (N_STATE - 1);
					}
					test_state_num = test_state_num + N_STATE;
				}
			}
		}
				
		//statistical accumulation
		FILE *test_file = fopen(test_route,"r");
		float input[MAX_TIME][N_DIMENSION] = {};
		
		fscanf(test_file,"%d",&T);
		T = T - 2 * CUT; //cut back and forth
		fscanf(test_file,"%d",&d); // drop 39

		//cut forth
		t = 0;
		for(i = 0; i < CUT; i++){
			for(d = 0; d < N_DIMENSION; d++){
				fscanf(test_file,"%e",&(input[t][d]));
			}
		}
		
		//get input
		for(t = 0; t < T; t++){
			for(d = 0; d < N_DIMENSION; d++){
				fscanf(test_file,"%e",&(input[t][d]));
			}
		}
		
		float b_k[TEST_MAX_PRO * N_STATE][N_PDF][MAX_TIME] = {}; //log
		float b[TEST_MAX_PRO * N_STATE][MAX_TIME] = {}; //log
		float alpha[TEST_MAX_PRO * N_STATE][MAX_TIME] = {}; //log
		float beta[TEST_MAX_PRO * N_STATE][MAX_TIME] = {}; //log
		float l[N_PDF] = {}; 
		float p = 0, temp = 0, max = 0;
		
		//initailize alpha & beta with 1 (cause 1 never happen, it mean that probability is 0)
		for(i = 0; i < TEST_MAX_PRO * N_STATE; i++){
			for(t = 0; t < MAX_TIME; t++){
				alpha[i][t] = 1;
				beta[i][t] = 1;
			}
		}

		//calculate b_k[][][]
		for(j = 0; j < test_state_num; j++){
			for(t = 0; t < T; t++){
				for(k = 0; k < N_PDF; k++){
				b_k[j][k][t] = log_N_dis(input, t, state[j]->pdf[k].mean, state[j]->pdf[k].var);
				}
			}
		}
		
		//calculate b[][]
		for(j = 0; j < test_state_num; j++){
			for(t = 0; t < T; t++){
				max = 0;
				for(k = 0; k < N_PDF; k++){
					if(state[j]->pdf[k].weight == 0) continue;
					l[k] = log(state[j]->pdf[k].weight) + b_k[j][k][t];
					if(max == 0) max = l[k];
					if(max < l[k]) max = l[k];
				}
				temp = 0;
				for(k = 0; k < N_PDF; k++){
					if(state[j]->pdf[k].weight == 0) continue;
					temp = temp + exp(l[k]-max);
				}
				if(temp == 0) continue;
				b[j][t] = max + log(temp);
			}
		}
		
		//forward
		//alpha 1(i)
		for(i = 0; i < test_state_num; i++){
			if(tp[0][i+1] == 0) continue;
			alpha[i][0] = log(tp[0][i+1]) + b[i][0];
		}
		//recursion
		for(t = 0; t < T - 1; t++){
			for(j = 0; j < test_state_num; j++){
				max = 0;
				for(i = 0; i < test_state_num; i++){
					if(alpha[i][t] == 1 || tp[i+1][j+1] == 0) continue;
					l[i] = alpha[i][t] + log(tp[i+1][j+1]);
					if(max == 0) max = l[i];
					if(max < l[i]) max = l[i];
				}
				temp = 0;
				for(i = 0; i < test_state_num; i++){
					if(alpha[i][t] == 1 || tp[i+1][j+1] == 0) continue;
					temp = temp + exp(l[i]-max);
				}
				if(temp == 0) continue;
				alpha[j][t+1] = max + log(temp);
				alpha[j][t+1] = alpha[j][t+1] + b[j][t+1];
			}
		}
		//backward
		//beta T(i) = 1
		for(i = 0; i < test_state_num; i++){
			beta[i][T-1] = log(1.0);
		}
		//recursion
		for(t = T - 2; t >= 0; t--){
			for(i = 0; i < test_state_num; i++){
				max = 0;
				for(j = 0; j < test_state_num; j++){
					if(tp[i+1][j+1] == 0 || beta[j][t+1] == 1) continue;
					l[j] = log(tp[i+1][j+1]) + b[j][t+1] + beta[j][t+1];
					if(max == 0) max = l[j];
					if(max < l[j]) max = l[j];
				}
				temp = 0;
				for(j = 0; j < test_state_num; j++){
					if(tp[i+1][j+1] == 0 || beta[j][t+1] == 1) continue;
					temp = temp + exp(l[j]-max);
				}
				if(temp == 0) continue;
				beta[i][t] = max + log(temp);
			}
		}
		
		//calculate p
		max = 0;
		for(i = 0; i < test_state_num; i++){
			if(alpha[i][T-1] == 1 || tp[i+1][test_state_num + 1] == 0) continue;
			l[i] = alpha[i][T-1] + log(tp[i+1][test_state_num + 1]);
			if(max == 0) max = l[i];
			if(max < l[i]) max = l[i];
		}
		temp = 0;
		for(i = 0; i < test_state_num; i++){
			if(tp[i+1][test_state_num + 1] == 0) continue;
			temp = temp + exp(l[i] - max);
		}
		p = max + log(temp);

		//for each time t
		for(t = 0; t < T; t++){
			//for each state i
			for(i = 0; i < test_state_num; i++){
				//sum gamma
				if(!(alpha[i][t] == 1 || beta[i][t] == 1)){
					gamma_sum[Nth_state[i]] = gamma_sum[Nth_state[i]] + exp(alpha[i][t] + beta[i][t] - p - MY_K);
				}
				//if t = 0
				if(t == 0){
					if(!(alpha[i][0] == 1 || beta[i][0] == 1)){
						xi_sum[0][Nth_state[i] + 1] = xi_sum[0][Nth_state[i] + 1] + exp(alpha[i][0] + beta[i][0] - p - MY_K);
					}
				}
				//if t != T
				if(t != T-1){
					for(j = 0; j < test_state_num; j++){
						if(alpha[i][t] == 1 || tp[i+1][j+1] == 0 || beta[j][t+1] == 1) continue;
						xi_sum[Nth_state[i]+1][Nth_state[j]+1] = xi_sum[Nth_state[i]+1][Nth_state[j]+1] +
						exp(alpha[i][t] + b[j][t+1] + beta[j][t+1] - p - MY_K) * tp[i+1][j+1];
					}
				}
				//for each Gausian k
				for(k = 0; k < N_PDF; k++){
					if(!(alpha[i][t] == 1 || beta[i][t] == 1 || state[i]->pdf[k].weight == 0)){
						temp = alpha[i][t] + beta[i][t] - p + log(state[i]->pdf[k].weight) + b_k[i][k][t] - b[i][t] - MY_K;
						gamma_k_sum[Nth_state[i]][k] = gamma_k_sum[Nth_state[i]][k] + exp(temp);
						for(d = 0; d < N_DIMENSION; d++){
							sum[Nth_state[i]][k][d] = sum[Nth_state[i]][k][d] + exp(temp) * input[t][d];
							squared_sum[Nth_state[i]][k][d] = squared_sum[Nth_state[i]][k][d] + exp(temp) * pow(input[t][d],2);						}
					}
				}
			}
		}
		fclose(test_file);
		printf("%s\n",test_route);
	}
	
	//for each state i
	float M_a[PHONES_NUM * N_STATE + 2][PHONES_NUM * N_STATE + 2];
	for(i = 0; i < PHONES_NUM * N_STATE; i++){
		M_a[0][i] = xi_sum[0][i+1] * exp(MY_K) / tau;
		//for each state j
		for(j = 0; j < PHONES_NUM * N_STATE; j++){
			if(gamma_sum[i] == 0) continue;
			M_a[i+1][j+1] = xi_sum[i+1][j+1] / gamma_sum[i];
		}
		//a[][] = 1 - sigma a[][]
		M_a[i+1][PHONES_NUM * N_STATE + 1] = 1;
		for(j = 0; j < PHONES_NUM * N_STATE; j++){
			M_a[i+1][PHONES_NUM * N_STATE + 1] = M_a[i+1][PHONES_NUM * N_STATE + 1] - M_a[i+1][j+1];
		}
		//for each gaussian
		for(k = 0; k < N_PDF; k++){
			if(gamma_sum[i] != 0) new_phones[i/N_STATE].state[i%N_STATE].pdf[k].weight = gamma_k_sum[i][k] / gamma_sum[i];
			for(d = 0; d < N_DIMENSION; d++){
				if(gamma_k_sum[i][k] == 0) continue;
				new_phones[i/N_STATE].state[i%N_STATE].pdf[k].mean[d] = sum[i][k][d] / gamma_k_sum[i][k];
				new_phones[i/N_STATE].state[i%N_STATE].pdf[k].var[d] = 
				(squared_sum[i][k][d] / gamma_k_sum[i][k]) - pow(new_phones[i/N_STATE].state[i%N_STATE].pdf[k].mean[d],2);
			}
		}
	}
		
	for(k = 0; k < PHONES_NUM; k++){
		new_phones[k].tp[0][1] = 1.0;
		for(i = 0; i < N_STATE; i++){
			for(j = 0; j < N_STATE; j++){
				new_phones[k].tp[i+1][j+1] = M_a[k*3 + i+1][k*3 + j+1];
			}
		}
		new_phones[k].tp[i][N_STATE+1] = 1;
		for(j = 0; j < N_STATE; j++){
			new_phones[k].tp[i][N_STATE+1] = new_phones[k].tp[i][N_STATE+1] - new_phones[k].tp[i][j+1];
		}
	}

	for(k = 0; k < PHONES_NUM; k++){
		//if pronunciation sequence is f-sp-f and sp is skipped it looks like move from f-3 to f-1. in real, it move f-3 to out of f
		if(k == PHONES_NUM - 1) continue; // sil
		new_phones[k].tp[N_STATE][N_STATE+1] = new_phones[k].tp[N_STATE][N_STATE+1] + new_phones[k].tp[N_STATE][1];
		new_phones[k].tp[N_STATE][1] = 0.0;
	}
	
	// sp
	new_phones[17].tp[1][2] = 1 - new_phones[17].tp[1][1];
	new_phones[17].tp[N_STATE][N_STATE+1] = 0;
	
	// sp skip probability. i can't do it...
	new_phones[17].tp[0][1] = phones[17].tp[0][1];
	new_phones[17].tp[0][2] = phones[17].tp[0][2];
	
	fclose(trn_mono);
	
	//make hmm header file
	FILE *header = fopen(HEADER_NEXT,"w");
	fprintf(header,"#define N_STATE		%d\n",N_STATE);
	fprintf(header,"#define N_PDF		%d\n",N_PDF);
	fprintf(header,"#define N_DIMENSION	%d\n",N_DIMENSION);
	fprintf(header,"\n");
	fprintf(header,"typedef struct {\n");
	fprintf(header,"  float weight;\n");
	fprintf(header,"  float mean[N_DIMENSION];\n");
	fprintf(header,"  float var[N_DIMENSION];\n");
	fprintf(header,"} pdfType;\n");
	fprintf(header,"\n");
	fprintf(header,"typedef struct {\n");
	fprintf(header,"  pdfType pdf[N_PDF];\n");
	fprintf(header,"} stateType;\n");
	fprintf(header,"\n");
	fprintf(header,"typedef struct {\n");
	fprintf(header,"  char *name;\n");
	fprintf(header,"  float tp[N_STATE+2][N_STATE+2];\n");
	fprintf(header,"  stateType state[N_STATE];\n");
	fprintf(header,"} hmmType;\n");
	fprintf(header,"\n");

	//print phones
	int a, b, c;
	fprintf(header,"hmmType phones[] = {\n");
	for(a = 0; a < sizeof(new_phones)/sizeof(hmmType); a++){
		//print name
		fprintf(header,"  { \"%s\", // HMM\n",new_phones[a].name);
		
		//print transition probability
		fprintf(header,"    { // transition probability\n");
		for(b = 0; b < N_STATE + 2; b++){
			fprintf(header,"      { ");
			for(c = 0; c < N_STATE + 2; c++){
				fprintf(header,"%e",new_phones[a].tp[b][c]);
				if(c != N_STATE + 1) fprintf(header,", ");
			}
			fprintf(header," }");
			if(b != N_STATE + 1) fprintf(header,",");
			fprintf(header,"\n");
		}
		fprintf(header,"    },\n");

		//print state
		fprintf(header,"    {\n");
		for(b = 0; b < N_STATE; b++){
			fprintf(header,"      {{// state %d\n", b + 1);			
			for(c = 0; c < N_PDF; c++){
				fprintf(header,"        { // pdf %d\n", c + 1);
				
				//print weight
				fprintf(header,"          %e,\n",new_phones[a].state[b].pdf[c].weight);						
				
				//print mean
				fprintf(header,"          { ");
				for(d = 0; d < N_DIMENSION; d++){
					fprintf(header,"%e",new_phones[a].state[b].pdf[c].mean[d]);
					if(d != N_DIMENSION - 1) fprintf(header,",");
					fprintf(header," ");
				}
				fprintf(header,"},\n");

				//print var
				fprintf(header,"          { ");
				for(d = 0; d < N_DIMENSION; d++){
					fprintf(header,"%e",new_phones[a].state[b].pdf[c].var[d]);
					if(d != N_DIMENSION - 1) fprintf(header,",");
					fprintf(header," ");
				}
				fprintf(header,"}\n");
				fprintf(header,"        }");		
				if(c != N_PDF - 1) fprintf(header,",");
				fprintf(header,"\n");
			}
			fprintf(header,"      }}");
			if(b != N_STATE - 1) fprintf(header,",");
			fprintf(header,"\n");	
		}
		fprintf(header,"    }\n");
		fprintf(header,"  }");
		if(a != sizeof(new_phones)/sizeof(hmmType) - 1) fprintf(header,",");
		fprintf(header,"\n");
	}
	fprintf(header,"};\n");
	fclose(header);
	return 0;
}

//========end main==========

hmmType *pro(char *name){
	int i;
	for(i = 0; i < sizeof(phones)/sizeof(hmmType); i++){
		if(!strcmp(name,phones[i].name)){
			return &phones[i];
		}
	}
}
int pro_index(char *name){
	int i;
	for(i = 0; i < sizeof(phones)/sizeof(hmmType); i++){
		if(!strcmp(name,new_phones[i].name)){
			return i;
		}
	}
}

float log_N_dis(float input[][N_DIMENSION], int t, float mean[], float var[]){
	float pro = -log(2.0*M_PI)*N_DIMENSION/2.0;
	int n = 0;
	for(n = 0; n < N_DIMENSION; n++){
		pro = pro - log(var[n])/2.0 - ((1.0/2.0)*(input[t][n]-mean[n])*(input[t][n]-mean[n]))/var[n];
	}
	return pro;
}
